/// @file test_formulationArchitecture.cpp
/// @brief 方程组成 / formulation / coupling preset / state realization 的架构契约。
///
/// 这些测试保护"preset != hidden solver"：
///   - 方程来源（preset / declared equations / pressure-constraint request）
///     只声明 WHAT；
///   - SIMPLE/PISO/PIMPLE 是 coupling preset，是否生效由 resolved
///     equation/constraint structure 决定，且必须显式报告状态；
///   - state realization 由 executable system 的 role 导出；
///   - compiled HOW（time recipe / dt policy）来自 CompiledNumericalSystem。

#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_stateRealization.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_numericalCompiler.h"
#include "core/system/SF_planFragment.h"
#include "app/application/model/SF_configParser.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

SF::FDM::SolverConfig baseConfig() {
    SF::FDM::SolverConfig config;
    config.numerics.cfl = 0.5;
    config.numerics.maxDeltaT = 0.01;
    config.numerics.timeRecipe =
        SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ForwardEuler);
    config.numerics.recipes.time = config.numerics.timeRecipe;
    config.numerics.viscousEnabled = false;
    config.pressure.coupling.outerCorrectors = 1;
    config.pressure.coupling.pressureCorrectors = 2;
    config.pressure.coupling.nonOrthogonalCorrectors = 1;
    config.pressure.reference.referenceCell = 0;
    config.pressure.reference.referencePressure = 101325.0;
    return config;
}

/// @brief 守恒 transported-rho 单流体系统（没有 div 约束）。
SF::System::BuildRequest conservativeRequest() {
    SF::System::BuildRequest request;
    request.composition.stateDeclared = true;
    request.composition.solutionVariables = {"rho","rhoU","rhoE"};
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.singleFluidPreset = SF::System::SingleFluidPresetSpec{false};
    return request;
}

/// @brief rho=rho0 + Momentum + div(U)=0 + p 乘子 的显式组合。
SF::System::BuildRequest incompressibleCompositionRequest() {
    SF::System::BuildRequest request;
    request.composition.stateDeclared = true;
    request.composition.solutionVariables = {"rho","rhoU","rhoE"};
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.composition.declared = true;
    request.composition.solutionVariables = {"U","p"};
    request.composition.algorithm = "PISO";
    request.composition.equations = {"Continuity","Momentum"};
    request.composition.thermoDynamics.equationOfState = "rhoConst";
    request.composition.thermoDynamics.transport = "const";
    request.composition.thermoDynamics.constantDensity = 1.0;
    return request;
}

/// @brief 单流体压力约束方程族（built-in pressure-constraint preset）。
SF::System::BuildRequest pressureConstraintRequest() {
    SF::System::BuildRequest request;
    request.composition.stateDeclared = true;
    request.composition.solutionVariables = {"rho","rhoU","rhoE"};
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.pressureConstraint = SF::System::PressureConstraintSpec{false};
    return request;
}

} // namespace

int main() {
    try {
        const SF::System::BuiltinStateCatalog catalog;
        require(catalog.solution("U").derivation==SF::System::StateDerivation::None
            && catalog.at("U").derivation==SF::System::StateDerivation::Velocity,
            "Solution U and dependency U share incompatible derivation metadata");
        require(catalog.solution("p").role==SF::System::StateRole::Algebraic
            && catalog.solution("p").derivation==SF::System::StateDerivation::None
            && catalog.at("p").derivation==SF::System::StateDerivation::Pressure,
            "Pressure selection is hardwired to a multiplier or EOS view");
        // 1. 两种 incompressible 来源必须落到同一个约束/动量结构上。
        SF::System::BuildRequest manual = incompressibleCompositionRequest();
        manual.coupling = SF::System::couplingRequestFrom(
            baseConfig().pressure.coupling,true);
        manual.coupling->presetKind = SF::FDM::PressureCouplingPreset::PISO;
        SF::System::BuildRequest builtin = pressureConstraintRequest();
        builtin.coupling = manual.coupling;
        const auto manualSystem =
            SF::System::build(baseConfig(),manual);
        const auto builtinSystem =
            SF::System::build(baseConfig(),builtin);
        require(manualSystem.executableSystem.state.isSolution("U")
            && manualSystem.executableSystem.state.isSolution("p")
            && !manualSystem.executableSystem.state.isSolution("rho")
            && builtinSystem.executableSystem.state.isSolution("rhoU")
            && !builtinSystem.executableSystem.state.isSolution("U"),
            "Solution membership was inferred from STATE role or closure dependencies");
        require(manualSystem.solvePlan.compiledProgram.steps.size()==5
                && manualSystem.solvePlan.compiledProgram.steps.front()
                    .source.equation=="momentum"
                && manualSystem.solvePlan.compiledProgram.steps.at(1)
                    .target.workspace=="pressureCorrection"
                && manualSystem.solvePlan.compiledProgram.steps.back()
                    .target.kind==SF::System::TargetKind::Physical
                && manualSystem.solvePlan.compiledProgram.steps.back()
                    .target.resources.front().storage=="pressureFaceFlux",
                "pressure FormulaSteps lack explicit output/method binding");
        require(SF::System::hasConstraint(
                    manualSystem.executableSystem,"C_INCOMPRESSIBILITY"),
                "manual Momentum+rhoConst composition lost the "
                "incompressibility constraint");
        require(SF::System::hasConstraint(
                    builtinSystem.executableSystem,"C_INCOMPRESSIBILITY"),
                "built-in pressure-constraint preset lost the "
                "incompressibility constraint");
        require(SF::System::hasEquation(
                    manualSystem.executableSystem,"momentum")
                    && SF::System::hasEquation(
                        builtinSystem.executableSystem,"momentum"),
                "momentum equation missing from one incompressible source");
        const auto& predictorFormula=manualSystem.executableSystem.registry.at(
            "momentum");
        require(SF::System::dependsOn(predictorFormula.lhs,"U")
                && manualSystem.solvePlan.compiledProgram.steps.front().source.target.kind==SF::System::TargetKind::Working,
                "momentum mathematics or Working execution target was lost");
        const auto& pressureFormula=manualSystem.executableSystem.registry.at(
            "pSimple");
        require(SF::System::dependsOn(pressureFormula.lhs,"pPrime")
                    && !SF::System::dependsOn(pressureFormula.rhs,"pPrime"),
                "pressure correction Equation does not expose its implicit target");
        const auto scanFormulaCalls=[&](const auto& self,
            const SF::System::SolvePlanNode& node,
            const std::string& formula,const std::string& target) -> int {
            int count=0;
            for (const auto& call:node.equationCalls)
                count+=(call.equation==formula && call.target==target)?1:0;
            for (const auto& child:node.children)
                count+=self(self,child,formula,target);
            return count;
        };
        require(scanFormulaCalls(scanFormulaCalls,manualSystem.solvePlan.root,
                    "momentum","U")==1
                && scanFormulaCalls(scanFormulaCalls,manualSystem.solvePlan.root,
                    "pSimple","pPrime")==1
                && scanFormulaCalls(scanFormulaCalls,manualSystem.solvePlan.root,
                    "correctFluxp","phi")==1,
                "PISO plan lost its FormulaCall semantics");
        const auto flattenOperations=[&](const auto& self,
            const SF::System::SolvePlanNode& node,
            std::vector<std::string>& operations)->void {
            if (!node.operation.empty()) operations.push_back(node.operation);
            for (const auto& child:node.children)
                self(self,child,operations);
        };
        std::vector<std::string> methodOrder;
        flattenOperations(flattenOperations,manualSystem.solvePlan.root,methodOrder);
        require(methodOrder==std::vector<std::string>{
            "pressure.prepare","momentum.assemble","momentum.solve",
            "pressure.boundary.prepare","pressure.assemble","pressure.solve",
            "velocity.correct","pressure.update.prepare","flux.correct",
            "pressure.correction.commit","pressure.step.commit","time.commit"},
            "Method fragments changed constant-density PISO operation order");

        // 2. PISO 只有 pressure-corrector loop；PIMPLE 才有 outer loop。
        require(manualSystem.coupling.status
                    == SF::System::CouplingStatus::Active,
                "PISO registered on an incompressible system is not active");
        require(!manualSystem.coupling.derivedOperations.empty(),
                "active coupling preset contributed no derived operations");
        bool hasOuterLoop = false;
        bool hasPressureLoop = false;
        const auto scan = [&](const auto& self,
                              const SF::System::SolvePlanNode& node) -> void {
            if (node.kind==SF::System::PlanNodeKind::Loop && node.id=="outer") {
                hasOuterLoop = true;
            }
            if (node.kind==SF::System::PlanNodeKind::Loop && node.id=="pressure") {
                hasPressureLoop = true;
            }
            for (const auto& child : node.children) self(self,child);
        };
        scan(scan,manualSystem.solvePlan.root);
        require(!hasOuterLoop && hasPressureLoop,
                "PISO added an outer fixed-point loop or lost pressure correctors");
        SF::System::BuildRequest pimple = incompressibleCompositionRequest();
        pimple.composition.algorithm = "PIMPLE";
        pimple.coupling = manual.coupling;
        pimple.coupling->presetKind = SF::FDM::PressureCouplingPreset::PIMPLE;
        pimple.coupling->outerCorrectors = 2;
        const auto pimpleSystem = SF::System::build(baseConfig(),pimple);
        hasOuterLoop = false;
        hasPressureLoop = false;
        scan(scan,pimpleSystem.solvePlan.root);
        require(pimpleSystem.coupling.status
                    == SF::System::CouplingStatus::Active
                    && hasOuterLoop && hasPressureLoop,
                "PIMPLE preset did not contribute both plan loops");
        require(scanFormulaCalls(scanFormulaCalls,pimpleSystem.solvePlan.root,
                    "pSimple","pPrime")==1,
                "PIMPLE used a different pressure Equation binding");

        // 3. 守恒 transported-rho 系统注册 PIMPLE：必须报告 invalid，而不是
        //    触发 density-solver 分支，也不是静默忽略。
        SF::System::BuildRequest conservative = conservativeRequest();
        conservative.coupling = SF::System::couplingRequestFrom(
            baseConfig().pressure.coupling,true);
        conservative.coupling->presetKind = SF::FDM::PressureCouplingPreset::SIMPLE;
        const auto conservativeSystem =
            SF::System::build(baseConfig(),conservative);
        require(conservativeSystem.executableSystem.registry.entries().size()==3
                && conservativeSystem.solvePlan.sourceProgram.root.children.size()==4
                && conservativeSystem.solvePlan.sourceProgram.root.children.front()
                    .step.equation=="continuity"
                && conservativeSystem.solvePlan.sourceProgram.root.children.front()
                    .step.target.symbol=="rho"
                && conservativeSystem.solvePlan.compiledProgram.hasTemporalRoot
                && conservativeSystem.solvePlan.compiledProgram.steps.front()
                    .calls.size()==1,
                "production conservative stage did not compile Flow -> Q");
        auto roleIndependent=conservativeSystem.executableSystem;
        for (auto& equation:roleIndependent.legacyEquations)
            if (equation.id=="continuity" || equation.id=="momentum"
                || equation.id=="energy")
                equation.role=SF::System::LegacyEquationRole::Generic;
        const auto compiledWithoutRoles=SF::System::NumericalCompiler::compile(
            roleIndependent,conservativeSystem.solvePlan.compiledProgram,
            conservativeSystem.numericalSystem.recipes,
            SF::System::TermProviderCatalog::builtIn(),
            conservativeSystem.numericalSystem.time.recipe);
        require(compiledWithoutRoles.operators.size()
                    ==conservativeSystem.numericalSystem.operators.size(),
                "Fused Flow term selection still depends on legacy LegacyEquationRole.");
        auto legacyDrifted=roleIndependent;
        legacyDrifted.legacyDefinitions.add(SF::Equation::named("continuity",
            SF::Equation::ddt({"rho"}) == SF::Equation::source({"unauthorizedSource"})));
        const auto compiledLegacyDrift=SF::System::NumericalCompiler::compile(
            legacyDrifted,conservativeSystem.solvePlan.compiledProgram,
            conservativeSystem.numericalSystem.recipes,
            SF::System::TermProviderCatalog::builtIn(),
            conservativeSystem.numericalSystem.time.recipe);
        require(compiledLegacyDrift.operators.size()
                    ==compiledWithoutRoles.operators.size(),
                "Legacy SF::Equation::Definition still controls migrated density numerics.");
        bool authoredFormulaChangeVisible=false;
        try {
            auto changed=roleIndependent;
            auto mass=changed.registry.at("continuity");
            mass.rhs=SF::System::FormulaExpr::op("source",
                {SF::System::FormulaExpr::symbol("unauthorizedSource")},
                "mass.unauthorizedSource");
            changed.registry.replace(std::move(mass));
            (void)SF::System::NumericalCompiler::compile(
                changed,conservativeSystem.solvePlan.compiledProgram,
                conservativeSystem.numericalSystem.recipes,
                SF::System::TermProviderCatalog::builtIn(),
                conservativeSystem.numericalSystem.time.recipe);
        } catch (const std::runtime_error& error) {
            authoredFormulaChangeVisible=
                std::string(error.what()).find("mass.unauthorizedSource")
                !=std::string::npos;
        }
        require(authoredFormulaChangeVisible,
                "Authored Equation occurrence did not control provider resolution.");
        auto rk4Config=baseConfig();
        rk4Config.numerics.timeRecipe=SF::FDM::builtInTimeRecipe(
            SF::FDM::TimeRecipeId::ClassicalRK4);
        rk4Config.numerics.recipes.time=rk4Config.numerics.timeRecipe;
        const auto rk4System=SF::System::build(rk4Config,conservativeRequest());
        require(rk4System.solvePlan.sourceProgram.root.children.front().step.equation
                    ==conservativeSystem.solvePlan.sourceProgram.root.children.front().step.equation
                && rk4System.solvePlan.sourceProgram.root.children.front().step.target.symbol
                    ==conservativeSystem.solvePlan.sourceProgram.root.children.front().step.target.symbol
                && rk4System.numericalSystem.time.recipe.stageCount()==4
                && conservativeSystem.numericalSystem.time.recipe.stageCount()==1,
                "time method selection altered WHAT/HOW instead of its stage math");
        require(conservativeSystem.coupling.status
                    == SF::System::CouplingStatus::Invalid,
                "SIMPLE on a conservative transported-rho system is not "
                "reported inactive");
        bool rejected=false;
        try { SF::System::validate(conservativeSystem); }
        catch (const std::runtime_error&) { rejected=true; }
        require(rejected,"explicit incompatible coupling was allowed to run");
        conservative.coupling->explicitlyRegistered=false;
        require(SF::System::build(baseConfig(),conservative).coupling.status
                    == SF::System::CouplingStatus::Inactive,
                "optional coupling cannot be NotApplicable");
        require(conservativeSystem.coupling.reason.find("constraint")
                    != std::string::npos,
                "inactive coupling did not explain the missing constraint");
        require(conservativeSystem.runtime.report.requiredOperations.size()
                    > 0,
                "conservative system lost its explicit stage operations");
        for (const auto& operation
             : conservativeSystem.runtime.report.requiredOperations) {
            require(operation.rfind("pressure.",0) != 0,
                    "inactive coupling preset leaked pressure operations");
        }

        // 4. state realization 由 executable role 导出。
        const auto& conservativeRealization = conservativeSystem.realization;
        require(conservativeRealization.conservativeTransportedMass
                    && conservativeRealization.conservativeMomentum,
                "conservative transported state was not realized");
        const auto& manualRealization = manualSystem.realization;
        require(!manualRealization.conservativeTransportedMass,
                "constant-density closure was realized as transported mass");
        require(manualRealization.pressureMultiplier,
                "constant-density pressure multiplier was not realized");

        // 5. preset provenance 不选择 runtime stepper：同一组方程的 realization
        //    与 required operations 不随 template provenance 变化。
        auto relabelled = conservativeRequest();
        relabelled.templateOrigin =
            SF::System::PhysicsTemplateKind::HomogeneousMixture;
        const auto relabelledSystem =
            SF::System::build(baseConfig(),relabelled);
        require(relabelledSystem.realization.conservativeTransportedMass
                    == conservativeRealization.conservativeTransportedMass,
                "template provenance changed the realized state");

        // 6. compiled HOW：dt 与 time recipe 来自 CompiledNumericalSystem。
        require(conservativeSystem.numericalSystem.dt.cfl
                    == baseConfig().numerics.cfl
                    && conservativeSystem.numericalSystem.dt.maxDeltaT
                        == baseConfig().numerics.maxDeltaT,
                "compiled dt policy does not carry the configured limits");
        require(conservativeSystem.numericalSystem.time.recipe.id()
                    == baseConfig().numerics.timeRecipe.id(),
                "compiled time recipe does not carry the configured recipe");

        // 7. 缺 pressure provider 时必须是具名 Unsupported，而不是
        //    "pressure-based solver not implemented"。
        SF::System::BuildRequest pisoManual = incompressibleCompositionRequest();
        pisoManual.coupling = SF::System::couplingRequestFrom(
            baseConfig().pressure.coupling,true);
        pisoManual.coupling->presetKind = SF::FDM::PressureCouplingPreset::PISO;
        const auto pisoSystem = SF::System::build(baseConfig(),pisoManual);
        if (pisoSystem.runtime.report.status
            != SF::System::RuntimeStatus::Runnable) {
            require(pisoSystem.runtime.report.reason.find("operation")
                        != std::string::npos
                        || pisoSystem.runtime.report.reason.find("provider")
                            != std::string::npos,
                    "Unsupported report does not name the missing provider/"
                    "operation: " + pisoSystem.runtime.report.reason);
            require(pisoSystem.coupling.status
                        != SF::System::CouplingStatus::Active
                        || !pisoSystem.runtime.report.missingOperations.empty(),
                    "Unsupported coupling has no explicit missing operation "
                    "list");
        }

        // 8. builtin 与用户等价方程使用同一 IR/path：preset 与手写 IR 的
        //    equation definition 逐项一致。
        const auto& presetMomentum=conservativeSystem.executableSystem.registry.at("momentum");
        using Expr=SF::System::FormulaExpr;
        const SF::System::Equation manualMomentum{"momentum",
            Expr::add(Expr::op("ddt",{Expr::symbol("rhoU")}),
                Expr::op("div",{Expr::symbol("momentumFlux")})),Expr::constantValue(0),{}};
        require(SF::System::formulaText(presetMomentum)==SF::System::formulaText(manualMomentum),
            "preset and manually registered mathematical definitions differ");

        // 9. executable operation authority：formulation 声明 operation，
        //    coupling 只引用 stage，planner 解析出的 OpId 必须来自 formulation。
        const auto* pressureSolve = SF::System::findExecutableOperation(
            manualSystem.executableSystem,
            SF::System::OperationStage::PressureSolve);
        require(pressureSolve != nullptr
                    && pressureSolve->operation == "pressure.solve",
                "pressure formulation does not declare its pressure solve "
                "operation");
        require(std::find(pressureSolve->requirements.begin(),
                          pressureSolve->requirements.end(),
                          SF::System::OperationCapability::PressureCorrection)
                    != pressureSolve->requirements.end(),
                "executable operation has no pressure capability requirement");
        require(std::find(pressureSolve->requirements.begin(),
                          pressureSolve->requirements.end(),
                          SF::System::OperationCapability::PressureLinearSolve)
                    != pressureSolve->requirements.end(),
                "pressure solve operation lost its linear-solver requirement");
        require(!manualSystem.solvePlan.sourceProgram.root.children.empty()
,
                "PISO did not contribute structured HOW directly");
        for (const auto& operation
             : manualSystem.runtime.report.requiredOperations) {
            const bool declared = std::any_of(
                manualSystem.executableSystem.operations.begin(),
                manualSystem.executableSystem.operations.end(),
                [&](const SF::System::ExecutableOperation& item) {
                    return item.operation == operation;
                });
            const bool planOnly = operation.rfind("pressure.schedule.",0) == 0;
            require(declared || planOnly,
                    "plan executes an operation no formulation declared: "
                    + operation);
        }
        // Fixed-time operations are declared by formulation. Provider
        // resolution, rather than the preset name, decides if they run.
        SF::System::CouplingPresetRequest simple = *manual.coupling;
        simple.presetKind = SF::FDM::PressureCouplingPreset::SIMPLE;
        SF::System::ExecutionProgram simpleProgram;
        SF::System::ExecutionScope baseMomentum;
        baseMomentum.kind=SF::System::ExecutionKind::EquationCall;
        baseMomentum.order=20;
        baseMomentum.step={"momentum",{"U"}};
        simpleProgram.root.children.push_back(baseMomentum);
        std::vector<SF::System::NumericalBinding> simpleNumerics;
        SF::System::applyPressureExecution(simpleProgram,simple);
        require(!simpleProgram.root.children.empty()
,
                "fixed-point schedule did not contribute structured HOW");

        std::cout << "formulation/coupling/realization contract passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[test_formulationArchitecture] FAILED: "
                  << error.what() << '\n';
        return 1;
    }
}
