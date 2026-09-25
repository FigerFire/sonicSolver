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
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_stateRealization.h"
#include "solver/system/SF_pressureCoupling.h"
#include "core/system/SF_planFragment.h"
#include "app/application/model/SF_configParser.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

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
    config.pressure.coupling.outerCorrectors = 2;
    config.pressure.coupling.pressureCorrectors = 2;
    config.pressure.coupling.nonOrthogonalCorrectors = 1;
    config.pressure.reference.referenceCell = 0;
    config.pressure.reference.referencePressure = 101325.0;
    return config;
}

/// @brief 守恒 transported-rho 单流体系统（没有 div 约束）。
SF::System::BuildRequest conservativeRequest() {
    SF::System::BuildRequest request;
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.singleFluidPreset = SF::System::SingleFluidPresetSpec{false};
    return request;
}

/// @brief rho=rho0 + Momentum + div(U)=0 + p 乘子 的显式组合。
SF::System::BuildRequest incompressibleCompositionRequest() {
    SF::System::BuildRequest request;
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.composition.declared = true;
    request.composition.equations = {"Continuity","Momentum"};
    request.composition.thermoDynamics.equationOfState = "rhoConst";
    request.composition.thermoDynamics.transport = "const";
    request.composition.thermoDynamics.constantDensity = 1.0;
    return request;
}

/// @brief 单流体压力约束方程族（built-in pressure-constraint preset）。
SF::System::BuildRequest pressureConstraintRequest() {
    SF::System::BuildRequest request;
    request.templateOrigin = SF::System::PhysicsTemplateKind::SingleFluid;
    request.pressureConstraint = SF::System::PressureConstraintSpec{false};
    return request;
}

} // namespace

int main() {
    try {
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
        require(SF::System::hasConstraint(
                    manualSystem.executableSystem,"C_INCOMPRESSIBILITY"),
                "manual Momentum+rhoConst composition lost the "
                "incompressibility constraint");
        require(SF::System::hasConstraint(
                    builtinSystem.executableSystem,"C_INCOMPRESSIBILITY"),
                "built-in pressure-constraint preset lost the "
                "incompressibility constraint");
        require(SF::System::hasEquation(
                    manualSystem.executableSystem,"E_MOMENTUM")
                    && SF::System::hasEquation(
                        builtinSystem.executableSystem,"E_MOMENTUM"),
                "momentum equation missing from one incompressible source");

        // 2. coupling preset 贡献 formulation + plan：PISO 必须给出重复层级。
        require(manualSystem.coupling.status
                    == SF::System::CouplingStatus::Active,
                "PISO registered on an incompressible system is not active");
        require(!manualSystem.coupling.derivedOperations.empty(),
                "active coupling preset contributed no derived operations");
        bool hasOuterLoop = false;
        bool hasPressureLoop = false;
        const auto scan = [&](const auto& self,
                              const SF::System::SolvePlanNode& node) -> void {
            if (node.id.find("outerCorrectors") != std::string::npos) {
                hasOuterLoop = hasOuterLoop
                    || node.repetitions
                        == baseConfig().pressure.coupling.outerCorrectors;
            }
            if (node.id.find("pressureCorrectors") != std::string::npos) {
                hasPressureLoop = true;
            }
            for (const auto& child : node.children) self(self,child);
        };
        scan(scan,manualSystem.solvePlan.root);
        require(hasOuterLoop && hasPressureLoop,
                "PIMPLE preset did not contribute its plan loops");

        // 3. 守恒 transported-rho 系统注册 PIMPLE：必须报告 inactive，而不是
        //    触发 density-solver 分支，也不是静默忽略。
        SF::System::BuildRequest conservative = conservativeRequest();
        conservative.coupling = SF::System::couplingRequestFrom(
            baseConfig().pressure.coupling,true);
        conservative.coupling->presetKind = SF::FDM::PressureCouplingPreset::SIMPLE;
        const auto conservativeSystem =
            SF::System::build(baseConfig(),conservative);
        require(conservativeSystem.coupling.status
                    == SF::System::CouplingStatus::Inactive,
                "SIMPLE on a conservative transported-rho system is not "
                "reported inactive");
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
        const auto& presetMomentum = SF::System::equationDefinition(
            conservativeSystem,"E_MOMENTUM");
        const SF::Equation::Definition manualMomentum = SF::Equation::named(
            "E_MOMENTUM",
            SF::Equation::ddt({"rhoU"})
                + SF::Equation::div({"momentumFlux"})
                == SF::Equation::Symbol{"zero"});
        require(presetMomentum.left.terms.size()
                    == manualMomentum.left.terms.size()
                    && presetMomentum.right.terms.size()
                        == manualMomentum.right.terms.size(),
                "preset and manual IR definitions have different structure");

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
        const SF::System::PlanFragment pisoFragment =
            SF::System::couplingPlanFragment(*manual.coupling);
        std::vector<std::string> unresolved;
        require(SF::System::couplingPlanResolved(
                    pisoFragment,manualSystem.executableSystem,&unresolved),
                "PISO fragment references stages the formulation never "
                "declared");
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
        // fixed-point schedule：formulation 不提供 dedicated predictor，
        // fragment 只给出具名缺失标记，planner 不得自行发明 operation。
        SF::System::CouplingPresetRequest simple = *manual.coupling;
        simple.presetKind = SF::FDM::PressureCouplingPreset::SIMPLE;
        const SF::System::PlanFragment simpleFragment =
            SF::System::couplingPlanFragment(simple);
        std::vector<std::string> simpleMissing;
        require(!SF::System::couplingPlanResolved(
                    simpleFragment,manualSystem.executableSystem,
                    &simpleMissing),
                "fixed-point schedule claimed to be fully resolved");
        require(std::find(simpleMissing.begin(),simpleMissing.end(),
                          "pressure.schedule.predictor.solve")
                    != simpleMissing.end(),
                "fixed-point missing marker lost its dedicated predictor id");
        require(std::find(simpleMissing.begin(),simpleMissing.end(),
                          "pressure.schedule.momentum.solve")
                    == simpleMissing.end(),
                "fixed-point schedule silently reuses the physical-time "
                "momentum solve");

        std::cout << "formulation/coupling/realization contract passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[test_formulationArchitecture] FAILED: "
                  << error.what() << '\n';
        return 1;
    }
}
