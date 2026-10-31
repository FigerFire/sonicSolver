#include "SF_scalarMethod.h"
/// @file SF_systemPrinter.cpp
/// @brief 打印 case 最终数学身份，不参与求解和状态修改。

#include "SF_systemPrinter.h"
#include "core/system/SF_operationIds.h"
#include "SF_eulerianAssembly.h"
#include "SF_eulerianTurbulence.h"
#include "SF_immersedMethods.h"
#include "SF_couplingStatus.h"

#include "SF_config.h"
#include "SF_executionContract.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <vector>

namespace SF::System {
namespace {

const char* roleName(StateRole role) {
    switch (role) {
        case StateRole::Primary: return "primary";
        case StateRole::Transported: return "transported";
        case StateRole::Algebraic: return "algebraic";
        case StateRole::Multiplier: return "multiplier";
        case StateRole::Derived: return "derived";
    }
    return "unknown";
}

const char* storageName(StorageBinding binding) {
    switch (binding) {
        case StorageBinding::PackedDistributed: return "packed-distributed";
        case StorageBinding::NamedDistributed: return "named-distributed";
        case StorageBinding::ProviderDistributed: return "provider-distributed";
        case StorageBinding::TransientWorkspace: return "transient-workspace";
        case StorageBinding::SpecializedExecutor: return "specialized-executor";
    }
    return "unknown";
}

bool isPressureConstraint(const std::string& id) {
    return id == "pSimple" || id == "E_SHARED_PRESSURE";
}

bool isAuxiliaryUnknown(const StateSymbol& unknown) {
    return unknown.id == "p" || unknown.id == "deltaP"
        || (unknown.id.size() >= 3
            && unknown.id.compare(unknown.id.size() - 3, 3, "Aux") == 0);
}

bool isConstraintUnknown(const StateSymbol& unknown) {
    return unknown.location == VariableLocation::BodyConstraint
        || unknown.location == VariableLocation::SurfaceConstraint;
}

bool isSolidUnknown(const StateSymbol& unknown) {
    return unknown.location == VariableLocation::SolidGlobal;
}

std::string unknownDescription(const StateSymbol& unknown) {
    std::ostringstream output;
    output << "  " << unknown.id << "  " << unknown.name
           << "  [" << toString(unknown.location) << ", "
           << toString(unknown.ownership) << ", components="
           << unknown.components << ", role=" << roleName(unknown.role)
           << ", storage=" << storageName(unknown.storageBinding);
    if (!unknown.storageKey.empty()) {
        output << ":" << unknown.storageKey;
        if (unknown.componentOffset != 0) {
            output << "@" << unknown.componentOffset;
        }
    }
    output << ", initial=" << unknown.initializationRequired << ", boundary=" << unknown.boundaryRequired
           << ", runtime-storage=" << unknown.runtimeStorageRequired << ", origin=" << toString(unknown.origin.kind) << ":" << unknown.origin.source;
    if (unknown.constantValue) output << ", value=" << *unknown.constantValue;
    output << "]\n";
    return output.str();
}

void printUnknownSection(
        std::ostringstream& output,
        const char* title,
        const ResolvedSimulationSystem& system,
        bool (*predicate)(const StateSymbol&)) {
    output << "\n" << title << "\n";
    bool printed = false;
    for (const auto& unknown : system.executableSystem.state.symbols()) {
        if (predicate(unknown)) {
            output << unknownDescription(unknown);
            printed = true;
        }
    }
    if (!printed) output << "  (none)\n";
}

bool isStateUnknown(const StateSymbol& unknown) {
    return !isAuxiliaryUnknown(unknown) && !isConstraintUnknown(unknown)
        && !isSolidUnknown(unknown);
}

std::string constraintId(const std::string& id) {
    if (isPressureConstraint(id)) return "C_INCOMPRESSIBILITY";
    return id;
}

void printMathematicalCoverage(std::ostringstream& out,const ResolvedSimulationSystem& system) {
    const auto& equations=system.executableSystem;
    const auto& program=system.solvePlan.compiledProgram;
    out << "\nTIME / CONFIGURED INTENT AND ACTUAL EXECUTION\n  configured recipe: "
        << FDM::toString(system.numericalSelection.recipes.time.id()) << "\n";
    if (program.hasTemporalRoot) {
        out << "  actual treatment: compiled recipe stages, count=" << system.numericalSystem.time.recipe.stageCount()
            << "; common physical dt, one terminal clock commit\n";
        for (int i=0;i<system.numericalSystem.time.recipe.stageCount();++i)
            out << "    stage " << i << " abscissa=" << system.numericalSystem.time.recipe.stage(i).abscissa << "\n";
    } else out << "  actual treatment: authored fixed-time HOW with compiled numerical input/output operations; no recipe RK stages\n"
        << "  corrector/outer/pseudo-time loops do not advance the physical clock\n";
    out << "\nMATHEMATICAL EXECUTION COVERAGE\n";
    for (const auto& equation:equations.registry.entries()) {
        out << "  equation " << equation.id << " origin=" << toString(equation.origin.kind)
            << ":" << equation.origin.source << "\n    WHAT: " << formulaText(equation) << "\n"
            << "    mathematical signature: " << Math::text(equation.lhs.signature()) << "\n";
        const auto consumers=equationConsumers(equations,program,equation.id);
        if (consumers.empty()) out << "    BLOCKED: no compiled consumer\n";
        for (const auto* call:consumers) {
            out << "    HOW: " << call->source.occurrence << " equation=" << call->source.equation
                << " target=" << targetText(call->source.target) << "\n    WHICH: " << call->equationMethod
                << " provider=" << call->backendProvider << "\n";
            if (!call->temporalMethod.empty()) out << "      actual treatment: " << call->temporalMethod << "\n";
            for (const auto& read:call->stateUses) out << "      STATE read " << targetText({read.symbol,read.kind})
                << " version=" << toString(read.version) << (read.halo ? " halo-required" : " owner/view") << "\n";
            for (const auto& resource:call->target.resources) out << "      STATE target storage=" << resource.storage
                << " slice=" << resource.componentOffset << "+" << resource.components << "\n";
            for (const auto& need:call->requirements) out << "      requires " << need << "\n";
        }
        out << "    actual compiled leaves (shared fusion/constraint leaves may have several mathematical members):\n";
        for (const auto* node:equationOperations(equations,system.solvePlan,equation.id))
            out << "      " << node->operation << " provider=" << node->provider << " occurrence=" << node->occurrence << "\n";
        for (const auto& change:system.rawSystem.modifications)
            if (change.kind==ModificationKind::Extend && change.targetId==equation.id)
                out << "    added term origin=" << toString(change.origin.kind) << ":" << change.origin.source
                    << " canonical contribution=" << change.detail << "\n";
    }
    for (const auto& constraint:equations.constraints) {
        out << "  constraint " << constraint.id << " origin=" << toString(constraint.origin.kind)
            << ":" << constraint.origin.source << "\n    WHAT: " << constraint.equation
            << "\n    STATE multiplier: " << (constraint.multiplierUnknown.empty()?"none":constraint.multiplierUnknown)
            << "\n    realized by mathematical members:";
        for (const auto& id:constraint.equations) out << " " << id;
        if (constraint.equations.empty()) out << " BLOCKED: no realization references";
        out << "\n";
    }
    for (const auto& closure:equations.closures) {
        out << "  closure: " << closure.description << " origin=" << toString(closure.origin.kind)
            << ":" << closure.origin.source << "\n    equation references:";
        for (const auto& id:closure.equations) out << " " << id;
        out << "\n    STATE references:";
        for (const auto& id:closure.states) out << " " << id;
        out << "\n";
        for (const auto& call:program.steps) {
            bool related=false;
            for (const auto& id:closure.equations) {
                const auto consumers=equationConsumers(equations,program,id);
                related=related || std::find(consumers.begin(),consumers.end(),&call)!=consumers.end();
            }
            for (const auto& id:closure.states)
                related=related || std::any_of(call.stateUses.begin(),call.stateUses.end(),[&](const auto& use){return use.symbol==id;})
                    || std::any_of(call.stateEffects.begin(),call.stateEffects.end(),[&](const auto& effect){return effect.symbol==id;});
            if (related) out << "    producer/consumer occurrence=" << call.source.occurrence
                << " method=" << call.equationMethod << " provider=" << call.backendProvider << "\n";
        }
        if (closure.equations.empty() && closure.states.empty()) out << "    BLOCKED: no realization reference\n";
        out << "    declaration only; no extra PDE, array, loop or refresh implied\n";
    }
    out << "  Boundary/halo/reduction/publication/commit operations are infrastructure or numerical lifecycle; not additional PDEs.\n"
        << "  Exact scopes, repetition limits, termination and operation order follow in COMPILED SOLVE PLAN.\n";
}

void printNumericalSystem(
        std::ostringstream& output,
        const CompiledNumericalSystem& numerical,const CompiledExecutionProgram& program) {
    output << "\nCOMPILED NUMERICAL SYSTEM\n"
           << "  required halo width : " << numerical.requiredHaloWidth << "\n";
    if (numerical.pressureFaceCoupling != PressureFaceCoupling::None)
        output << "  pressure face coupling : RhieChow\n";
    if(numerical.pressureMomentum) {
        const auto& momentum=*numerical.pressureMomentum;
        output << "  compiled pressure momentum : step=" << momentum.occurrence
               << " equation=" << momentum.equation << " target=U time=" << momentum.timeOccurrence
               << " spatial/source occurrences=" << momentum.terms.size()
               << "; kernels consume current predictor U; nu=compiled mu/rhoConst; source returns acceleration\n";
    }
    if (numerical.pressureOperator) {
        output << "  configured pressure treatment (active only with a compiled method/resource consumer) : one physical step; first predictor R(U_n); "
               << "repeated outer predictors use U_n + dt R(U_k) (Legacy fixed-time iteration, "
               << "not strict forwardEuler or a generic implicit recipe)\n";
        output << "  pressure outer convergence : relativeDelta <= "
               << numerical.pressureOperator->relativeTolerance
               << " AND max|div(F_final)| <= "
               << numerical.pressureOperator->absoluteTolerance << "\n";
    }
    if (numerical.operators.empty()) {
        output << "  bound Equation occurrences : (none)\n";
    } else {
        output << "  bound Equation occurrences\n";
        for (const auto& term : numerical.operators) {
            output << "    " << term.formulaId << "@" << term.occurrence
                   << " -> " << term.output << " " << term.operatorName()
                   << "(" << term.primary;
            if (!term.secondary.empty()) output << "," << term.secondary;
            output << ") -> " << (term.recipe
                ? FDM::toString(term.recipe->id()) : "provider-defined")
                   << " provider=" << term.provider
                   << " owner=" << term.providerOwner
                   << " compiledData="
                   << (term.compiledDataAvailable ? "available" : "none")
                   << " side="
                   << (term.side==CompiledSpatialBinding::Side::Right ? "right" : "left");
            if (term.executionEquationId!=term.formulaId)
                output << " execution=" << term.executionEquationId;
            if (term.recipe)
                output << " [" << FDM::toString(term.recipe->temporalRole())
                       << ", halo=" << term.recipe->haloWidth() << "]";
            output << "\n";
            for(const auto& read:term.fieldSource.inputs) {
                output << "      bound source input " << read.symbol << " version=" << toString(read.version) << " (publication checked by STATE reader)\n";
                for(const auto& publisher:program.steps)if(publisher.target.symbol==read.symbol
                    && std::any_cast<CompiledAffineEquation>(&publisher.providerContract))
                    output << "        relation publisher step=" << publisher.source.occurrence << " equation=" << publisher.source.equation
                        << " operation=" << (read.version==StateVersion::Stage || !publisher.stageLocal ? publisher.backendOperation : publisher.temporalPublish) << "\n";
            }
            for(const auto& evaluation:term.fieldSource.evaluation)output << "      source evaluation " << evaluation << "\n";
        }
    }
    output << "  workspace requirements\n";
    if (numerical.workspaceRequirements.empty()) output << "    (none)\n";
    for (const auto& requirement : numerical.workspaceRequirements) {
        output << "    " << requirement << "\n";
    }
    output << "  provider requirements\n";
    if (numerical.providerRequirements.empty()) output << "    (none)\n";
    for (const auto& requirement : numerical.providerRequirements) {
        output << "    " << requirement << "\n";
    }
}

void printAlgebraicConstraints(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nALGEBRAIC CONSTRAINTS\n";
    bool printed = false;
    std::vector<std::string> emitted;
    const auto emit = [&](const std::string& id, const std::string& name) {
        const std::string displayId = constraintId(id);
        if (std::find(emitted.begin(), emitted.end(), displayId)
            != emitted.end()) return;
        emitted.push_back(displayId);
        output << "  " << displayId << "  " << name << "\n";
        printed = true;
    };
    for (const auto& constraint : system.executableSystem.constraints) {
        emit(constraint.id, constraint.name);
    }
    if (!printed) output << "  (none)\n";
}

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

bool hasConstraintVariables(const ResolvedSimulationSystem& system) {
    return std::any_of(
        system.executableSystem.state.symbols().begin(),
        system.executableSystem.state.symbols().end(),
        [](const StateSymbol& unknown) {
            return unknown.location == VariableLocation::BodyConstraint
                || unknown.location == VariableLocation::SurfaceConstraint;
        });
}

std::string requirementStatus(
        const ExecutionRequirement& requirement,
        const ResolvedSimulationSystem& system) {
    if (requirement.required) return requirement.available ? "OK" : "MISSING";
    if (requirement.name == "ConstraintGlobalDof"
        && hasConstraintVariables(system)) return "serial-only";
    return "not-required";
}

void printContributions(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nCONTRIBUTIONS\n";
    if (system.rawSystem.contributions.empty()) {
        output << "  (none)\n";
        return;
    }
    for (const auto& contribution : system.rawSystem.contributions) {
        output << "  " << toString(contribution.origin.kind) << ": "
               << contribution.id << " module=" << contribution.module;
        if (!contribution.name.empty()) output << "  " << contribution.name;
        output << "\n";
    }
    for (const auto& symbol:system.rawSystem.state.symbols())
        output << "  STATE " << symbol.id << " <- " << symbol.origin.source << "\n";
    for (const auto& closure:system.rawSystem.closures) if (!closure.outputs.empty()) {
        output << "  closure outputs:";
        for (const auto& id:closure.outputs) output << " " << id;
        output << " <- " << closure.origin.source << "\n";
    }
    for (const auto& modification : system.rawSystem.modifications) {
        output << "  modification: " << toString(modification.kind) << " "
               << modification.targetKind << " " << modification.targetId
               << "  origin=" << toString(modification.origin.kind) << ":" << modification.origin.source << "\n";
    }
}

void printRawSystem(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    if (system.executableSystem.immersed) {
        const auto& d=*system.executableSystem.immersed;
        output << "\nFROZEN IBM SELECTION\n  implementation=" << (d.enforcement==FDM::IBMEnforcement::GhostCell ? "ghostCellIBM" : FDM::toString(d.algorithm))
            << " enforcement=" << FDM::toString(d.enforcement) << " support=" << FDM::toString(d.support)
            << " normalization=" << FDM::toString(d.surfaceNormalization)
            << "\n  mathematical semantics: " << d.variational.stationaryFunctional
            << "\n  fluid port: rho / U -> rhoU / rhoE; conservative five-variable implementation only"
            << "\n  provider support and AST are validated independently of this selection.\n";
    }
    output << "\nRAW EQUATION SYSTEM\n  constraints:\n";
    if (system.rawSystem.constraints.empty()) output << "    (none)\n";
    for (const auto& constraint : system.rawSystem.constraints) {
        output << "    " << constraint.id << "  " << constraint.name << "\n";
    }
}

void printTransformations(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nSYSTEM TRANSFORMATIONS\n";
    if (system.transformations.empty()) output << "  (none)\n";
    for (const auto& transformation : system.transformations) {
        output << "  " << transformation.descriptor.id
               << "  " << transformation.descriptor.name
               << "  status=" << toString(transformation.state) << "\n"
               << "    reason: " << transformation.reason << "\n";
        if (!transformation.generatedEquations.empty()) {
            output << "    generated equations:";
            for (const auto& id : transformation.generatedEquations) {
                output << " " << id;
            }
            output << "\n";
        }
        if (!transformation.generatedOperators.empty()) {
            output << "    generated operators:";
            for (const auto& id : transformation.generatedOperators) {
                output << " " << id;
            }
            output << "\n";
        }
    }
}

void printCompiledBindings(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nCOMPILED EQUATION BINDINGS\n";
    if (system.executableSystem.compiledEquations.empty()) {
        output << "  (none)\n";
        return;
    }
    for (const auto& equation : system.executableSystem.compiledEquations) {
        output << "  " << equation.equationId
               << "  operator=" << equation.operatorBinding
               << "  matrix=" << (equation.assemblesMatrix ? "yes" : "no")
               << "  rhs=" << (equation.assemblesRhs ? "yes" : "no") << "\n";
        for (const auto& resource : equation.resources) {
            output << "    " << resource.symbol << " -> " << resource.storage
                   << "[" << resource.componentOffset << ":"
                   << resource.components << "]"
                   << "  access=" << toString(resource.access)
                   << "  boundary="
                   << (resource.boundaryFreshnessRequired ? "fresh" : "not-required")
                   << "  sync=" << toString(resource.synchronization) << "\n";
        }
    }
}

void printPlanNode(
        std::ostringstream& output,
        const SolvePlanNode& node,
        int depth) {
    output << std::string(static_cast<std::size_t>(depth*2),' ')
           << toString(node.kind) << "  " << node.id;
    if (!node.name.empty()) output << "  " << node.name;
    if (node.repetitions != 1) output << "  repeat=" << node.repetitions;
    if (!node.terminationSignal.empty()) {
        output << "  maxIterations=" << node.repetitions
               << "  termination=" << node.terminationSignal;
    }
    if(!node.residuals.empty()) {
        output << " res=" << (node.residualAll?"AND":"OR") << " min=" << node.minimumIterations
               << " max=" << node.repetitions << " onMaximum=" << (node.failOnMaximum?"fail":"continue-explicit");
        for(const auto& r:node.residuals)output << " [" << r.id << "@" << r.occurrence << " kind=" << toString(r.kind)
            << " norm=" << toString(r.norm) << " normalization=" << toString(r.normalization)
            << " atol=" << r.absoluteTolerance << " rtol=" << r.relativeTolerance << " scale=" << r.scale << "]";
    }
    if (!node.operation.empty()) {
        output << "  operation=" << node.operation << " provider=" << node.provider;
        if(node.operation==OpIds::TimeCommit)output << " role=physical-clock-advance";
        else if(node.kind==PlanNodeKind::Commit)output << " role=state/result-publication(no-clock)";
    }
    output << "\n";
    for (const CompiledMathRef& call:node.equationCalls) {
        output << std::string(static_cast<std::size_t>((depth+1)*2),' ')
               << call.equation << " -> " << call.target << "\n";
    }
    for (const auto& child : node.children) {
        printPlanNode(output,child,depth+1);
    }
}

void printCompiledPlan(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nCOMPILED SOLVE PLAN\n";
    if (system.runtime.report.status!=RuntimeStatus::Runnable)
        output << "  NOT EXECUTABLE: declared structure only; see blocking runtime requirements\n";
    printPlanNode(output,system.solvePlan.root,1);
    output << "\nRUNTIME STATUS\n"
           << "  " << (system.runtime.report.status == RuntimeStatus::Runnable
                ? "runnable" : system.runtime.report.status == RuntimeStatus::Invalid
                    ? "invalid" : "unsupported") << "\n";
    output << "\nREQUIRED OPERATIONS\n";
    if (system.runtime.report.requiredOperations.empty()) output << "  (none)\n";
    for (const auto& operation : system.runtime.report.requiredOperations) {
        output << "  " << operation << "\n";
    }
    output << "\nMISSING OPERATION PROVIDERS\n";
    if (system.runtime.report.missingOperations.empty()) output << "  (none)\n";
    for (const auto& operation : system.runtime.report.missingOperations) {
        output << "  " << operation << "\n";
    }
    if (!system.runtime.report.reason.empty()) {
        output << "\nRUNTIME CAPABILITY REASON\n"
               << "  " << system.runtime.report.reason << "\n";
    }
}

} // namespace

std::string describe(const ResolvedSimulationSystem& system) {
    std::ostringstream output;
    output << "\n============================================================\n"
           << "sonicSolver - Resolved Mathematical System\n"
           << "============================================================\n"
           ;
    if (system.numericalSystem.thermophysical) {
        const auto& t=*system.numericalSystem.thermophysical;
        output << "\nTHERMOPHYSICAL / FROZEN MODEL BINDING\n  source: " << t.source
               << "\n  equationOfState: " << t.selection.equationOfState
               << "\n  thermo: " << (t.selection.thermo.empty()?"not selected":t.selection.thermo)
               << "\n  transport: " << (t.selection.transport.empty()?"not selected":t.selection.transport) << "\n";
        if (t.selection.equationOfState=="perfectGas")
            output << "  gamma=" << t.gamma << " R=" << t.gasConstant << " cp=" << t.cp()
                   << " cv=" << t.cv() << " k=" << t.conductivity() << "\n";
        else output << "  rho=" << t.selection.constantDensity << " (no PerfectGas caloric closure)\n";
        output << "  mu=" << t.dynamicViscosity << " Pr=" << t.prandtl << "\n";
        for (const auto& [name,source]:t.provenance) output << "  " << name << " origin: " << source << "\n";
    }
    if(system.numericalSystem.conservative) {
        const auto& c=*system.numericalSystem.conservative;
        output<<"\nCONSERVATIVE EXECUTION INPUTS\n  static origin: "<<c.origin
            <<"\n  IBM stencil policy="<<static_cast<int>(c.ibmBoundary)<<" ILW order="<<c.ilwOrder
            <<" interface flux policy="<<static_cast<int>(c.interfaceFlux)
            <<"\n  molecular compatibility mu="<<c.molecularViscosity<<" Pr="<<c.prandtl
            <<" gamma="<<c.caloricGamma<<" R="<<c.gasConstant
            <<"\n  actual convection gamma/EOS: bound FluidStateModel capability; no runtime config selection"
            <<"\n  actual T/mu/conductivity: bound Field thermodynamic state, plus selected transport resource"
            <<"\n  physical BC: bound boundary pipeline/applicator; pressure-corrector settings remain a separate legacy input\n";
    }
    printMathematicalCoverage(output,system);
    output << "\nWHAT / EQUATIONS\n";
    for (const Equation& formula:system.executableSystem.registry.entries()) {
        output << "  " << formula.id << ": " << formulaText(formula) << " [" << formula.origin.source << "]\n";
    }
    output << "\nBOUNDARY / STENCIL CONTRACTS\n";
    for (const auto& boundary:system.executableSystem.boundaryClosures) {
        output << "  " << boundary.id << " -> " << boundary.provider
               << " [spatial stage time; no EquationCall]\n    reads:";
        for (const auto& read:boundary.reads) output << " " << read;
        output << "\n    writes:";
        for (const auto& write:boundary.writes) output << " " << write;
        output << "\n    order:";
        for (const auto& step:boundary.order) output << " -> " << step;
        for (const auto& item:boundary.boundCapabilities) output << "\n    bound capability: " << item.name
            << " equation=" << item.equation << " fluid-port=" << item.fluidPort << " boundary=" << item.boundary;
        output << "\n";
    }
    output << "\nIMMERSED / COMPILED BLOCK CONTRACTS\n";
    for (const auto& call:system.solvePlan.compiledProgram.steps) {
        const auto* contract=std::any_cast<std::shared_ptr<const CompiledImmersedContract>>(&call.providerContract);
        if (!contract || !*contract) continue;
        const auto& c=**contract;
        output << "  " << call.source.occurrence << ": " << FDM::toString(c.algorithm)
               << " / " << FDM::toString(c.enforcement) << " -> " << c.equation << " / " << c.target
               << " [" << call.backendProvider << "; complete predictor -> correction -> commit; targetTime=time+dt]\n";
        if (c.support==FDM::IBMConstraintSupport::Surface)
            output << "    surface transfer: " << FDM::toString(c.surfaceNormalization) << "; adjoint spreading\n";
        if (c.solid==FDM::IBMSolidModel::SelfPropelledRigid && c.enforcement!=FDM::IBMEnforcement::MonolithicKKT)
            output << "    virtual-fluid generalized mass; active rigid DOFs: "
                   << (c.rigidMotionMode==FDM::IBMRigidMotionMode::Rotate?"rotation":"translation") << "\n";
    }
    output << "\nSTATE — selected solution variables\n  origin : "
           << system.executableSystem.state.selectionOrigin() << "\n  use :";
    for (const auto& id:system.executableSystem.state.solutionVariables()) output << " " << id;
    output << "\nSTATE — model contributed transported variables\n  use :";
    for (const auto& symbol:system.executableSystem.state.symbols())
        if (symbol.role==StateRole::Transported
            && !system.executableSystem.state.isSolution(symbol.id)) output << " " << symbol.id;
    output << "\nSTATE — active dependencies\n  require :";
    for (const auto& symbol:system.executableSystem.state.symbols())
        if (!system.executableSystem.state.isSolution(symbol.id)) output << " " << symbol.id;
    output << "\n";
    output << "\nSTATE / BASE VARIABLES (active case symbols only)\n";
    for (const auto& symbol:system.executableSystem.state.symbols())
        output << unknownDescription(symbol);
    for (const auto& symbol:system.executableSystem.state.symbols()) if (!symbol.dependencies.empty()) {
        output << "    dependency " << symbol.id << (symbol.evaluation==StateEvaluation::Lazy ? " [lazy/versioned view]" : " [materialized]") << " <-";
        for (const auto& dependency:symbol.dependencies) output << " " << dependency;
        output << "\n";
    }
    for (const auto& placement:system.solvePlan.sourceProgram.requirements) {
        output << "    HOW placement: " << placement.occurrence << " scope=" << (placement.scope.empty() ? "root" : placement.scope);
        for (const auto& after:placement.after) output << " after=" << after;
        for (const auto& before:placement.before) output << " before=" << before;
        output << "\n";
    }
    output << "\nSTATE VIEWS / COMPILED DEMANDS\n"
           << "  requested views (storage authority; no execution order):\n";
    for (const auto& view:system.solvePlan.compiledProgram.stateViews) {
        output << "    " << toString(view.kind) << "(" << view.symbol;
        if (view.stage>=0) output << "," << view.stage;
        output << ") -> " << view.storage << " [" << view.components << "] "
               << (view.owner==StateViewOwner::PhysicalState ? "physical-state"
                   : view.owner==StateViewOwner::CompilerWorkspace ? "compiler-workspace" : "numerical-provider")
               << "\n";
    }
    for (const auto& symbol:system.executableSystem.state.symbols()) if (symbol.producer) {
        output << "    derived producer: " << symbol.id << " <- " << symbol.producer->id
               << " / " << (symbol.evaluation==StateEvaluation::Lazy ? "original versioned lazy closure" : "explicit HOW Update") << "\n";
        output << "      relation: " << formulaText(Equation{"derived."+symbol.id,symbol.expression(),symbol.producer->expression}) << "\n";
        for (const auto& input:symbol.dependencies) output << "      input: " << input << " [same requested view version]\n";
        output << "      supported views:";
        for (auto kind:symbol.producer->versions) output << " " << toString(kind);
        output << "\n";
    }
    output << "\nSTATE / CALL INPUTS (single authored version authority)\n";
    for(const auto& binding:system.executableSystem.state.callInputs()) {
        output<<"  "<<(binding.reference==ReferenceKind::Step?"Step ":"Equation ")<<binding.call<<" ["<<binding.origin.source<<"]";
        for(const auto& read:binding.reads)output<<" "<<targetText({read.symbol,read.kind})<<":"<<toString(read.version);
        output<<"\n";
    }
    for(const auto& residual:system.solvePlan.compiledProgram.residuals)output<<"  res("<<residual.id<<") occurrence="<<residual.occurrence<<" target="<<residual.target<<" metric="<<toString(residual.kind)<<" norm="<<toString(residual.norm)<<" normalization="<<toString(residual.normalization)<<" provider="<<residual.provider<<"; checked after full body; reference reset at loop entry\n";
    output << "\nHOW / EXECUTION\n";
    const auto printHow=[&](const auto& self,const ExecutionScope& node,int depth,ExecutionFrequency inherited,bool iteration)->void {
        auto frequency=node.frequency==ExecutionFrequency::Inherit?inherited:node.frequency;
        if(node.frequency==ExecutionFrequency::Inherit && node.kind==ExecutionKind::EquationCall) {
            const auto call=std::find_if(system.solvePlan.compiledProgram.steps.begin(),system.solvePlan.compiledProgram.steps.end(),[&](const auto& value){return value.source.occurrence==node.step.occurrence;});
            if(call!=system.solvePlan.compiledProgram.steps.end() && (call->temporalResidual || call->stageLocal))frequency=ExecutionFrequency::Stage;
        }
        if(node.kind==ExecutionKind::StageLoop)frequency=ExecutionFrequency::Stage;
        iteration=iteration || node.kind==ExecutionKind::Loop;
        output << std::string(static_cast<std::size_t>(depth)*2,' ')
               << node.order << " " << toString(node.kind);
        if (node.kind==ExecutionKind::EquationCall) {
            output << " " << node.step.equation;
            if (!node.step.occurrence.empty()) output << "@" << node.step.occurrence;
            output << " -> " << targetText(node.step.target);
        } else if (!node.id.empty()) {
            output << " " << node.id;
        }
        if(node.kind==ExecutionKind::BlockSolve){output<<" [";for(const auto& e:node.equations)output<<" "<<e;output<<" ] => [";for(const auto& u:node.unknowns)output<<" "<<u;output<<" ]";}
        if (node.kind==ExecutionKind::Loop)
            output << " limit=" << node.repetitions;
        if (!node.terminationSignal.empty())
            output << " until=" << node.terminationSignal;
        if(!node.residuals.empty()){output<<" res="<<(node.residualAll?"AND":"OR");for(const auto& id:node.residuals)output<<" "<<id;output<<" minimum="<<node.minimumIterations<<" onMaximum="<<(node.failOnMaximum?"fail":"continue-explicit");}
        for(std::size_t i=0;i<node.commitSources.size();++i)output<<" COPY "<<targetText(node.commitSources[i])<<" -> "<<targetText(node.commitTargets[i])<<" (clock unchanged)";
        output << " frequency=" << (iteration?"Iteration within ":"") << toString(frequency) << " [" << node.origin.source << "]\n";
        for (const auto& child:node.children) self(self,child,depth+1,frequency,iteration);
    };
    printHow(printHow,system.solvePlan.compiledProgram.root,1,ExecutionFrequency::Step,false);
    if (system.solvePlan.compiledProgram.steps.empty())
        output << "  Legacy numerical backend / migration pending\n";
    output << "\nWHICH / NUMERICS\n";
    output << "  time: " << FDM::toString(system.numericalSelection.recipes.time.id()) << "\n";
    for (const auto& binding:system.numericalSelection.bindings) {
        for (const auto& parameter:binding.parameters)
            output << "    " << binding.equation << "." << parameter.first << "=" << parameter.second << "\n";
        output << "  " << binding.equation;
        if (!binding.occurrence.empty()) output << "@" << binding.occurrence;
        output << " -> " << binding.method << " origin=" << toString(binding.origin.kind) << ":" << binding.origin.source << "\n";
    }
    if (const auto* pressure=requiredPressureInputs(programOf(system))) {
        output << "  required pressure storage: U[3]/velocity owner=Field; p[1]/pressure owner=application; rho=rhoConst\n"
               << "    Working/Correction, pressureFaceFlux/canonical COPY and matrix/linear solver owner=PressureOperators\n"
               << "    one orthogonal source zone, one patch per rank; fixed U/p/rho/phi/deltaP symbols; no implied phase/KKT or arbitrary scalar composition\n"
               << "    reference cell=" << pressure->reference.referenceCell << " value=" << pressure->reference.referencePressure
               << " (runtime gauge; application initial fields are applied separately)\n"
               << "    mu=" << pressure->dynamicViscosity << " relaxation(U,p)=" << pressure->velocityRelaxation << "," << pressure->pressureRelaxation
               << " boundary entries(U,p)=" << pressure->velocityBoundary.size() << "," << pressure->pressureBoundary.size() << "\n"
               << "  linear: HYPRE (frozen pressure linear configuration); face coupling: RhieChow\n"
               << "    resource requirements are not proof of successful runtime binding\n";
    }
    if(const auto* phase=requiredPhaseInputs(programOf(system))) {
        output << "  required phase resources: PhaseSystem primary per-phase mass/momentum/enthalpy, N-1 alpha closure and shared pressure\n"
            << "    application owns phase/turbulence/source/boundary objects; EulerianStepper borrows them\n"
            << "    old-time/outer snapshots, pressure correction, paired canonical phase flux and linear systems owner=Eulerian numerical workspace\n"
            << "    relax(momentum,p)=" << phase->momentumRelaxation << "," << phase->pressureRelaxation
            << " reference(cell,value)=" << phase->reference.referenceCell << "," << phase->reference.referencePressure
            << " pressure boundary entries=" << phase->pressureBoundary.size() << "\n"
            << "    initialPressure must still equal referencePressure; distinct initial and correction gauge semantics\n"
            << "    complete fixed mathematical contracts; one local patch; no implied mixed phase/packed/constant-pressure assembly\n"
            << "    resource requirements are not proof of successful runtime binding\n";
    }
    output << "  frozen method -> runtime provider bindings:\n";
    for (const auto& call:system.solvePlan.compiledProgram.steps)
        output << "    " << call.source.equation << "@" << call.source.occurrence
               << " method=" << call.equationMethod << " provider=" << call.backendProvider << "\n";
    const auto packed=requiredPackedComponents(system.executableSystem.state,system.solvePlan.compiledProgram,"conservative");
    if(packed)output << "  requested backing: conservative components=" << packed
        << " requiredOwner=Field borrowed by STATE; numerical workspace allocated during operation binding\n"
        << "    operations bind by required OpId; field ports by compiled physical resources; binding success is a runtime check; phase assembly remains transitional\n";
    output << "  compiled physical resource requirements (binding validated before execution):\n";
    for (const auto& call:system.solvePlan.compiledProgram.steps)
        for (const auto& resource:call.physicalResources) {
            output << "    " << resource.id << " owner=" << resource.owner
                << " target=" << (resource.fields?resource.fields->target:call.target.symbol) << " access=" << (resource.writable?"write":"read");
            if(resource.fields) {
                output<<" states=[";for(const auto& id:resource.fields->states)output<<" "<<id;
                output<<" ] boundaries=[";for(const auto& id:resource.fields->boundaries)output<<" "<<id;
                output<<" ] sources=[";for(const auto& id:resource.fields->sources)output<<" "<<id;
                output<<" ] boundaryRequired="<<resource.fields->boundaryRequired;
            }
            output<<"\n";
        }
    output << "  explicit provider fusion contracts:\n";
    std::vector<std::string> shownFusion;
    for (const auto& call:system.solvePlan.compiledProgram.steps) {
        if (call.fusionKey.empty() || std::find(shownFusion.begin(),shownFusion.end(),call.fusionKey)!=shownFusion.end()) continue;
        shownFusion.push_back(call.fusionKey);
        output << "    " << call.fusionKey << " {";
        for (const auto& member:call.fusionMembers) output << " " << member.equation << "->" << member.target;
        output << " } -> " << call.backendOperation << " provider=" << call.backendProvider
               << " temporal=" << call.temporalMethod << "\n";
    }
    if(!system.solvePlan.compiledProgram.temporalParticipants.empty()) {
        output << "\nSYNCHRONOUS TEMPORAL GROUP\n"
               << "  one dt / recipe / physical clock / terminal Commit\n"
               << "  snapshots -> all prepare -> coefficient publication -> ordered RHS -> all advance -> validate all -> publish -> final coefficients -> time.commit\n";
        for(const auto& p:system.solvePlan.compiledProgram.temporalParticipants) {
            output<<"  participant "<<p.identity<<" provider="<<p.provider<<"\n    calls:";
            for(const auto& c:p.calls)output<<" "<<c.equation<<"->"<<c.target;
            for(const auto& r:p.stageRequires)output<<"\n    same-stage requires: "<<r;
            for(const auto& r:p.stageProvides)output<<"\n    RHS provides: "<<r;
            if(!p.coefficientPrepare.empty())output<<"\n    coefficients: "<<p.coefficientPrepare<<" / "<<p.coefficientPublish;
            for(const auto& input:p.coefficientInputs)output<<"\n    coefficient preparation input: "<<input.symbol<<" version="<<toString(input.version)<<" after Ready / before authored RHS";
            output<<"\n    RHS evaluations per stage: 1 (complete ordered numerical-owner group)";
            if(!p.coefficientProvides.empty()) {
                output<<"\n    coefficient publisher: "<<toString(p.coefficientVersion)<<" resource="<<p.coefficientProvides.front()<<" states:";
                for(const auto& symbol:p.coefficientStates)output<<" "<<symbol;
            }
            for(const auto& use:p.coefficientUses)output<<"\n    coefficient read: "<<use.symbol<<" version="<<toString(use.version)<<" producer="<<use.producer<<" resource="<<use.resource;
            output<<"\n    phases: "<<p.snapshot<<" / "<<p.prepareStage<<" / "<<p.rhs<<" / "<<p.advance<<" / "<<p.publish<<"\n";
        }
    }
    output << "\nCOMPILED OCCURRENCES / STORAGE / PROVIDERS\n";
    for (const auto& step:system.solvePlan.compiledProgram.steps) {
        if(const auto* block=std::any_cast<CompiledLinearBlock>(&step.providerContract)){output<<"  [";for(const auto& e:block->equations)output<<" "<<e;output<<" ] => [";for(const auto& u:block->targets)output<<" "<<u;output<<" ]\n";}else output<<"  "<<step.source.equation<<" -> "<<targetText(step.source.target)<<"\n";
        output << "    equation method: " << step.equationMethod << "\n"
               << "    selected runtime provider: " << step.backendProvider << "\n";
        for(const auto& origin:step.numericalOrigins)
            output<<"    WHICH selection layer: "<<toString(origin.kind)<<":"<<origin.source<<"\n";
        for(const auto& op:step.numericalOperators)
            output<<"    effective operator: "<<op.formula<<"/"<<op.occurrence<<" -> "<<op.provider<<"\n";
        for(const auto& use:step.coefficientUses)output<<"    coefficient binding: "<<use.symbol<<" version="<<toString(use.version)<<" producer="<<use.producer<<" resource="<<use.resource<<"\n";
        if (const auto* contract=eulerianAssemblyContract(step)) {
            output << "    compiled assembly contract: " << toString(contract->relation);
            if (contract->phaseSlot) output << " phaseSlot=" << *contract->phaseSlot << " phase=" << contract->phaseName;
            output << "\n    supported source extensions:";
            for (auto extension:contract->extensions) output << " " << toString(extension);
            output << "\n";
        }
        if (const auto* contract=eulerianTurbulenceContract(step)) {
            output<<"    compiled turbulence contract: "<<FDM::toString(contract->model)
                <<(contract->transport ? " phase-mass weighted transport" : " algebraic closure")<<"\n    selected phases:";
            for (const auto& phase:contract->phases) output<<" "<<phase.name<<"(slot="<<phase.slot<<")";
            output<<"\n";
        }
        if (!step.fragment.children.empty()) {
            output << "    compiled method fragment:\n";
            const auto printMethod=[&](const auto& self,
                const SolvePlanNode& node,int depth)->void {
                output << std::string(static_cast<std::size_t>(depth)*2,' ')
                       << toString(node.kind);
                if (!node.operation.empty()) output << " -> " << node.operation << " provider=" << node.provider;
                output << "\n";
                for (const auto& child:node.children) self(self,child,depth+1);
            };
            printMethod(printMethod,step.fragment,3);
        } else {
            output << "    compiled backend: " << step.backendOperation << "\n";
        }
        if (!step.temporalMethod.empty())
            output << "    temporal method: " << step.temporalMethod << "\n";
        if(const auto* block=std::any_cast<CompiledLinearBlock>(&step.providerContract)) {
            output<<"    realization: whole steady scalar block; no temporal recipe participates\n    columns: unknown-major, physical-node index; GlobalDof ScalarTransport identity\n    boundary: existing nodal fixedValue identity row / zeroGradient mirror-ghost elimination\n";
            for(std::size_t i=0;i<block->equations.size();++i)output<<"      equation="<<block->equations[i]<<" row target="<<block->rowTargets[i]<<"; cross columns extracted from the same WHAT residual\n";
            output<<"    linear method: HYPRE FlexGMRES / BoomerAMG max="<<block->linear.maxIterations<<" relative="<<block->linear.relativeTolerance<<" absolute="<<block->linear.absoluteTolerance<<"\n    order: assemble -> whole solve -> original-matrix residual acceptance -> all-target publication; physical clock unchanged\n";
            for(const auto& read:block->valueReads)output<<"      block input="<<read.symbol<<" version="<<toString(read.version)<<"\n";
        }
        if(const auto* closure=std::any_cast<CompiledStageClosure>(&step.providerContract)) {
            output<<"    realization: one nonlocal closure group, not a time-integrated unknown\n"
                <<"    authored Stage operation: "<<step.backendOperation<<"; final Current refresh: "<<step.temporalPublish<<"\n    borrowed closure inputs:";
            for(const auto& id:closure->targets)output<<" "<<id;
            output<<"\n    shared outputs:";for(const auto& id:step.coefficientStates)output<<" "<<id;
            output<<"\n    original gradient/wall-distance/configured SST implementation; scalar boundary resources bind independently\n";
        }
        if(const auto* affine=std::any_cast<CompiledAffineEquation>(&step.providerContract)) {
            output<<"    realization: local affine residual A*target+R=0; target=-R/A; zero/non-finite A is an error\n"
                <<"    execution scope: "<<(step.stageLocal?"each declared recipe stage; publish relation after physical unknowns":"authored physical occurrence; no time integration")<<"\n";
            if(step.stageLocal)output<<"    physical relation publication: "<<step.temporalPublish<<" (same equation, committed physical inputs)\n";
            output<<"    output binding: "<<storageName(system.executableSystem.state.at(step.target.symbol).storageBinding)
                <<"; publication by HOW step="<<step.source.occurrence<<"; no hidden coefficient refresh\n";
            if(affine->boundaryPreparation==WorkspaceBoundary::NearestInterior)
                output<<"    spatial publication: internal final values -> nearest-interior extension over allocated ghost region -> readable; serial single patch, not MPI halo\n";
            else output<<"    spatial publication: internal nodes only; no implicit boundary preparation\n";
            for(const auto& read:affine->valueReads)output<<"      affine input="<<read.symbol<<" version="<<toString(read.version)<<"\n";
        }
        if(const auto* scalar=std::any_cast<CompiledScalarEquation>(&step.providerContract)) {
            output<<"    compiled RHS normalization: "<<scalar->normalized<<" (authored order preserved)\n";
            output<<"    target-relative structure: "<<(scalar->nonlinear?"nonlinear explicit values; no implicit Jacobian":"affine candidate")<<"\n";
            for(const auto& d:scalar->diffusions) {
                output<<"      diffusion occurrence="<<d.occurrence<<" input="<<d.input
                    <<(scalar->conservativeTransport && !scalar->primitiveStorage?"/rho":"")<<" version=Stage coefficient=";
                if(!d.coefficientState.empty())output<<d.coefficientState;
                else if(scalar->diffusivityState.empty())output<<d.coefficient;else output<<scalar->diffusivityState<<"(Stage)";
                output<<" components="<<d.components<<" implementation="
                    <<(scalar->arithmeticFaces?"Central2 arithmetic faces; nodal coefficient; physical-node stencil":scalar->conservativeTransport?"rho-weighted Central2 concentration":"Central2 Cartesian")<<" halo=1\n";
            }
            for(const auto& read:scalar->valueReads) {
                output<<"      numerical input="<<read.symbol<<" version="<<toString(read.version);
                for(const auto& publisher:system.solvePlan.compiledProgram.steps)
                    if(publisher.stageLocal && publisher.target.symbol==read.symbol && read.version==StateVersion::Stage)
                        output<<" publisher="<<publisher.source.occurrence<<" equation="<<publisher.source.equation
                            <<" operation="<<publisher.backendOperation;
                output<<"\n";
            }
            for(const auto& binding:scalar->lowering)output<<"      value lowering: "<<binding<<"\n";
        }
        if (!step.operatorBindings.empty()) {
            output << "    spatial/source providers:";
            for (const auto& provider:step.operatorBindings)
                output << " " << provider;
            output << "\n";
        }
        for (const auto& read:step.viewReads) output << "    STATE view read: " << read.symbol << " kind=" << toString(read.kind) << " storage=" << read.storage << " owner=NumericalProvider\n";
        for (const auto& use:step.stateUses) output << "    STATE read: " << targetText({use.symbol,use.kind}) << " version=" << toString(use.version) << (use.halo ? " halo-required" : " owner/view") << (use.perIteration ? " frozen-per-enclosing-iteration" : "") << "\n";
        for (const auto& effect:step.stateEffects) output << "    STATE write: " << effect.symbol << (effect.publish ? " publication/snapshot" : " invalidates dependent materialized views") << (effect.halo ? " with halo publication" : " without halo publication") << "\n";
        for (const auto& requirement:step.capabilityRequirements) output << "    capability: " << requirement.required << " when=" << requirement.when << " equation=" << requirement.equation << " fluid-port=" << requirement.fluidPort << " boundary=" << requirement.boundary << " : " << requirement.reason << "\n";
        output
               << "    output symbols:";
        for (const auto& symbol:step.writes)
            output << " " << symbol;
        output << "\n";
        if (!step.target.workspace.empty())
            output << "    output workspace: "
                   << step.target.workspace << "\n";
        if (!step.workspaceRequires.empty()) {
            output << "    requires workspace:";
            for (const auto& id:step.workspaceRequires) output << " " << id;
            output << "\n";
        }
        for (const auto& id:step.workspaceInvalidates) output << "    invalidates workspace: " << id << "\n";
        for (const auto& id:step.commitRequires) output << "    Commit requires publication: " << id << "\n";
        if (!step.workspaceProvides.empty()) {
            output << "    provides workspace:";
            for (const auto& id:step.workspaceProvides) output << " " << id;
            output << "\n";
        }
        if (!step.sourceMathInputs.empty()) {
            output << "    source mathematics (Equation):";
            for (const auto& id:step.sourceMathInputs) output << " " << id;
            output << "\n";
        }
        for (const auto& requirement:step.requirements)
            output << "    requires: " << requirement << "\n";
        if (!step.fragment.children.empty())
            output << "    execution: selected provider [Implemented numerical backend]"
                      " (method owns operation order)\n";
        else if (step.equationMethod=="DirectEvaluation")
            output << "    execution: generic direct Equation method\n";
        else if (step.temporalResidual)
            output << "    execution: selected temporal stage provider\n";
        else
            output << "    execution: registered relation provider\n";
    }
    if (std::any_of(system.solvePlan.compiledProgram.steps.begin(),
            system.solvePlan.compiledProgram.steps.end(),
            [](const CompiledEquationCall& step) {
                return step.equationMethod=="PressureMomentum";
            }))
        output << "  execution authority: ordered ExecutionProgram.root\n";

    output << "SYSTEM / STATE / BACKEND DETAILS\n"

           << "  template origin : "
           << (system.classification.defaultFluidPresetIncluded
               ? toString(system.classification.templateOrigin) : "none (explicit contributions)")
           << "\n"
           << "  density behavior: "
           << system.classification.densityBehavior << "\n"
           << "  thermo compress.: "
           << system.classification.thermodynamicCompressibility
           << "\n";

    // REGISTERED MODELS / COUPLING：注册状态必须显式可见，禁止静默忽略。
    output << "\nREGISTERED MODELS\n";
    if (system.coupling.preset.empty()) {
        output << "  coupling        : inactive\n"
               << "  reason          : no coupling preset registered\n";
    } else {
        output << "  coupling " << system.coupling.preset << " : "
               << toString(system.coupling.status) << "\n"
               << "  reason          : " << system.coupling.reason << "\n";
        if (!system.coupling.requirements.empty()) {
            output << "  requirements    :";
            for (const auto& requirement : system.coupling.requirements) {
                output << " " << requirement << ";";
            }
            output << "\n";
        }
        if (!system.coupling.derivedEquations.empty()) {
            output << "  derived equations:";
            for (const auto& equation : system.coupling.derivedEquations) {
                output << " " << equation;
            }
            output << "\n";
        }
        if (!system.coupling.derivedOperations.empty()) {
            output << "  derived operations:";
            for (const auto& operation : system.coupling.derivedOperations) {
                output << " " << operation;
            }
            output << "\n";
        }
    }

    // STATE role summary is independent of HOW order.
    output << "\nSTATE REALIZATION\n"
           << "  conservative transported mass : "
           << (system.realization.conservativeTransportedMass ? "yes" : "no")
           << "\n"
           << "  conservative momentum         : "
           << (system.realization.conservativeMomentum ? "yes" : "no")
           << "\n"
           << "  pressure as multiplier        : "
           << (system.realization.pressureMultiplier ? "yes" : "no")
           << "\n"
           << "  thermodynamic pressure closure: "
           << (system.realization.thermodynamicPressure ? "yes" : "no")
           << "\n"
           << "  phase transported state       : "
           << (system.realization.phaseTransportedState ? "yes" : "no")
           << "\n";
    const auto printRoleGroup = [&output](const char* label,
                                          const auto& groups) {
        if (groups.empty()) return;
        output << "  " << label << " :";
        for (const auto& group : groups) {
            output << " " << group.id << "[" << group.components
                   << "@" << group.componentOffset << "," << group.storageKey;
            if (group.constantValue) output << "=" << *group.constantValue;
            output << "]";
        }
        output << "\n";
    };
    printRoleGroup("transported",system.realization.transported);
    printRoleGroup("derived",system.realization.derived);
    printRoleGroup("multipliers",system.realization.multipliers);

    const auto& time=system.numericalSystem.time.recipe;
    output << "\nCONFIGURED TIME RECIPE (tableau only; actual execution above)\n"
           << "  id       : " << FDM::toString(time.id()) << "\n"
           << "  family   : " << FDM::toString(time.family()) << "\n"
           << "  order    : " << time.order() << "\n"
           << "  topology : " << FDM::toString(time.topology()) << "\n"
           << "  stages   : " << time.stageCount() << "\n";
    for (int index=0;index<time.stageCount();++index) {
        const auto& stage=time.stage(index);
        output << "    stage " << index+1 << ": c=" << stage.abscissa
               << " base=" << stage.baseWeight
               << " increment=" << stage.incrementWeight << "\n";
    }
    if (time.id()==FDM::TimeRecipeId::ClassicalRK4) {
        const auto weights=time.finalWeights();
        output << "    final weights: " << weights[0] << "," << weights[1]
               << "," << weights[2] << "," << weights[3]
               << " / " << time.finalDivisor() << "\n";
    }

    printNumericalSystem(output,system.numericalSystem,system.solvePlan.compiledProgram);

    printContributions(output,system);
    printRawSystem(output,system);
    printTransformations(output,system);
    printCompiledBindings(output,system);

    output << "\nEXECUTABLE EQUATION SYSTEM\n";

    output << "\nEXECUTABLE OPERATIONS\n";
    if (system.executableSystem.operations.empty()) output << "  (none)\n";
    for (const auto& operation : system.executableSystem.operations) {
        output << "  " << operation.operation << "\n";
        for (const auto requirement : operation.requirements) {
            output << "    requires: " << toString(requirement) << "\n";
        }
    }
    output << "\nOPERATION BINDINGS\n";
    if (system.runtime.operationBindings.empty()) output << "  (none)\n";
    for (const auto& binding : system.runtime.operationBindings) {
        output << "  " << binding.operation << "\n"
               << "    provider : "
               << (binding.provider.empty() ? "none" : binding.provider) << "\n"
               << "    status   : "
               << (binding.status == BindingStatus::Resolved
                   ? "Resolved" : binding.status == BindingStatus::Invalid
                       ? "Invalid" : "Unsupported") << "\n";
        if (!binding.reason.empty()) {
            output << "    reason   : " << binding.reason << "\n";
        }
    }

    printUnknownSection(output, "STATE VARIABLES", system, isStateUnknown);
    printUnknownSection(output, "AUXILIARY UNKNOWNS", system,
                        isAuxiliaryUnknown);
    printUnknownSection(output, "SOLID VARIABLES", system, isSolidUnknown);
    printUnknownSection(output, "CONSTRAINT VARIABLES", system,
                        isConstraintUnknown);
    printAlgebraicConstraints(output, system);

    printCompiledPlan(output,system);
    output << "\nFORMULA CALL EXECUTION\n";
    bool hasEquationCalls=false;
    const auto printCallBindings=[&](const auto& self,const SolvePlanNode& node)
        -> void {
        for (const CompiledMathRef& call:node.equationCalls) {
            hasEquationCalls=true;
            output << "  " << call.equation << " -> " << node.operation
                   << " provider="
                   << (node.provider.empty() ? "unbound" : node.provider) << "\n";
        }
        for (const auto& child:node.children) self(self,child);
    };
    printCallBindings(printCallBindings,system.solvePlan.root);
    if (hasEquationCalls) {
        output << "  lowering: frozen operation bindings; AST compilation belongs to the selected numerical provider\n";
    } else {
        output << "  (none)\n";
    }
    output << "\nEXECUTION REQUIREMENTS\n";
    for (const auto& requirement : system.runtime.requirements) {
        output << "  " << requirement.name << "  "
               << requirementStatus(requirement, system) << "\n";
    }
    output << "\nWORKSPACE REQUIREMENTS\n";
    if (system.runtime.workspaceRequirements.empty()) {
        output << "  (none)\n";
    } else {
        for (const auto& workspace : system.runtime.workspaceRequirements) {
            output << "  " << workspace.id
                   << "  [" << toString(workspace.location) << ", "
                   << toString(workspace.ownership) << ", components="
                   << workspace.components << "]\n";
        }
    }
    output << "\nNUMERICAL PROVIDERS\n";
    if (system.runtime.providerRequirements.empty()) output << "  (none)\n";
    for (const auto& provider : system.runtime.providerRequirements) {
        output << "  " << provider.id << "  "
               << provider.responsibility << "\n";
    }
    output << "\nRUNTIME SERVICES\n";
    if (system.runtime.runtimeServiceRequirements.empty()) output << "  (none)\n";
    for (const auto& service : system.runtime.runtimeServiceRequirements) {
        output << "  " << service.id << "  " << service.reason << "\n";
    }
    output << "============================================================";
    return output.str();
}

} // namespace SF::System
