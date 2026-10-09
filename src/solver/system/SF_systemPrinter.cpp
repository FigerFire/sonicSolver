/// @file SF_systemPrinter.cpp
/// @brief 打印 case 最终数学身份，不参与求解和状态修改。

#include "SF_systemPrinter.h"
#include "SF_eulerianAssembly.h"
#include "SF_eulerianTurbulence.h"
#include "SF_immersedMethods.h"
#include "SF_couplingStatus.h"

#include "SF_config.h"
#include "SF_couplingStatus.h"

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

bool isAlgebraicEquation(const EquationDescriptor& equation) {
    return equation.category != EquationCategory::PhysicalEquation;
}

bool isAuxiliaryUnknown(const StateSymbol& unknown) {
    return unknown.id == "p" || unknown.id == "pPrime"
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

std::string constraintName(const EquationDescriptor& equation) {
    if (isPressureConstraint(equation.id)) return "incompressibility";
    return equation.name;
}

std::string termText(const SF::Equation::Term& term) {
    switch (term.kind) {
        case SF::Equation::TermKind::Transient:
            return "ddt("+term.primary.name+")";
        case SF::Equation::TermKind::Divergence:
            return "div("+term.primary.name+")";
        case SF::Equation::TermKind::Gradient:
            return "grad("+term.primary.name+")";
        case SF::Equation::TermKind::Diffusion:
            return "diffusion("+term.secondary.name+","+term.primary.name+")";
        case SF::Equation::TermKind::Source:
            return "source("+term.primary.name+")";
        case SF::Equation::TermKind::Constraint:
            return "constraint("+term.primary.name+")";
        case SF::Equation::TermKind::AlgebraicRelation:
            return "algebraic("+term.primary.name+")";
    }
    return "unknownTerm";
}

std::string expressionText(const SF::Equation::Expression& expression) {
    std::string result;
    for (const auto& term : expression.terms) {
        if (!result.empty()) result += " + ";
        result += termText(term);
    }
    return result.empty() ? "0" : result;
}

const char* termKindName(SF::Equation::TermKind kind) {
    switch (kind) {
        case SF::Equation::TermKind::Transient: return "transient";
        case SF::Equation::TermKind::Divergence: return "divergence";
        case SF::Equation::TermKind::Gradient: return "gradient";
        case SF::Equation::TermKind::Diffusion: return "diffusion";
        case SF::Equation::TermKind::Source: return "source";
        case SF::Equation::TermKind::Constraint: return "constraint";
        case SF::Equation::TermKind::AlgebraicRelation: return "algebraic";
    }
    return "unknown";
}

void printNumericalSystem(
        std::ostringstream& output,
        const CompiledNumericalSystem& numerical) {
    output << "\nCOMPILED NUMERICAL SYSTEM\n"
           << "  required halo width : " << numerical.requiredHaloWidth << "\n";
    if (numerical.pressureFaceCoupling != PressureFaceCoupling::None)
        output << "  pressure face coupling : RhieChow\n";
    if (numerical.pressureOperator) {
        output << "  pressure temporal treatment : one physical step; first predictor R(U_n); "
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
                   << " -> " << term.output << " " << term.operatorName
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

void printPhysicalEquations(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nLEGACY COMPATIBILITY PHYSICAL EQUATIONS\n";
    bool printed = false;
    for (const auto& equation : system.executableSystem.legacyEquations) {
        if (isAlgebraicEquation(equation)) continue;
        output << "  " << equation.id << "  " << equation.name;
        if (!equation.kind.empty()) output << "  {" << equation.kind << "}";
        const auto& formula=system.executableSystem.registry.at(equation.id);
        output << "\n    " << formulaText(formula) << "\n";
        printed = true;
    }
    if (!printed) output << "  (none)\n";
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
    for (const auto& equation : system.executableSystem.legacyEquations) {
        if (equation.kind == "constraint" || isPressureConstraint(equation.id)) {
            emit(equation.id, constraintName(equation));
        }
    }
    for (const auto& constraint : system.executableSystem.constraints) {
        emit(constraint.id, constraint.name);
    }
    if (!printed) output << "  (none)\n";
}

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

bool hasAlgebraicContent(
        const SolveBlock& block,
        const ResolvedSimulationSystem& system) {
    for (const auto& id : block.equations) {
        if (isPressureConstraint(id)) return true;
        const auto equation = std::find_if(
            system.executableSystem.legacyEquations.begin(),
            system.executableSystem.legacyEquations.end(),
            [&](const EquationDescriptor& value) { return value.id == id; });
        if (equation != system.executableSystem.legacyEquations.end()
            && equation->category != EquationCategory::AlgorithmicDerivedEquation
            && isAlgebraicEquation(*equation)) {
            return true;
        }
    }
    return !block.constraints.empty();
}

std::string algebraicSystemId(const SolveBlock& block) {
    if (block.id.find("KKT") != std::string::npos) return "KKT_FLUID_IBM";
    if (block.id == kPressureScheduleId) return "PRESSURE_CORRECTION";
    return block.id;
}

std::string algebraicRowLabel(
        const std::string& id,
        const ResolvedSimulationSystem& system) {
    if (isPressureConstraint(id)) return "pressure/continuity constraint";
    for (const auto& constraint : system.executableSystem.constraints) {
        if (constraint.id == id) return constraint.name;
    }
    for (const auto& equation : system.executableSystem.legacyEquations) {
        if (equation.id == id) return equation.name;
    }
    return id;
}

void printAlgebraicSystems(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    std::vector<const SolveBlock*> blocks;
    for (const auto& block : system.solvePlan.blocks) {
        if (hasAlgebraicContent(block, system)) blocks.push_back(&block);
    }
    if (blocks.empty()) return;

    output << "\nALGEBRAIC SYSTEM\n";
    for (const SolveBlock* block : blocks) {
        output << "  " << algebraicSystemId(*block) << "\n";
        std::vector<std::string> rows;
        const auto appendRow = [&](const std::string& id) {
            const std::string key = constraintId(id);
            const bool duplicate = std::any_of(
                rows.begin(),rows.end(),
                [&](const std::string& existing) {
                    return constraintId(existing) == key;
                });
            if (!duplicate) rows.push_back(id);
        };
        for (const auto& id : block->equations) appendRow(id);
        for (const auto& id : block->constraints) {
            appendRow(id);
        }

        std::stable_sort(rows.begin(), rows.end(), [&](const std::string& left,
                                                        const std::string& right) {
            const auto rank = [](const std::string& id) {
                if (isPressureConstraint(id)) return 1;
                return 3;
            };
            return rank(left) < rank(right);
        });
        int row = 1;
        for (const auto& id : rows) {
            output << "    row " << row++ << " : "
                   << algebraicRowLabel(id, system) << "\n";
        }
    }
}

void printSolveBlocks(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nSOLVE BLOCKS\n";
    for (const auto& block : system.solvePlan.blocks) {
        output << "  " << block.id << "  " << block.name
               << "  strategy=" << block.strategy
               << "  kind=" << FDM::toString(block.strategyKind) << "\n";
    }
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
               << contribution.id;
        if (!contribution.name.empty()) output << "  " << contribution.name;
        output << "\n";
    }
    for (const auto& modification : system.rawSystem.modifications) {
        output << "  modification: " << toString(modification.kind) << " "
               << modification.targetKind << " " << modification.targetId
               << "  origin=" << toString(modification.origin.kind) << "\n";
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
    output << "\nRAW EQUATION SYSTEM\n  Legacy compatibility descriptors:\n";
    if (system.rawSystem.legacyEquations.empty()) output << "  (none)\n";
    for (const auto& equation : system.rawSystem.legacyEquations) {
        output << "  " << equation.id << "  " << equation.name
               << "  category=" << toString(equation.category)
               << "  origin=" << toString(equation.origin.kind) << "\n";
    }
    output << "  constraints:\n";
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
    if (!node.operation.empty()) {
        output << "  operation=" << node.operation << " provider=" << node.provider;
        output << " legacyAdapter=" << (node.legacyAdapter ? "true" : "false");
        if (node.legacyAdapter) output << " [Legacy adapter]";
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
    output << "\nWHAT / EQUATIONS\n";
    for (const Equation& formula:system.executableSystem.registry.entries()) {
        output << "  " << formula.id << ": " << formulaText(formula) << "\n";
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
    output << "\nHOW / EXECUTION\n";
    const auto printHow=[&](const auto& self,const ExecutionScope& node,int depth)->void {
        output << std::string(static_cast<std::size_t>(depth)*2,' ')
               << node.order << " " << toString(node.kind);
        if (node.kind==ExecutionKind::EquationCall) {
            output << " " << node.step.equation;
            if (!node.step.occurrence.empty()) output << "@" << node.step.occurrence;
            output << " -> " << targetText(node.step.target);
        } else if (!node.id.empty()) {
            output << " " << node.id;
        }
        if (node.kind==ExecutionKind::Loop)
            output << " limit=" << node.repetitions;
        if (!node.terminationSignal.empty())
            output << " until=" << node.terminationSignal;
        output << " [" << node.origin.source << "]\n";
        for (const auto& child:node.children) self(self,child,depth+1);
    };
    printHow(printHow,system.solvePlan.sourceProgram.root,1);
    if (system.solvePlan.compiledProgram.steps.empty())
        output << "  Legacy numerical backend / migration pending\n";
    if (!system.solvePlan.sourceProgram.legacyEntries.empty()) {
        output << "  Legacy execution authority (pending entries are not scheduled twice):\n";
        for (const auto& entry:system.solvePlan.sourceProgram.legacyEntries)
            printHow(printHow,entry,2);
    }
    output << "\nWHICH / NUMERICS\n";
    output << "  time: " << FDM::toString(system.numericalSelection.recipes.time.id()) << "\n";
    for (const auto& binding:system.numericalSelection.bindings) {
        for (const auto& parameter:binding.parameters)
            output << "    " << binding.equation << "." << parameter.first << "=" << parameter.second << "\n";
        output << "  " << binding.equation;
        if (!binding.occurrence.empty()) output << "@" << binding.occurrence;
        output << " -> " << binding.method << "\n";
    }
    if (system.numericalSelection.pressureOperator) {
        output << "  linear: HYPRE (selected pressure linear configuration)\n"
               << "  face coupling: RhieChow\n";
    }
    if (!system.numericalSelection.legacySpatialInputs.empty())
        output << "  Legacy numerical backend: explicit spatial input adapter\n";
    output << "  frozen method -> runtime provider bindings:\n";
    for (const auto& call:system.solvePlan.compiledProgram.steps)
        output << "    " << call.source.equation << "@" << call.source.occurrence
               << " method=" << call.equationMethod << " provider=" << call.backendProvider << "\n";
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
               << "  snapshots -> all prepare -> all RHS -> all advance -> validate all -> publish -> time.commit\n";
        for(const auto& p:system.solvePlan.compiledProgram.temporalParticipants) {
            output<<"  participant "<<p.identity<<" provider="<<p.provider<<"\n    calls:";
            for(const auto& c:p.calls)output<<" "<<c.equation<<"->"<<c.target;
            output<<"\n    phases: "<<p.snapshot<<" / "<<p.prepareStage<<" / "<<p.rhs<<" / "<<p.advance<<" / "<<p.publish<<"\n";
        }
    }
    output << "\nCOMPILED OCCURRENCES / STORAGE / PROVIDERS\n";
    for (const auto& step:system.solvePlan.compiledProgram.steps) {
        output << "  " << step.source.equation << " -> "
               << targetText(step.source.target) << "\n"
               << "    equation method: " << step.equationMethod << "\n"
               << "    selected runtime provider: " << step.backendProvider << "\n";
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
        if (!step.operatorBindings.empty()) {
            output << "    spatial/source providers:";
            for (const auto& provider:step.operatorBindings)
                output << " " << provider;
            output << "\n";
        }
        for (const auto& use:step.stateUses) output << "    STATE read: " << use.symbol << " version=" << toString(use.version) << (use.halo ? " halo-required" : " owner/view") << (use.perIteration ? " frozen-per-enclosing-iteration" : "") << "\n";
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
    output << "\nTIME RECIPE\n"
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

    printNumericalSystem(output,system.numericalSystem);

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
    printPhysicalEquations(output, system);
    printAlgebraicConstraints(output, system);
    printAlgebraicSystems(output, system);

    printSolveBlocks(output, system);
    printCompiledPlan(output,system);
    output << "\nFORMULA CALL EXECUTION\n";
    bool hasFormulaCalls=false;
    const auto printCallBindings=[&](const auto& self,const SolvePlanNode& node)
        -> void {
        for (const CompiledMathRef& call:node.equationCalls) {
            hasFormulaCalls=true;
            output << "  " << call.equation << " -> " << node.operation
                   << " provider="
                   << (node.provider.empty() ? "unbound" : node.provider) << "\n";
        }
        for (const auto& child:node.children) self(self,child);
    };
    printCallBindings(printCallBindings,system.solvePlan.root);
    if (hasFormulaCalls) {
        output << "  lowering: operation-bound production numerical provider; "
                  "generic FormulaCompiler is not yet this path\n";
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
