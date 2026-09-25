/// @file SF_systemPrinter.cpp
/// @brief 打印 case 最终数学身份，不参与求解和状态修改。

#include "SF_systemPrinter.h"
#include "SF_couplingStatus.h"

#include "SF_config.h"
#include "SF_couplingStatus.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <vector>

namespace SF::System {
namespace {

const char* roleName(UnknownRole role) {
    switch (role) {
        case UnknownRole::Primary: return "primary";
        case UnknownRole::Transported: return "transported";
        case UnknownRole::Algebraic: return "algebraic";
        case UnknownRole::Multiplier: return "multiplier";
        case UnknownRole::Derived: return "derived";
    }
    return "unknown";
}

const char* storageName(StorageBinding binding) {
    switch (binding) {
        case StorageBinding::PackedDistributed: return "packed-distributed";
        case StorageBinding::NamedDistributed: return "named-distributed";
        case StorageBinding::TransientWorkspace: return "transient-workspace";
        case StorageBinding::SpecializedExecutor: return "specialized-executor";
    }
    return "unknown";
}

bool isPressureConstraint(const std::string& id) {
    return id == "E_PRESSURE" || id == "E_SHARED_PRESSURE";
}

bool isAlgebraicEquation(const EquationDescriptor& equation) {
    return equation.category != EquationCategory::PhysicalEquation;
}

bool isAuxiliaryUnknown(const UnknownDescriptor& unknown) {
    return unknown.id == "p" || unknown.id == "pPrime"
        || (unknown.id.size() >= 3
            && unknown.id.compare(unknown.id.size() - 3, 3, "Aux") == 0);
}

bool isConstraintUnknown(const UnknownDescriptor& unknown) {
    return unknown.location == VariableLocation::BodyConstraint
        || unknown.location == VariableLocation::SurfaceConstraint;
}

bool isSolidUnknown(const UnknownDescriptor& unknown) {
    return unknown.location == VariableLocation::SolidGlobal;
}

std::string unknownDescription(const UnknownDescriptor& unknown) {
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
    output << "]\n";
    return output.str();
}

void printUnknownSection(
        std::ostringstream& output,
        const char* title,
        const ResolvedSimulationSystem& system,
        bool (*predicate)(const UnknownDescriptor&)) {
    output << "\n" << title << "\n";
    bool printed = false;
    for (const auto& unknown : system.executableSystem.unknowns) {
        if (predicate(unknown)) {
            output << unknownDescription(unknown);
            printed = true;
        }
    }
    if (!printed) output << "  (none)\n";
}

bool isStateUnknown(const UnknownDescriptor& unknown) {
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

std::string termText(const Equation::Term& term) {
    switch (term.kind) {
        case Equation::TermKind::Transient:
            return "ddt("+term.primary.name+")";
        case Equation::TermKind::Divergence:
            return "div("+term.primary.name+")";
        case Equation::TermKind::Gradient:
            return "grad("+term.primary.name+")";
        case Equation::TermKind::Diffusion:
            return "diffusion("+term.secondary.name+","+term.primary.name+")";
        case Equation::TermKind::Source:
            return "source("+term.primary.name+")";
        case Equation::TermKind::Constraint:
            return "constraint("+term.primary.name+")";
        case Equation::TermKind::AlgebraicRelation:
            return "algebraic("+term.primary.name+")";
    }
    return "unknownTerm";
}

std::string expressionText(const Equation::Expression& expression) {
    std::string result;
    for (const auto& term : expression.terms) {
        if (!result.empty()) result += " + ";
        result += termText(term);
    }
    return result.empty() ? "0" : result;
}

const char* termKindName(Equation::TermKind kind) {
    switch (kind) {
        case Equation::TermKind::Transient: return "transient";
        case Equation::TermKind::Divergence: return "divergence";
        case Equation::TermKind::Gradient: return "gradient";
        case Equation::TermKind::Diffusion: return "diffusion";
        case Equation::TermKind::Source: return "source";
        case Equation::TermKind::Constraint: return "constraint";
        case Equation::TermKind::AlgebraicRelation: return "algebraic";
    }
    return "unknown";
}

void printNumericalSystem(
        std::ostringstream& output,
        const CompiledNumericalSystem& numerical) {
    output << "\nCOMPILED NUMERICAL SYSTEM\n"
           << "  required halo width : " << numerical.requiredHaloWidth << "\n";
    if (numerical.terms.empty()) {
        output << "  bound terms         : (none)\n";
    } else {
        output << "  bound terms\n";
        for (const auto& term : numerical.terms) {
            output << "    " << term.equationId << "[" << term.ordinal << "] "
                   << termKindName(term.kind) << "(" << term.primary;
            if (!term.secondary.empty()) output << "," << term.secondary;
            output << ") -> " << FDM::toString(term.recipe.id())
                   << " [" << FDM::toString(term.recipe.temporalRole())
                   << ", halo=" << term.recipe.haloWidth() << "]\n";
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
    output << "\nPHYSICAL EQUATIONS\n";
    bool printed = false;
    for (const auto& equation : system.executableSystem.equations) {
        if (isAlgebraicEquation(equation)) continue;
        output << "  " << equation.id << "  " << equation.name;
        if (!equation.kind.empty()) output << "  {" << equation.kind << "}";
        const auto& definition = equationDefinition(system,equation.id);
        output << "\n    " << expressionText(definition.left)
               << " = " << expressionText(definition.right) << "\n";
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
    for (const auto& equation : system.executableSystem.equations) {
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
        if (id == "E_IBM_STATIONARITY" || isPressureConstraint(id)) return true;
        const auto equation = std::find_if(
            system.executableSystem.equations.begin(),
            system.executableSystem.equations.end(),
            [&](const EquationDescriptor& value) { return value.id == id; });
        if (equation != system.executableSystem.equations.end()
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
    if (id == "E_IBM_STATIONARITY") return "momentum stationarity";
    if (isPressureConstraint(id)) return "pressure/continuity constraint";
    if (id == "C_IBM_NO_SLIP") return "immersed no-slip constraint";
    for (const auto& constraint : system.executableSystem.constraints) {
        if (constraint.id == id) return constraint.name;
    }
    for (const auto& equation : system.executableSystem.equations) {
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
        const bool hasStationarity = contains(
            block->equations, "E_IBM_STATIONARITY");
        for (const auto& id : block->equations) {
            // The monolithic descriptor also carries the predictor momentum
            // id for ownership; stationarity is its actual KKT row.
            if (hasStationarity && id == "E_MOMENTUM") continue;
            appendRow(id);
        }
        for (const auto& id : block->constraints) {
            appendRow(id);
        }

        std::stable_sort(rows.begin(), rows.end(), [&](const std::string& left,
                                                        const std::string& right) {
            const auto rank = [](const std::string& id) {
                if (id == "E_IBM_STATIONARITY") return 0;
                if (isPressureConstraint(id)) return 1;
                if (id == "C_IBM_NO_SLIP") return 2;
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
        system.executableSystem.unknowns.begin(),
        system.executableSystem.unknowns.end(),
        [](const UnknownDescriptor& unknown) {
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
    output << "\nRAW EQUATION SYSTEM\n";
    if (system.rawSystem.equations.empty()) output << "  (none)\n";
    for (const auto& equation : system.rawSystem.equations) {
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
    if (!node.operation.empty()) {
        output << "  operation=" << node.operation;
    }
    output << "\n";
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
                ? "runnable" : "unsupported") << "\n";
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
           << "FLOW (WHAT)\n"
           << "  template origin : "
           << toString(system.classification.templateOrigin)
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

    // STATE REALIZATION：由 executable equations 的 role 导出。
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
                   << "@" << group.componentOffset << "," << group.storageKey
                   << "]";
        }
        output << "\n";
    };
    printRoleGroup("transported",system.realization.transported);
    printRoleGroup("derived",system.realization.derived);
    printRoleGroup("multipliers",system.realization.multipliers);

    output << "\nTIME RECIPE\n"
           << "  id       : " << FDM::toString(system.timeRecipe.id()) << "\n"
           << "  family   : " << FDM::toString(system.timeRecipe.family()) << "\n"
           << "  order    : " << system.timeRecipe.order() << "\n"
           << "  topology : " << FDM::toString(system.timeRecipe.topology()) << "\n"
           << "  stages   : " << system.timeRecipe.stageCount() << "\n";

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
                   ? "Resolved" : "Unsupported") << "\n";
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
