/// @file SF_systemPrinter.cpp
/// @brief 打印 case 最终数学身份，不参与求解和状态修改。

#include "SF_systemPrinter.h"

#include "SF_config.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <vector>

namespace SF::System {
namespace {

bool isPressureConstraint(const std::string& id) {
    return id == "E_PRESSURE" || id == "E_SHARED_PRESSURE";
}

bool isAlgebraicEquation(const EquationDescriptor& equation) {
    return equation.kind == "constraint"
        || isPressureConstraint(equation.id)
        || equation.id == "E_IBM_STATIONARITY";
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
           << unknown.components << "]\n";
    return output.str();
}

void printUnknownSection(
        std::ostringstream& output,
        const char* title,
        const ResolvedSimulationSystem& system,
        bool (*predicate)(const UnknownDescriptor&)) {
    output << "\n" << title << "\n";
    bool printed = false;
    for (const auto& unknown : system.unknowns) {
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

void printPhysicalEquations(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nPHYSICAL EQUATIONS\n";
    bool printed = false;
    for (const auto& equation : system.equations) {
        if (isAlgebraicEquation(equation)) continue;
        output << "  " << equation.id << "  " << equation.name;
        if (!equation.kind.empty()) output << "  {" << equation.kind << "}";
        output << "\n";
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
    for (const auto& equation : system.equations) {
        if (equation.kind == "constraint" || isPressureConstraint(equation.id)) {
            emit(equation.id, constraintName(equation));
        }
    }
    for (const auto& constraint : system.constraints) {
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
            system.equations.begin(), system.equations.end(),
            [&](const EquationDescriptor& value) { return value.id == id; });
        if (equation != system.equations.end() && isAlgebraicEquation(*equation)) {
            return true;
        }
    }
    return !block.constraints.empty();
}

std::string algebraicSystemId(const SolveBlock& block) {
    if (block.id.find("KKT") != std::string::npos) return "KKT_FLUID_IBM";
    if (block.id == "S_PRESSURE") return "PRESSURE_CORRECTION";
    return block.id;
}

std::string algebraicRowLabel(
        const std::string& id,
        const ResolvedSimulationSystem& system) {
    if (id == "E_IBM_STATIONARITY") return "momentum stationarity";
    if (isPressureConstraint(id)) return "pressure/continuity constraint";
    if (id == "C_IBM_NO_SLIP") return "immersed no-slip constraint";
    for (const auto& constraint : system.constraints) {
        if (constraint.id == id) return constraint.name;
    }
    for (const auto& equation : system.equations) {
        if (equation.id == id) return equation.name;
    }
    return id;
}

void printAlgebraicSystems(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    std::vector<const SolveBlock*> blocks;
    for (const auto& block : system.solveBlocks) {
        if (hasAlgebraicContent(block, system)) blocks.push_back(&block);
    }
    if (blocks.empty()) return;

    output << "\nALGEBRAIC SYSTEM\n";
    for (const SolveBlock* block : blocks) {
        output << "  " << algebraicSystemId(*block) << "\n";
        std::vector<std::string> rows;
        const bool hasStationarity = contains(
            block->equations, "E_IBM_STATIONARITY");
        for (const auto& id : block->equations) {
            // The monolithic descriptor also carries the predictor momentum
            // id for ownership; stationarity is its actual KKT row.
            if (hasStationarity && id == "E_MOMENTUM") continue;
            if (!contains(rows, id)) rows.push_back(id);
        }
        for (const auto& id : block->constraints) {
            if (!contains(rows, id)) rows.push_back(id);
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

void printSolveStages(
        std::ostringstream& output,
        const ResolvedSimulationSystem& system) {
    output << "\nSOLVE STAGES\n";
    for (const auto& block : system.solveBlocks) {
        output << "  " << block.id << "  " << block.name
               << "  strategy=" << block.strategy << "\n";
    }
}

bool hasConstraintVariables(const ResolvedSimulationSystem& system) {
    return std::any_of(
        system.unknowns.begin(), system.unknowns.end(),
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

} // namespace

std::string describe(const ResolvedSimulationSystem& system) {
    std::ostringstream output;
    output << "\n============================================================\n"
           << "sonicSolver - Resolved Mathematical System\n"
           << "============================================================\n"
           << "FLOW\n"
           << "  formulation     : " << FDM::toString(system.formulation) << "\n"
           << "  template origin : " << toString(system.templateOrigin) << "\n"
           << "  time integrator : " << system.timeIntegrator << "\n";

    printUnknownSection(output, "STATE VARIABLES", system, isStateUnknown);
    printUnknownSection(output, "AUXILIARY UNKNOWNS", system,
                        isAuxiliaryUnknown);
    printUnknownSection(output, "SOLID VARIABLES", system, isSolidUnknown);
    printUnknownSection(output, "CONSTRAINT VARIABLES", system,
                        isConstraintUnknown);
    printPhysicalEquations(output, system);
    printAlgebraicConstraints(output, system);
    printAlgebraicSystems(output, system);

    if (!system.immersedAlgorithm.empty()) {
        output << "\nIBM\n"
               << "  algorithm      : " << system.immersedAlgorithm << "\n"
               << "  reference      : " << system.immersedReference << "\n"
               << "  support        : " << system.immersedSupport << "\n"
               << "  representation : " << system.immersedRepresentation << "\n"
               << "  enforcement    : " << system.immersedEnforcement << "\n"
               << "  solid          : " << system.immersedSolid << "\n";
        if (!system.immersedFunctional.empty()) {
            output << "  functional     : "
                   << system.immersedFunctional << "\n";
        }
    }
    printSolveStages(output, system);
    output << "\nEXECUTION REQUIREMENTS\n";
    for (const auto& requirement : system.requirements) {
        output << "  " << requirement.name << "  "
               << requirementStatus(requirement, system) << "\n";
    }
    output << "============================================================";
    return output.str();
}

} // namespace SF::System
