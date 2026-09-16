/// @file SF_equationContribution.cpp
/// @brief Unified equation contribution registration.

#include "SF_equationContribution.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::System {

void EquationSystemBuilder::addUnknown(UnknownDescriptor unknown) {
    if (unknown.id.empty()) throw std::runtime_error("Unknown contribution requires an id.");
    if (hasUnknown(system_,unknown.id)) {
        throw std::runtime_error("Duplicate unknown contribution '"+unknown.id+"'.");
    }
    system_.unknowns.push_back(std::move(unknown));
}

void EquationSystemBuilder::addEquation(
        EquationDescriptor descriptor, Equation::Definition definition) {
    if (descriptor.id != definition.name) {
        throw std::runtime_error(
            "Executable equation id does not match descriptor '"+descriptor.id+"'.");
    }
    if (hasEquation(system_,descriptor.id)) {
        throw std::runtime_error("Duplicate equation contribution '"+descriptor.id+"'.");
    }
    system_.equations.push_back(std::move(descriptor));
    system_.equationDefinitions.add(std::move(definition));
}

void EquationSystemBuilder::addTerm(
        const std::string& equationId, Equation::Term term) {
    system_.equationDefinitions.addRightTerm(equationId,std::move(term));
}

void EquationSystemBuilder::addConstraint(ConstraintDescriptor constraint) {
    if (hasConstraint(system_,constraint.id)) {
        throw std::runtime_error("Duplicate constraint contribution '"+constraint.id+"'.");
    }
    system_.constraints.push_back(std::move(constraint));
}

void EquationSystemBuilder::addSolveBlock(SolveBlock block) {
    if (hasSolveBlock(system_,block.id)) {
        throw std::runtime_error("Duplicate solve-block contribution '"+block.id+"'.");
    }
    system_.solveBlocks.push_back(std::move(block));
}

void EquationSystemBuilder::addClosure(std::string closure) {
    if (closure.empty()) throw std::runtime_error("Closure contribution requires a name.");
    system_.closures.push_back(std::move(closure));
}

void EquationSystemBuilder::require(ExecutionRequirement requirement) {
    if (hasRequirement(system_,requirement.name)) {
        throw std::runtime_error(
            "Duplicate execution requirement '"+requirement.name+"'.");
    }
    system_.requirements.push_back(std::move(requirement));
}

} // namespace SF::System
