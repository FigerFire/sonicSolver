/// @file SF_systemValidator.cpp
/// @brief 启动阶段统一拒绝未声明未知量、悬空约束和缺失执行能力。

#include "SF_systemValidator.h"

#include <set>
#include <stdexcept>

namespace SF::System {

void validate(const ResolvedSimulationSystem& system) {
    std::set<std::string> unknowns;
    for (const auto& unknown : system.unknowns) {
        if (unknown.id.empty() || unknown.name.empty()
            || unknown.components <= 0
            || !unknowns.insert(unknown.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate unknown.");
        }
    }
    std::set<std::string> equations;
    for (const auto& equation : system.equations) {
        if (equation.id.empty() || equation.name.empty()
            || !equations.insert(equation.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate equation.");
        }
        for (const std::string& unknown : equation.solvedUnknowns) {
            if (unknowns.find(unknown) == unknowns.end()) {
                throw std::runtime_error(
                    "Equation '"+equation.id
                    +"' references undeclared unknown '"+unknown+"'.");
            }
        }
    }
    std::set<std::string> constraints;
    for (const auto& constraint : system.constraints) {
        if (constraint.id.empty() || constraint.equation.empty()
            || (!constraint.multiplierUnknown.empty()
                && unknowns.find(constraint.multiplierUnknown) == unknowns.end())
            || !constraints.insert(constraint.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid multiplier constraint.");
        }
    }
    std::set<std::string> blocks;
    for (const auto& block : system.solveBlocks) {
        if (block.id.empty() || block.strategy.empty()
            || !blocks.insert(block.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate solve block.");
        }
        for (const auto& equation : block.equations) {
            if (equations.find(equation) == equations.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared equation '"+equation+"'.");
            }
        }
        for (const auto& constraint : block.constraints) {
            if (constraints.find(constraint) == constraints.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared constraint '"+constraint+"'.");
            }
        }
        for (const auto& unknown : block.unknowns) {
            if (unknowns.find(unknown) == unknowns.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared unknown '"+unknown+"'.");
            }
        }
    }
    for (const auto& requirement : system.requirements) {
        if (requirement.required && !requirement.available) {
            throw std::runtime_error(
                "Resolved system requires unavailable capability '"
                +requirement.name+"': "+requirement.detail+".");
        }
    }
}

} // namespace SF::System
