/// @file SF_assemblyPlan.cpp
/// @brief Executable equation lowering implementation.

#include "solver/equation/SF_assemblyPlan.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Equation {

bool AssemblyPlan::contains(TermKind kind) const {
    const auto has = [kind](const std::vector<const Term*>& terms) {
        return std::any_of(terms.begin(),terms.end(),
            [kind](const Term* term) { return term->kind == kind; });
    };
    return has(left) || has(right);
}

AssemblyPlanRegistry::AssemblyPlanRegistry(
        const System& definitions,
        const std::vector<std::string>& equationIds) {
    plans_.reserve(equationIds.size());
    for (const std::string& id : equationIds) {
        if (contains(id)) {
            throw std::runtime_error(
                "AssemblyPlanRegistry duplicate equation '" + id + "'.");
        }
        const Definition& definition = definitions.at(id);
        AssemblyPlan plan;
        plan.definition = &definition;
        for (const Term& term : definition.left.terms) {
            plan.left.push_back(&term);
        }
        for (const Term& term : definition.right.terms) {
            plan.right.push_back(&term);
        }
        plans_.push_back({id,std::move(plan)});
    }
}

const AssemblyPlan& AssemblyPlanRegistry::at(
        std::string_view equationId) const {
    const auto found = std::find_if(
        plans_.begin(),plans_.end(),[equationId](const Entry& entry) {
            return entry.id == equationId;
        });
    if (found == plans_.end()) {
        throw std::runtime_error(
            "No AssemblyPlan is registered for equation '"
            + std::string(equationId) + "'.");
    }
    return found->plan;
}

bool AssemblyPlanRegistry::contains(std::string_view equationId) const {
    return std::any_of(plans_.begin(),plans_.end(),
        [equationId](const Entry& entry) { return entry.id == equationId; });
}

std::vector<AssemblyPlan> makeAssemblyPlan(
        const System& definitions,
        const std::vector<std::string>& equationIds) {
    const AssemblyPlanRegistry registry(definitions,equationIds);
    std::vector<AssemblyPlan> result;
    result.reserve(equationIds.size());
    for (const std::string& id : equationIds) {
        result.push_back(registry.at(id));
    }
    return result;
}

} // namespace SF::Equation
