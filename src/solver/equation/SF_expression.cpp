/// @file SF_expression.cpp
/// @brief 强类型方程项表达式的存储、校验与遍历实现。

#include "core/system/SF_expression.h"

#include <algorithm>

namespace SF::Equation {

void Definition::validate() const {
    if (name.empty()) {
        throw std::runtime_error("Equation definition requires a non-empty name.");
    }
    int transientCount = 0;
    for (const Term& term : left.terms) {
        if (term.primary.name.empty()) {
            throw std::runtime_error(
                "Equation '" + name + "' contains an unnamed left-hand term.");
        }
        transientCount += term.kind == TermKind::Transient ? 1 : 0;
    }
    for (const Term& term : right.terms) {
        if (term.primary.name.empty()) {
            throw std::runtime_error(
                "Equation '" + name + "' contains an unnamed right-hand term.");
        }
        transientCount += term.kind == TermKind::Transient ? 1 : 0;
    }
    if (transientCount > 1) {
        throw std::runtime_error(
            "Equation '" + name + "' contains more than one ddt term.");
    }
}

void System::add(Definition definition) {
    definition.validate();
    const auto duplicate = std::find_if(
        equations_.begin(), equations_.end(),
        [&](const Definition& existing) {
            return existing.name == definition.name;
        });
    if (duplicate != equations_.end()) {
        throw std::runtime_error(
            "Equation system already contains '" + definition.name + "'.");
    }
    equations_.push_back(std::move(definition));
}

const Definition& System::at(const std::string& name) const {
    const auto found = std::find_if(
        equations_.begin(), equations_.end(),
        [&](const Definition& equation) { return equation.name == name; });
    if (found == equations_.end()) {
        throw std::runtime_error(
            "Equation system does not contain '" + name + "'.");
    }
    return *found;
}

void System::addRightTerm(const std::string& name, Term term) {
    auto found = std::find_if(
        equations_.begin(),equations_.end(),
        [&](const Definition& equation) { return equation.name == name; });
    if (found == equations_.end()) {
        throw std::runtime_error(
            "Equation system does not contain '"+name+"'.");
    }
    if (term.primary.name.empty()) {
        throw std::runtime_error(
            "Equation contribution for '"+name+"' has no symbol.");
    }
    found->right.terms.push_back(std::move(term));
    found->validate();
}

bool System::contains(TermKind kind) const {
    for (const Definition& equation : equations_) {
        const auto containsKind = [kind](const Expression& expression) {
            return std::any_of(
                expression.terms.begin(), expression.terms.end(),
                [kind](const Term& term) { return term.kind == kind; });
        };
        if (containsKind(equation.left) || containsKind(equation.right)) {
            return true;
        }
    }
    return false;
}

} // namespace SF::Equation
