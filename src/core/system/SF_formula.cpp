/// @file SF_formula.cpp
/// @brief Equation AST construction and source-independent composition.

#include "SF_formula.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace SF::System {
namespace {

FormulaExpr binary(FormulaExpr::Kind kind, FormulaExpr left, FormulaExpr right) {
    return {kind,{},0.0,{std::move(left),std::move(right)},{}};
}

void validateExpr(const FormulaExpr& expression) {
    using Kind = FormulaExpr::Kind;
    const auto size = expression.arguments.size();
    switch (expression.kind) {
        case Kind::Symbol:
            if (expression.name.empty() || size) break;
            return;
        case Kind::Constant:
            if (std::isfinite(expression.constant) && !size) return;
            break;
        case Kind::Negate:
            if (size == 1) break;
            throw std::runtime_error("Equation negation requires one operand.");
        case Kind::Add:
        case Kind::Subtract:
        case Kind::Multiply:
        case Kind::Divide:
            if (size == 2) break;
            throw std::runtime_error("Equation binary operator requires two operands.");
        case Kind::Operator:
            if (!expression.name.empty() && size) break;
            throw std::runtime_error("Equation mathematical operator requires a name and operands.");
    }
    if (expression.kind == Kind::Symbol || expression.kind == Kind::Constant)
        throw std::runtime_error("Equation symbol or constant has invalid data.");
    for (const auto& argument : expression.arguments) validateExpr(argument);
}

void validateOccurrences(const FormulaExpr& expression,
                         std::unordered_set<std::string>& paths) {
    if (expression.kind==FormulaExpr::Kind::Operator
        && !expression.occurrence.empty()
        && !paths.insert(expression.occurrence).second)
        throw std::runtime_error("Equation has duplicate operator occurrence '"
                                 +expression.occurrence+"'.");
    for (const auto& child:expression.arguments)
        validateOccurrences(child,paths);
}

void validateFormula(const Equation& formula) {
    validateExpr(formula.lhs);
    validateExpr(formula.rhs);
    std::unordered_set<std::string> paths;
    validateOccurrences(formula.lhs,paths);
    validateOccurrences(formula.rhs,paths);
}

void encode(std::ostringstream& output, const FormulaExpr& expression) {
    output << static_cast<int>(expression.kind) << '(';
    output << expression.name.size() << ':' << expression.name;
    if (expression.kind == FormulaExpr::Kind::Constant)
        output << std::setprecision(17) << expression.constant;
    output << '[';
    for (const auto& argument : expression.arguments) encode(output,argument);
    output << "])";
}

std::string textOf(const FormulaExpr& expression) {
    using Kind=FormulaExpr::Kind;
    if (expression.kind==Kind::Symbol) return expression.name;
    if (expression.kind==Kind::Constant) {
        std::ostringstream output;
        output << std::setprecision(17) << expression.constant;
        return output.str();
    }
    if (expression.kind==Kind::Negate)
        return "-("+textOf(expression.arguments[0])+")";
    if (expression.kind==Kind::Operator) {
        std::string result=expression.name+"(";
        for (std::size_t i=0;i<expression.arguments.size();++i) {
            if (i) result+=", ";
            result+=textOf(expression.arguments[i]);
        }
        return result+")";
    }
    const char* operation=expression.kind==Kind::Add?" + "
        : expression.kind==Kind::Subtract?" - "
        : expression.kind==Kind::Multiply?" * ":" / ";
    return "("+textOf(expression.arguments[0])+operation
        +textOf(expression.arguments[1])+")";
}

FormulaExpr termExpr(const SF::Equation::Term& term, std::size_t ordinal) {
    const auto symbol = term.primary.name == "zero"
        ? FormulaExpr::constantValue(0.0)
        : FormulaExpr::symbol(term.primary.name);
    const std::string occurrence = "term-" + std::to_string(ordinal);
    switch (term.kind) {
        case SF::Equation::TermKind::Transient:
            return FormulaExpr::op("ddt",{symbol},occurrence);
        case SF::Equation::TermKind::Divergence:
            return FormulaExpr::op("div",{symbol},occurrence);
        case SF::Equation::TermKind::Gradient:
            return FormulaExpr::op("grad",{symbol},occurrence);
        case SF::Equation::TermKind::Diffusion:
            return FormulaExpr::op("diffusion",
                {FormulaExpr::symbol(term.secondary.name),symbol},occurrence);
        case SF::Equation::TermKind::Source:
            return symbol;
        case SF::Equation::TermKind::Constraint:
            return FormulaExpr::op("constraint",{symbol},occurrence);
        case SF::Equation::TermKind::AlgebraicRelation:
            return FormulaExpr::op("algebraic",{symbol},occurrence);
    }
    throw std::runtime_error("Unknown equation term kind.");
}

FormulaExpr sideExpr(const SF::Equation::Expression& side, std::size_t& ordinal) {
    if (side.terms.empty()) return FormulaExpr::constantValue(0.0);
    FormulaExpr result = termExpr(side.terms.front(),ordinal++);
    for (std::size_t i=1;i<side.terms.size();++i)
        result = FormulaExpr::add(std::move(result),
                                  termExpr(side.terms[i],ordinal++));
    return result;
}

} // namespace

FormulaExpr FormulaExpr::symbol(std::string name) {
    return {Kind::Symbol,std::move(name),0.0,{},{}};
}
FormulaExpr FormulaExpr::constantValue(double value) {
    return {Kind::Constant,{},value,{},{}};
}
FormulaExpr FormulaExpr::add(FormulaExpr left,FormulaExpr right) {
    return binary(Kind::Add,std::move(left),std::move(right));
}
FormulaExpr FormulaExpr::subtract(FormulaExpr left,FormulaExpr right) {
    return binary(Kind::Subtract,std::move(left),std::move(right));
}
FormulaExpr FormulaExpr::multiply(FormulaExpr left,FormulaExpr right) {
    return binary(Kind::Multiply,std::move(left),std::move(right));
}
FormulaExpr FormulaExpr::divide(FormulaExpr left,FormulaExpr right) {
    return binary(Kind::Divide,std::move(left),std::move(right));
}
FormulaExpr FormulaExpr::negate(FormulaExpr value) {
    return {Kind::Negate,{},0.0,{std::move(value)},{}};
}
FormulaExpr FormulaExpr::op(std::string name,
                            std::vector<FormulaExpr> arguments,
                            std::string occurrence) {
    return {Kind::Operator,std::move(name),0.0,std::move(arguments),
            std::move(occurrence)};
}

void EquationRegistry::add(Equation formula) {
    if (formula.id.empty()) throw std::runtime_error("Equation ID is empty.");
    validateFormula(formula);
    const auto found=std::find_if(formulas_.begin(),formulas_.end(),
        [&](const Equation& current) { return current.id==formula.id; });
    if (found!=formulas_.end())
        throw std::runtime_error("Duplicate Equation ID '"+formula.id+"'.");
    formulas_.push_back(std::move(formula));
}

void EquationRegistry::replace(Equation formula) {
    validateFormula(formula);
    const auto found=std::find_if(formulas_.begin(),formulas_.end(),
        [&](const Equation& current) { return current.id==formula.id; });
    if (found==formulas_.end())
        throw std::runtime_error("Cannot replace missing Equation '"+formula.id+"'.");
    *found=std::move(formula);
}

void EquationRegistry::disable(std::string_view id) {
    const auto found=std::find_if(formulas_.begin(),formulas_.end(),
        [&](const Equation& current) { return current.id==id; });
    if (found==formulas_.end())
        throw std::runtime_error("Cannot disable missing Equation '"+std::string(id)+"'.");
    formulas_.erase(found);
}

const Equation& EquationRegistry::at(std::string_view id) const {
    const auto found=std::find_if(formulas_.begin(),formulas_.end(),
        [&](const Equation& current) { return current.id==id; });
    if (found==formulas_.end())
        throw std::runtime_error("Equation '"+std::string(id)+"' is not registered.");
    return *found;
}

std::string canonicalFormula(const Equation& formula) {
    std::ostringstream output;
    encode(output,formula.lhs);
    output << '=';
    encode(output,formula.rhs);
    return output.str();
}

std::string formulaText(const Equation& formula) {
    return textOf(formula.lhs)+" = "+textOf(formula.rhs);
}

bool dependsOn(const FormulaExpr& expression,std::string_view symbol) {
    if (expression.kind==FormulaExpr::Kind::Symbol
        && expression.name==symbol) return true;
    return std::any_of(expression.arguments.begin(),expression.arguments.end(),
        [&](const FormulaExpr& argument) { return dependsOn(argument,symbol); });
}

Equation formulaFromEquation(const SF::Equation::Definition& equation,
                            Provenance origin) {
    std::size_t ordinal=0;
    return {equation.name,sideExpr(equation.left,ordinal),
            sideExpr(equation.right,ordinal),std::move(origin)};
}

} // namespace SF::System
