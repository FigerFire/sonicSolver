#pragma once

/// @file SF_formula.h
/// @brief Mathematical expression AST; EquationRegistry is the WHAT authority.

#include "core/system/SF_provenance.h"
#include "core/system/SF_expression.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

/// @brief A node in a mathematical formula. Operator names describe mathematics,
/// not a numerical scheme or an execution mode.
struct FormulaExpr {
    enum class Kind { Symbol, Constant, Add, Subtract, Multiply, Divide, Negate, Operator };
    Kind kind = Kind::Symbol;
    std::string name;
    double constant = 0.0;
    std::vector<FormulaExpr> arguments;
    /// @brief Stable address for a discretized operator occurrence within a formula.
    std::string occurrence;

    static FormulaExpr symbol(std::string name);
    static FormulaExpr constantValue(double value);
    static FormulaExpr add(FormulaExpr left, FormulaExpr right);
    static FormulaExpr subtract(FormulaExpr left, FormulaExpr right);
    static FormulaExpr multiply(FormulaExpr left, FormulaExpr right);
    static FormulaExpr divide(FormulaExpr left, FormulaExpr right);
    static FormulaExpr negate(FormulaExpr value);
    static FormulaExpr op(std::string name, std::vector<FormulaExpr> arguments,
                          std::string occurrence = {});
};

/// @brief Stable equation identity and its mathematical AST. Targets and order belong to HOW.
struct Equation {
    std::string id;
    FormulaExpr lhs;
    FormulaExpr rhs;
    Provenance origin;
    /// Authored AST is the mathematical authority; a retained legacy
    /// SF::Equation::Definition may describe it but must not regenerate it.
    bool authored = false;
};

/// @brief Explicit add/replace/disable semantics for built-in and user formulas.
class EquationRegistry {
public:
    void add(Equation formula);
    void replace(Equation formula);
    void disable(std::string_view id);
    const Equation& at(std::string_view id) const;
    bool contains(std::string_view id) const {
        for (const auto& value:formulas_) if (value.id==id) return true;
        return false;
    }
    const std::vector<Equation>& entries() const { return formulas_; }
private:
    std::vector<Equation> formulas_;
};

/// @brief Structural equality key. Origin and display labels are excluded.
std::string canonicalFormula(const Equation& formula);
std::string formulaText(const Equation& formula);
bool dependsOn(const FormulaExpr& expression, std::string_view symbol);
/// @brief Transitional conversion of existing equation definitions into the
/// same mathematical AST consumed by FormulaCall compilation.
Equation formulaFromEquation(const SF::Equation::Definition& equation,
                            Provenance origin = {});

} // namespace SF::System
