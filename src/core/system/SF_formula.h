#pragma once

/// @file SF_formula.h
/// @brief Mathematical expression AST; EquationRegistry is the WHAT authority.

#include "core/system/SF_provenance.h"
#include "methods/math/SF_signature.h"
#include <map>

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

/// @brief A node in a mathematical formula. Operator names describe mathematics,
/// not a numerical scheme or an execution mode.
class FormulaExpr {
public:
    enum class Kind { Symbol, Constant, Add, Subtract, Multiply, Divide, Negate, Operator };
    FormulaExpr() = default; // A real literal zero; never an untyped symbol.
    Kind kind() const {return kind_;}
    const std::string& name() const {return name_;}
    double constant() const {return constant_;}
    const std::vector<FormulaExpr>& arguments() const {return arguments_;}
    const std::string& occurrence() const {return occurrence_;}
    const Math::Signature& signature() const {return signature_;}
    Math::Op operation() const {return operation_;}
    std::uint64_t extensionIdentity() const {return extensionIdentity_;}
    auto extensionRule() const {return extensionResult_;}
    FormulaExpr atOccurrence(std::string value) const {auto copy=*this;copy.occurrence_=std::move(value);return copy;}
    static FormulaExpr symbol(std::string name, Math::Signature signature);
    static FormulaExpr constantValue(double value);
    static FormulaExpr constantValue(double value, Math::Signature signature);
    static FormulaExpr zeroLike(const FormulaExpr& value);
    static FormulaExpr add(FormulaExpr left, FormulaExpr right);
    static FormulaExpr subtract(FormulaExpr left, FormulaExpr right);
    static FormulaExpr multiply(FormulaExpr left, FormulaExpr right);
    static FormulaExpr divide(FormulaExpr left, FormulaExpr right);
    static FormulaExpr negate(FormulaExpr value);
    static FormulaExpr op(Math::Op operation, std::vector<FormulaExpr> arguments, std::string occurrence = {});
    static FormulaExpr op(std::string name, std::vector<FormulaExpr> arguments, std::string occurrence = {});
    static FormulaExpr op(const Math::OperatorRule& rule, std::vector<FormulaExpr> arguments, std::string occurrence = {});
    static FormulaExpr power(FormulaExpr value,double exponent);
    /// @brief 重新构造并重新推导签名；原表达式不变。
    FormulaExpr withArguments(std::vector<FormulaExpr> values) const;
    friend FormulaExpr operator-(FormulaExpr a) {return negate(std::move(a));}
    friend FormulaExpr operator+(FormulaExpr a,FormulaExpr b) {return add(std::move(a),std::move(b));}
    friend FormulaExpr operator-(FormulaExpr a,FormulaExpr b) {return subtract(std::move(a),std::move(b));}
    friend FormulaExpr operator*(FormulaExpr a,FormulaExpr b) {return multiply(std::move(a),std::move(b));}
    friend FormulaExpr operator/(FormulaExpr a,FormulaExpr b) {return divide(std::move(a),std::move(b));}
private:
    friend struct Equation;
    static FormulaExpr make(Math::Op operation,std::vector<FormulaExpr> arguments,std::string occurrence={});
    Kind kind_=Kind::Constant;
    std::string name_;
    double constant_=0;
    std::vector<FormulaExpr> arguments_;
    std::string occurrence_;
    Math::Signature signature_=Math::signature(0,{});
    Math::Op operation_=Math::Op::Add;
    bool literalZero_=true;
    std::uint64_t extensionIdentity_=0;
    Math::Signature (*extensionResult_)(const std::vector<Math::Signature>&)=nullptr;
};

/// @brief 一个贡献中的已声明数学符号；纯逻辑量不要求物理数组。
class FormulaSymbols {
public:
    void add(std::string name,Math::Signature signature) {
        Math::requireKnown(signature);
        if(!symbols_.emplace(std::move(name),signature).second)throw std::runtime_error("Duplicate mathematical symbol");
    }
    FormulaExpr at(std::string_view name) const {
        const auto at=symbols_.find(std::string(name));
        if(at==symbols_.end())throw std::runtime_error("Unknown declared mathematical symbol: "+std::string(name));
        return FormulaExpr::symbol(at->first,at->second);
    }
private:
    std::map<std::string,Math::Signature> symbols_;
};

/// @brief Stable equation identity and its mathematical AST. Targets and order belong to HOW.
struct Equation {
    Equation()=default;
    Equation(std::string id,FormulaExpr lhs,FormulaExpr rhs,Provenance origin={});
    void validateTypes() const;
private:
    friend class EquationRegistry;
    void contextualizeZero();
public:
    std::string id;
    FormulaExpr lhs;
    FormulaExpr rhs;
    Provenance origin;
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
} // namespace SF::System
