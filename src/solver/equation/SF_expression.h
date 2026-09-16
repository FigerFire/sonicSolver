#pragma once

/// @file SF_expression.h
/// @brief 方程层的强类型符号表达式；只描述方程，不选择数值格式。

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SF::Equation {

/// @brief 方程中可被离散层绑定的物理量符号。
struct Symbol {
    std::string name;
};

/// @brief 方程项的数学类别。
enum class TermKind {
    Transient,
    Divergence,
    Gradient,
    Diffusion,
    ExplicitSource,
    ImplicitSource,
    Constraint,
    AlgebraicRelation
};

enum class EvaluationMode { Explicit, Implicit, Algebraic };
enum class TimeLevel { Current, Previous, Stage };

/// @brief 一个尚未选择数值格式的方程项。
struct Term {
    TermKind kind = TermKind::Transient;
    Symbol primary;
    Symbol secondary;
    EvaluationMode mode = EvaluationMode::Explicit;
    TimeLevel timeLevel = TimeLevel::Current;
    std::string closureProvider;
    std::string boundaryRequirement;
    std::string linearizationRequirement;
};

/// @brief 方程一侧的有序项列表。
struct Expression {
    std::vector<Term> terms;
};

/// @brief 完整的符号方程定义。
struct Definition {
    std::string name;
    Expression left;
    Expression right;

    /// @brief 验证项名称有效且最多含一个瞬态项。
    void validate() const;
};

/// @brief 构造瞬态项 `ddt(q)`。
inline Term ddt(Symbol q) {
    return {TermKind::Transient, std::move(q), {}};
}

/// @brief 构造守恒通量散度项 `div(flux)`。
inline Term div(Symbol flux) {
    return {TermKind::Divergence, std::move(flux), {}};
}

/// @brief 构造梯度项 `grad(q)`。
inline Term gradient(Symbol q) {
    return {TermKind::Gradient, std::move(q), {}};
}

/// @brief 构造扩散项 `diffusion(coefficient, unknown)`。
inline Term diffusion(Symbol coefficient, Symbol unknown) {
    return {TermKind::Diffusion, std::move(unknown), std::move(coefficient)};
}

/// @brief 构造显式源项。
inline Term source(Symbol value) {
    return {TermKind::ExplicitSource, std::move(value), {}};
}

/// @brief 构造隐式源项。
inline Term implicitSource(Symbol value, Symbol unknown) {
    Term term{TermKind::ImplicitSource, std::move(unknown), std::move(value)};
    term.mode = EvaluationMode::Implicit;
    return term;
}

/// @brief 构造边界或代数约束项；约束由边界/离散层执行，不伪装成体源。
inline Term constraint(Symbol value) {
    Term term{TermKind::Constraint, std::move(value), {}};
    term.mode = EvaluationMode::Algebraic;
    return term;
}

/// @brief 构造不经过 PDE 离散的代数关系。
inline Term algebraic(Symbol value) {
    Term term{TermKind::AlgebraicRelation, std::move(value), {}};
    term.mode = EvaluationMode::Algebraic;
    return term;
}

inline Expression operator+(Term left, Term right) {
    return {{std::move(left), std::move(right)}};
}

inline Expression operator+(Expression expression, Term term) {
    expression.terms.push_back(std::move(term));
    return expression;
}

inline Expression operator+(Term term, Expression expression) {
    expression.terms.insert(expression.terms.begin(), std::move(term));
    return expression;
}

inline Expression operator+(Expression left, Expression right) {
    left.terms.insert(left.terms.end(),
                      std::make_move_iterator(right.terms.begin()),
                      std::make_move_iterator(right.terms.end()));
    return left;
}

/// @brief 生成未命名定义；由 `named` 在注册前赋予稳定名称。
inline Definition operator==(Expression left, Term right) {
    return {{}, std::move(left), {{std::move(right)}}};
}

inline Definition operator==(Expression left, Expression right) {
    return {{}, std::move(left), std::move(right)};
}

inline Definition operator==(Term left, Term right) {
    return {{}, {{std::move(left)}}, {{std::move(right)}}};
}

/// @brief 右端裸物理量按显式源解释，使方程可写为 `... == mrf`。
inline Definition operator==(Expression left, Symbol right) {
    return std::move(left) == source(std::move(right));
}

inline Definition operator==(Term left, Symbol right) {
    return std::move(left) == source(std::move(right));
}

/// @brief 给表达式赋予用于配置、诊断和日志的稳定方程名。
inline Definition named(std::string name, Definition definition) {
    definition.name = std::move(name);
    definition.validate();
    return definition;
}

/// @brief 一组由 Solver Algorithm 按耦合顺序推进的方程定义。
class System {
public:
    void add(Definition definition);
    const std::vector<Definition>& equations() const { return equations_; }
    const Definition& at(const std::string& name) const;
    /// @brief 向已注册方程追加一个有序右端数学贡献。
    void addRightTerm(const std::string& name, Term term);
    /// @brief 判断系统任一方程是否含指定数学项。
    bool contains(TermKind kind) const;

private:
    std::vector<Definition> equations_;
};

} // namespace SF::Equation
