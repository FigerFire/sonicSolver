/// @file SF_formulaCompiler.cpp
/// @brief Structural AST compilation; numerical mathematics belongs to providers.

#include "SF_formulaCompiler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::System {
namespace {

using Kind = FormulaExpr::Kind;

void appendProvider(std::vector<std::string>& used,const std::string& id) {
    if (std::find(used.begin(),used.end(),id)==used.end()) used.push_back(id);
}

FormulaValueKernel valueKernel(
        const Equation& formula,const FormulaExpr& expression,
        const FormulaOperatorCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,
        std::vector<std::string>& used,const std::vector<std::string>* slots=nullptr) {
    switch (expression.kind) {
        case Kind::Constant: {
            const double value=expression.constant;
            return [value](const FormulaValues&,int,int) { return value; };
        }
        case Kind::Symbol: {
            const auto name=expression.name;
            if (slots) {
                const auto found=std::find(slots->begin(),slots->end(),name);
                if (found==slots->end()) throw std::runtime_error("Unbound AST value symbol: "+name);
                const auto slot=static_cast<std::size_t>(found-slots->begin());
                return [slot](const FormulaValues& values,int cell,int component) {
                    return values.boundReads.at(slot)(cell,component);
                };
            }
            return [name](const FormulaValues& values,int cell,int component) {
                return values.read(name,cell,component);
            };
        }
        case Kind::Negate: {
            auto child=valueKernel(formula,expression.arguments[0],catalog,bindings,used,slots);
            return [child=std::move(child)](const FormulaValues& values,int cell,int component) {
                return -child(values,cell,component);
            };
        }
        case Kind::Add:
        case Kind::Subtract:
        case Kind::Multiply:
        case Kind::Divide: {
            auto left=valueKernel(formula,expression.arguments[0],catalog,bindings,used,slots);
            auto right=valueKernel(formula,expression.arguments[1],catalog,bindings,used,slots);
            const auto kind=expression.kind;
            return [left=std::move(left),right=std::move(right),kind](
                    const FormulaValues& values,int cell,int component) {
                const double a=left(values,cell,component);
                const double b=right(values,cell,component);
                if (kind==Kind::Add) return a+b;
                if (kind==Kind::Subtract) return a-b;
                if (kind==Kind::Multiply) return a*b;
                if (b==0.0) throw std::runtime_error("Equation division by zero.");
                return a/b;
            };
        }
        case Kind::Operator: {
            const auto& provider=catalog.resolve(formula,expression,bindings);
            if (!provider.compileValue)
                throw std::runtime_error("Equation operator '"+expression.name
                    +"' provider '"+provider.id+"' cannot evaluate explicitly.");
            appendProvider(used,provider.id);
            return provider.compileValue(expression);
        }
    }
    throw std::runtime_error("Unknown Equation expression kind.");
}

FormulaLinearTerm scale(FormulaLinearTerm term,double value) {
    for (auto& coefficient:term.coefficients) coefficient.value*=value;
    term.constant*=value;
    return term;
}

FormulaLinearTerm combine(FormulaLinearTerm left,
                          const FormulaLinearTerm& right,double sign) {
    for (const auto& coefficient:right.coefficients)
        left.add(coefficient.column,sign*coefficient.value);
    left.constant+=sign*right.constant;
    return left;
}

FormulaLinearKernel linearKernel(
        const Equation& formula,const FormulaExpr& expression,
        const std::string& target,const FormulaOperatorCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,
        std::vector<std::string>& used) {
    if (!dependsOn(expression,target)) {
        auto known=valueKernel(formula,expression,catalog,bindings,used);
        return [known=std::move(known)](
                const FormulaValues& values,int cell,int component) {
            return FormulaLinearTerm{{},known(values,cell,component)};
        };
    }
    switch (expression.kind) {
        case Kind::Symbol:
            return [target](const FormulaValues& values,int cell,int component) {
                FormulaLinearTerm result;
                result.add(values.dof(target,cell,component),1.0);
                return result;
            };
        case Kind::Negate: {
            auto child=linearKernel(formula,expression.arguments[0],target,
                                    catalog,bindings,used);
            return [child=std::move(child)](
                    const FormulaValues& values,int cell,int component) {
                return scale(child(values,cell,component),-1.0);
            };
        }
        case Kind::Add:
        case Kind::Subtract: {
            auto left=linearKernel(formula,expression.arguments[0],target,
                                   catalog,bindings,used);
            auto right=linearKernel(formula,expression.arguments[1],target,
                                    catalog,bindings,used);
            const double sign=expression.kind==Kind::Add?1.0:-1.0;
            return [left=std::move(left),right=std::move(right),sign](
                    const FormulaValues& values,int cell,int component) {
                return combine(left(values,cell,component),
                               right(values,cell,component),sign);
            };
        }
        case Kind::Multiply:
        case Kind::Divide: {
            const bool leftDepends=dependsOn(expression.arguments[0],target);
            const bool rightDepends=dependsOn(expression.arguments[1],target);
            if ((leftDepends && rightDepends)
                || (expression.kind==Kind::Divide && rightDepends))
                throw std::runtime_error("Implicit Equation contains a nonlinear or "
                    "target-dependent denominator; declare a linearized Equation.");
            const std::size_t linearIndex=rightDepends?1:0;
            const std::size_t factorIndex=rightDepends?0:1;
            auto factor=valueKernel(formula,expression.arguments[factorIndex],
                                    catalog,bindings,used);
            auto linear=linearKernel(formula,expression.arguments[linearIndex],
                                     target,catalog,bindings,used);
            const bool divide=expression.kind==Kind::Divide;
            return [linear=std::move(linear),factor=std::move(factor),divide](
                    const FormulaValues& values,int cell,int component) {
                const double value=factor(values,cell,component);
                if (divide && value==0.0)
                    throw std::runtime_error("Implicit Equation divides by zero.");
                return scale(linear(values,cell,component),
                             divide?1.0/value:value);
            };
        }
        case Kind::Operator: {
            const auto& provider=catalog.resolve(formula,expression,bindings);
            if (!provider.compileLinear)
                throw std::runtime_error("Implicit Equation operator '"
                    +expression.name+"' provider '"+provider.id
                    +"' has no linear assembly capability.");
            appendProvider(used,provider.id);
            return provider.compileLinear(expression,target);
        }
        case Kind::Constant:
            break;
    }
    throw std::runtime_error("Implicit Equation cannot lower this AST node.");
}

bool isTargetDdt(const FormulaExpr& expression,const std::string& target) {
    return expression.kind==Kind::Operator && expression.name=="ddt"
        && expression.arguments.size()==1
        && expression.arguments[0].kind==Kind::Symbol
        && expression.arguments[0].name==target;
}

bool hasOccurrence(const FormulaExpr& expression,const std::string& path) {
    if (expression.kind==Kind::Operator && expression.occurrence==path)
        return true;
    return std::any_of(expression.arguments.begin(),expression.arguments.end(),
        [&](const FormulaExpr& child) { return hasOccurrence(child,path); });
}

FormulaExpr withoutDdt(const FormulaExpr& expression,
                       const std::string& target,int& count) {
    if (isTargetDdt(expression,target)) {
        ++count;
        return FormulaExpr::constantValue(0.0);
    }
    if (expression.kind==Kind::Add && expression.arguments.size()==2) {
        return FormulaExpr::add(
            withoutDdt(expression.arguments[0],target,count),
            withoutDdt(expression.arguments[1],target,count));
    }
    return expression;
}

} // namespace

FormulaValueKernel compileFormulaValue(const Equation& equation,const FormulaExpr& expression,
        const FormulaOperatorCatalog& providers,std::vector<std::string>& used,
        const std::vector<std::string>& boundSymbols) {
    return valueKernel(equation,expression,providers,{},used,boundSymbols.empty()?nullptr:&boundSymbols);
}

void FormulaOperatorCatalog::add(FormulaOperatorProvider provider) {
    if (provider.id.empty() || provider.mathematicalOperator.empty()
        || !provider.matches)
        throw std::runtime_error("Equation numerical provider is incomplete.");
    const auto found=std::find_if(providers_.begin(),providers_.end(),
        [&](const auto& value) { return value.id==provider.id; });
    if (found!=providers_.end())
        throw std::runtime_error("Duplicate Equation numerical provider '"
            +provider.id+"'.");
    providers_.push_back(std::move(provider));
}

const FormulaOperatorProvider& FormulaOperatorCatalog::resolve(
        const Equation& formula,const FormulaExpr& expression,
        const std::vector<FormulaOperatorBinding>& bindings) const {
    const auto* selected=selectFormulaBinding(
        formula,expression,expression.occurrence,bindings);
    const FormulaOperatorProvider* resolved=nullptr;
    for (const auto& provider:providers_) {
        if (provider.mathematicalOperator!=expression.name
            || !provider.matches(expression)
            || (selected && selected->provider!=provider.id))
            continue;
        if (resolved)
            throw std::runtime_error("Equation operator '"+formula.id+"/"
                +expression.occurrence+"' matches multiple numerical providers.");
        resolved=&provider;
    }
    if (!resolved)
        throw std::runtime_error("Equation operator '"+formula.id+"/"
            +expression.occurrence+"' has no compatible numerical provider"
            +(selected ? " for binding '"+selected->provider+"'" : "")+".");
    return *resolved;
}

void CompiledFormulaCall::writeSolution(
        FormulaValues& values,const LinearAlgebra::SolveResult& result) const {
    if (call.mode!=Legacy::FormulaMode::Implicit || !assemble)
        throw std::runtime_error("Equation call is not an implicit assembly.");
    const auto expected=static_cast<std::size_t>(values.ownedCells)
        *static_cast<std::size_t>(values.components);
    if (!result.converged || result.solution.size()!=expected)
        throw std::runtime_error("Implicit Equation solve failed or returned a "
                                 "solution with the wrong layout.");
    for (int cell=0;cell<values.ownedCells;++cell)
        for (int component=0;component<values.components;++component) {
            const double value=result.solution[static_cast<std::size_t>(cell)
                *values.components+component];
            if (!std::isfinite(value))
                throw std::runtime_error("Implicit Equation solution is non-finite.");
            values.write(call.target,cell,component,value);
        }
}

CompiledFormulaCall FormulaCompiler::compile(
        const EquationRegistry& formulas,const Legacy::FormulaCall& request,
        const FormulaOperatorCatalog& providers,
        const std::vector<FormulaOperatorBinding>& bindings) {
    const Equation& formula=formulas.at(request.equation);
    for (const auto& binding:bindings) {
        if (binding.formula==formula.id && !binding.occurrence.empty()
            && !hasOccurrence(formula.lhs,binding.occurrence)
            && !hasOccurrence(formula.rhs,binding.occurrence))
            throw std::runtime_error("Equation operator binding '"+formula.id
                +"/"+binding.occurrence+"' does not address an AST operator.");
    }
    CompiledFormulaCall compiled;
    compiled.call=request;
    if (compiled.call.mode==Legacy::FormulaMode::Assign) {
        if (formula.lhs.kind!=Kind::Symbol
            || (!request.target.empty() && request.target!=formula.lhs.name))
            throw std::runtime_error("Assignment Equation requires a direct "
                                     "left-hand target symbol.");
        compiled.call.target=formula.lhs.name;
        auto rhs=valueKernel(formula,formula.rhs,providers,bindings,
                             compiled.numericalProviders);
        const auto target=compiled.call.target;
        compiled.assign=[rhs=std::move(rhs),target](FormulaValues& values) {
            if (values.ownedCells<0 || values.components<1 || !values.read
                || !values.write)
                throw std::runtime_error("Assignment Equation has invalid storage bindings.");
            std::vector<double> updated(static_cast<std::size_t>(values.ownedCells)
                *values.components);
            for (int cell=0;cell<values.ownedCells;++cell)
                for (int component=0;component<values.components;++component) {
                    const double value=rhs(values,cell,component);
                    if (!std::isfinite(value))
                        throw std::runtime_error("Assignment Equation produced a non-finite value.");
                    updated[static_cast<std::size_t>(cell)*values.components+component]=value;
                }
            for (int cell=0;cell<values.ownedCells;++cell)
                for (int component=0;component<values.components;++component)
                    values.write(target,cell,component,
                        updated[static_cast<std::size_t>(cell)
                            *values.components+component]);
        };
        compiled.backend="formula.assignment";
        return compiled;
    }
    if (request.target.empty())
        throw std::runtime_error("Explicit/implicit FormulaCall has no target.");
    if (request.mode==Legacy::FormulaMode::Explicit) {
        int count=0;
        const auto rest=withoutDdt(formula.lhs,request.target,count);
        const auto hasDdt=[&](const auto& self,const FormulaExpr& expr) -> bool {
            if (expr.kind==Kind::Operator && expr.name=="ddt") return true;
            return std::any_of(expr.arguments.begin(),expr.arguments.end(),
                [&](const FormulaExpr& child) { return self(self,child); });
        };
        if (count!=1 || hasDdt(hasDdt,formula.rhs))
            throw std::runtime_error("Explicit Equation requires exactly one "
                                     "positive ddt(target) on the left.");
        auto left=valueKernel(formula,rest,providers,bindings,
                              compiled.numericalProviders);
        auto right=valueKernel(formula,formula.rhs,providers,bindings,
                               compiled.numericalProviders);
        compiled.evaluateRhs=[left=std::move(left),right=std::move(right)](
                const FormulaValues& values,int cell,int component) {
            return right(values,cell,component)
                -left(values,cell,component);
        };
        compiled.backend="formula.explicit-residual";
        return compiled;
    }
    if (!dependsOn(formula.lhs,request.target)
        && !dependsOn(formula.rhs,request.target))
        throw std::runtime_error("Implicit Equation does not contain target '"
                                 +request.target+"'.");
    auto left=linearKernel(formula,formula.lhs,request.target,
                           providers,bindings,compiled.numericalProviders);
    auto right=linearKernel(formula,formula.rhs,request.target,
                            providers,bindings,compiled.numericalProviders);
    const auto target=request.target;
    compiled.assemble=[left=std::move(left),right=std::move(right),target](
            const FormulaValues& values) {
        if (values.ownedCells<1 || values.components<1 || !values.read
            || !values.dof)
            throw std::runtime_error("Implicit Equation has invalid state/DOF bindings.");
        LinearAlgebra::GlobalDofSystem system;
        system.rows.reserve(static_cast<std::size_t>(values.ownedCells)
                            *values.components);
        for (int cell=0;cell<values.ownedCells;++cell)
            for (int component=0;component<values.components;++component) {
                const auto lhs=left(values,cell,component);
                const auto rhs=right(values,cell,component);
                LinearAlgebra::GlobalDofRow row(values.dof(target,cell,component));
                for (const auto& item:lhs.coefficients)
                    row.add(item.column,item.value);
                for (const auto& item:rhs.coefficients)
                    row.add(item.column,-item.value);
                const double rightHandSide=rhs.constant-lhs.constant;
                if (!std::isfinite(rightHandSide))
                    throw std::runtime_error("Implicit Equation RHS is non-finite.");
                row.setRightHandSide(rightHandSide);
                row.setInitialGuess(values.read(target,cell,component));
                system.rows.push_back(std::move(row));
            }
        return system;
    };
    compiled.backend="formula.linear-assembly";
    return compiled;
}

} // namespace SF::System
