/// @file SF_formulaCompiler.cpp
/// @brief Structural AST compilation; numerical mathematics belongs to providers.

#include "SF_formulaCompiler.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"

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
        std::vector<std::string>& used,const std::vector<std::string>* slots,const StateRegistry* state,
        std::vector<std::string>* trace,const std::string& path) {
    const auto record=[&](std::string backend){if(trace)trace->push_back(path+" ["+expression.occurrence()+"] -> "+std::move(backend));};
    switch (expression.kind()) {
        case Kind::Constant: {
            record("local.constant");
            const double value=expression.constant();
            return [value](const FormulaValues&,int,int) { return value; };
        }
        case Kind::Symbol: {
            const auto name=expression.name();
            if(state && state->contains(name) && state->at(name).constantValue) {
                record("fixed coefficient "+name);const double value=*state->at(name).constantValue;
                if(!std::isfinite(value))throw std::runtime_error("Invalid: non-finite fixed coefficient "+name);
                return [value](const FormulaValues&,int,int){return value;};
            }
            record("bound symbol "+name);
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
            record("local.negate");
            auto child=valueKernel(formula,expression.arguments()[0],catalog,bindings,used,slots,state,trace,path+"/0");
            return [child=std::move(child)](const FormulaValues& values,int cell,int component) {
                return -child(values,cell,component);
            };
        }
        case Kind::Add:
        case Kind::Subtract:
        case Kind::Multiply:
        case Kind::Divide: {
            auto left=valueKernel(formula,expression.arguments()[0],catalog,bindings,used,slots,state,trace,path+"/0");
            auto right=valueKernel(formula,expression.arguments()[1],catalog,bindings,used,slots,state,trace,path+"/1");
            record(std::string("local.")+Math::name(expression.operation()));
            const auto kind=expression.kind();
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
            if(expression.operation()==Math::Op::Pow) {
                record("local.pow");
                auto base=valueKernel(formula,expression.arguments()[0],catalog,bindings,used,slots,state,trace,path+"/0");
                const double exponent=expression.arguments()[1].constant();
                return [base=std::move(base),exponent](const FormulaValues& values,int cell,int component){return std::pow(base(values,cell,component),exponent);};
            }
            const auto& provider=catalog.resolve(formula,expression,bindings);
            if (!provider.compileValue)
                throw std::runtime_error("Equation operator '"+expression.name()
                    +"' provider '"+provider.id+"' cannot evaluate explicitly.");
            record(provider.id);
            appendProvider(used,provider.id);
            return provider.compileValue(expression);
        }
    }
    throw std::runtime_error("Unknown Equation expression kind.");
}


} // namespace

FormulaStructure inspectFormula(const FormulaExpr& expression,std::string_view target) {
    FormulaStructure result;
    const auto visit=[&](const auto& self,const FormulaExpr& node,std::string path)->unsigned {
        const auto index=result.nodes.size();result.nodes.push_back({&node,std::move(path),0});
        if(node.kind()==Kind::Symbol && std::find(result.symbols.begin(),result.symbols.end(),node.name())==result.symbols.end())result.symbols.push_back(node.name());
        if(node.kind()!=Kind::Symbol && node.kind()!=Kind::Constant
            && std::find(result.operations.begin(),result.operations.end(),node.operation())==result.operations.end())result.operations.push_back(node.operation());
        std::vector<unsigned> degree;
        for(std::size_t i=0;i<node.arguments().size();++i)degree.push_back(self(self,node.arguments()[i],result.nodes[index].path+"/"+std::to_string(i)));
        unsigned value=node.kind()==Kind::Symbol && node.name()==target?1:0;
        for(auto d:degree)value=std::max(value,d);
        if(node.kind()==Kind::Multiply)value=std::min(2u,degree[0]+degree[1]);
        if(node.kind()==Kind::Divide && degree[1])value=2;
        if(node.kind()==Kind::Operator && value) {
            if(node.operation()==Math::Op::Pow)value=node.arguments()[1].constant()==0?0:node.arguments()[1].constant()==1?degree[0]:2;
            else if(node.operation()==Math::Op::Diffusion && degree[0])value=2;
            else if(node.operation()==Math::Op::Extension || node.operation()==Math::Op::Source || node.operation()==Math::Op::Sin)value=2;
        }
        result.nodes[index].targetDegree=value;return value;
    };
    result.targetDegree=visit(visit,expression,"root");return result;
}

FormulaExpr explicitFormulaRhs(const Equation& equation,std::string_view target) {
    int count=0;bool hasRest=false;
    const auto remove=[&](const auto& self,const FormulaExpr& node,int sign)->FormulaExpr {
        if(node.kind()==Kind::Operator && node.operation()==Math::Op::Ddt) {
            if(sign!=1 || node.arguments()[0].kind()!=Kind::Symbol || node.arguments()[0].name()!=target)
                throw std::runtime_error("Unsupported: explicit RHS requires positive additive ddt(target); no mass inversion or unknown conversion");
            ++count;return FormulaExpr::zeroLike(node);
        }
        if(node.kind()==Kind::Add || node.kind()==Kind::Subtract) {
            auto left=self(self,node.arguments()[0],sign),right=self(self,node.arguments()[1],node.kind()==Kind::Subtract?-sign:sign);
            return node.withArguments({std::move(left),std::move(right)});
        }
        if(node.kind()==Kind::Negate)return node.withArguments({self(self,node.arguments()[0],-sign)});
        const auto structure=inspectFormula(node,target);
        if(std::find(structure.operations.begin(),structure.operations.end(),Math::Op::Ddt)!=structure.operations.end())
            throw std::runtime_error("Unsupported: time derivative is not an additive unit-mass ddt(target)");
        hasRest=true;return node;
    };
    const auto rest=remove(remove,equation.lhs,1);
    const auto right=inspectFormula(equation.rhs,target);
    if(count!=1 || std::find(right.operations.begin(),right.operations.end(),Math::Op::Ddt)!=right.operations.end())
        throw std::runtime_error("Unsupported: explicit RHS requires exactly one positive ddt(target); algebraic/constraint relations retain their semantics");
    return hasRest?FormulaExpr::subtract(equation.rhs,rest):equation.rhs;
}

FormulaValueKernel compileFormulaValue(const Equation& equation,const FormulaExpr& expression,
        const FormulaOperatorCatalog& providers,std::vector<std::string>& used,
        const std::vector<std::string>& boundSymbols,const StateRegistry* state,std::vector<std::string>* trace) {
    return valueKernel(equation,expression,providers,{},used,boundSymbols.empty()?nullptr:&boundSymbols,state,trace,"rhs");
}

void FormulaOperatorCatalog::add(FormulaOperatorProvider provider) {
    if (provider.id.empty() || !provider.matches
        || (provider.mathematicalOperator==Math::Op::Extension && (provider.extensionName.empty() || !provider.extensionRule)))
        throw std::runtime_error("Equation numerical provider is incomplete.");
    const auto found=std::find_if(providers_.begin(),providers_.end(),
        [&](const auto& value) { return value.id==provider.id; });
    if (found!=providers_.end())
        throw std::runtime_error("Duplicate Equation numerical provider '"
            +provider.id+"'.");
    if(provider.mathematicalOperator==Math::Op::Extension)for(const auto& existing:providers_)
        if(existing.mathematicalOperator==Math::Op::Extension &&
            (existing.extensionName==provider.extensionName || Math::extensionIdentity(existing.extensionName)==Math::extensionIdentity(provider.extensionName)) &&
            (existing.extensionName!=provider.extensionName || existing.extensionRule!=provider.extensionRule))
            throw std::runtime_error("Conflicting numerical extension mathematical rule: "+provider.extensionName);
    providers_.push_back(std::move(provider));
}

const FormulaOperatorProvider& FormulaOperatorCatalog::resolve(
        const Equation& formula,const FormulaExpr& expression,
        const std::vector<FormulaOperatorBinding>& bindings) const {
    for(const auto& binding:bindings)if(binding.formula==formula.id && !binding.occurrence.empty()) {
        const auto left=inspectFormula(formula.lhs),right=inspectFormula(formula.rhs);
        const auto contains=[&](const FormulaStructure& structure){return std::any_of(structure.nodes.begin(),structure.nodes.end(),[&](const auto& node){return node.expression->kind()==FormulaExpr::Kind::Operator && node.expression->occurrence()==binding.occurrence;});};
        if(!contains(left) && !contains(right))throw std::runtime_error("Numerical binding does not address an AST operator: "+formula.id+"/"+binding.occurrence);
    }
    const auto* selected=selectFormulaBinding(
        formula,expression,expression.occurrence(),bindings);
    const FormulaOperatorProvider* resolved=nullptr;
    for (const auto& provider:providers_) {
        if (provider.mathematicalOperator==Math::Op::Extension && (provider.extensionName!=expression.name() || provider.extensionRule!=expression.extensionRule())) continue;
        if (provider.mathematicalOperator!=expression.operation()
            || !provider.matches(expression)
            || (selected && selected->provider!=provider.id))
            continue;
        if (resolved)
            throw std::runtime_error("Equation operator '"+formula.id+"/"
                +expression.occurrence()+"' matches multiple numerical providers.");
        resolved=&provider;
    }
    if (!resolved)
        throw std::runtime_error("Unsupported: Equation operator '"+formula.id+"/"
            +expression.occurrence()+"' has no compatible numerical provider"
            +(selected ? " for binding '"+selected->provider+"'" : "")+".");
    return *resolved;
}


} // namespace SF::System

namespace SF::Discretization {
namespace detail {
inline bool central2Pattern(const System::FormulaExpr& expression) {
    using Kind=System::FormulaExpr::Kind;
    if (expression.kind()!=Kind::Operator || expression.operation()!=Math::Op::Div
        || expression.arguments().size()!=1) return false;
    const auto& product=expression.arguments()[0];
    if (product.kind()!=Kind::Multiply || product.arguments().size()!=2
        || product.arguments()[0].kind()!=Kind::Symbol) return false;
    const auto& gradient=product.arguments()[1];
    return gradient.kind()==Kind::Operator && gradient.operation()==Math::Op::Grad
        && gradient.arguments().size()==1
        && gradient.arguments()[0].kind()==Kind::Symbol;
}

inline std::pair<std::string,std::string> central2Symbols(
        const System::FormulaExpr& expression) {
    if (!central2Pattern(expression))
        throw std::runtime_error("Central2 diffusion requires div(nu * grad(X)).");
    const auto& product=expression.arguments()[0];
    return {product.arguments()[0].name(),product.arguments()[1].arguments()[0].name()};
}


}
System::FormulaOperatorProvider central2ScalarDiffusion(double D) {
    System::FormulaOperatorProvider p;p.id="scalar.central2.diffusion";p.mathematicalOperator=Math::Op::Diffusion;
    p.matches=[](const auto& e){return e.arguments().size()==2;};
    p.compileValue=[D](const auto&){return scalarCentral2Kernel(D);};return p;
}
System::FormulaOperatorProvider conservativeScalarCentral2(double D,int slot,bool primitive) {
    System::FormulaOperatorProvider p;p.id="scalar.conservative.central2";p.mathematicalOperator=Math::Op::Diffusion;
    p.matches=[](const auto&){return true;};p.compileValue=[D,slot,primitive](const auto&){return conservativeCentral2Kernel(D,slot,primitive);};return p;
}
System::FormulaOperatorProvider central2FormulaDiffusion(FormulaCartesianGrid grid) {
    System::FormulaOperatorProvider p;p.id="central2.cartesian.diffusion";p.mathematicalOperator=Math::Op::Div;
    p.matches=detail::central2Pattern;
    p.compileValue=[grid](const auto& e){return central2CartesianValue(grid,detail::central2Symbols(e));};
    p.compileLinear=[grid](const auto& e,const std::string& target){
        const auto symbols=detail::central2Symbols(e);
        if(symbols.second!=target)throw std::runtime_error("Central2 implicit target does not match grad(X)");
        return central2CartesianLinear(grid,symbols);
    };return p;
}
}
