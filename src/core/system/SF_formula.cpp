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

void validateExpr(const FormulaExpr& expression) {
    using Kind = FormulaExpr::Kind;
    const auto size = expression.arguments().size();
    switch (expression.kind()) {
        case Kind::Symbol:
            if (expression.name().empty() || size) break;
            return;
        case Kind::Constant:
            if (std::isfinite(expression.constant()) && !size) return;
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
            if (!expression.name().empty() && size) break;
            throw std::runtime_error("Equation mathematical operator requires a name and operands.");
    }
    if (expression.kind() == Kind::Symbol || expression.kind() == Kind::Constant)
        throw std::runtime_error("Equation symbol or constant has invalid data.");
    for (const auto& argument : expression.arguments()) validateExpr(argument);
}

void validateOccurrences(const FormulaExpr& expression,
                         std::unordered_set<std::string>& paths) {
    if (expression.kind()==FormulaExpr::Kind::Operator
        && !expression.occurrence().empty()
        && !paths.insert(expression.occurrence()).second)
        throw std::runtime_error("Equation has duplicate operator occurrence '"
                                 +expression.occurrence()+"'.");
    for (const auto& child:expression.arguments())
        validateOccurrences(child,paths);
}

void validateExtensionRules(const std::vector<const Equation*>& equations) {
    std::vector<const FormulaExpr*> rules;
    const auto visit=[&](const auto& self,const FormulaExpr& node)->void {
        if(node.kind()==FormulaExpr::Kind::Operator && node.operation()==Math::Op::Extension) {
            for(const auto* rule:rules)if(rule->name()==node.name() || rule->extensionIdentity()==node.extensionIdentity()) {
                if(rule->name()!=node.name() || rule->extensionRule()!=node.extensionRule())
                    throw std::runtime_error("Conflicting mathematical extension identity/rule: "+node.name());
            }
            rules.push_back(&node);
        }
        for(const auto& child:node.arguments())self(self,child);
    };
    for(const auto* equation:equations){visit(visit,equation->lhs);visit(visit,equation->rhs);}
}

void validateFormula(const Equation& formula) {
    formula.validateTypes();
    validateExpr(formula.lhs);
    validateExpr(formula.rhs);
    std::unordered_set<std::string> paths;
    validateOccurrences(formula.lhs,paths);
    validateOccurrences(formula.rhs,paths);
}

void encode(std::ostringstream& output, const FormulaExpr& expression) {
    output << static_cast<int>(expression.kind()) << '(';
    output << expression.name().size() << ':' << expression.name();
    output << '{' << Math::text(expression.signature()) << '}';
    if (expression.kind() == FormulaExpr::Kind::Constant)
        output << std::setprecision(17) << expression.constant();
    output << '[';
    for (const auto& argument : expression.arguments()) encode(output,argument);
    output << "])";
}

std::string textOf(const FormulaExpr& expression) {
    using Kind=FormulaExpr::Kind;
    if (expression.kind()==Kind::Symbol) return expression.name();
    if (expression.kind()==Kind::Constant) {
        std::ostringstream output;
        output << std::setprecision(17) << expression.constant();
        return output.str();
    }
    if (expression.kind()==Kind::Negate)
        return "-("+textOf(expression.arguments()[0])+")";
    if (expression.kind()==Kind::Operator) {
        std::string result=expression.name()+"(";
        for (std::size_t i=0;i<expression.arguments().size();++i) {
            if (i) result+=", ";
            result+=textOf(expression.arguments()[i]);
        }
        return result+")";
    }
    const char* operation=expression.kind()==Kind::Add?" + "
        : expression.kind()==Kind::Subtract?" - "
        : expression.kind()==Kind::Multiply?" * ":" / ";
    return "("+textOf(expression.arguments()[0])+operation
        +textOf(expression.arguments()[1])+")";
}

} // namespace

FormulaExpr FormulaExpr::symbol(std::string name,Math::Signature signature) {
    if(name.empty())throw std::runtime_error("Empty mathematical symbol");
    Math::requireKnown(signature);FormulaExpr value;
    value.kind_=Kind::Symbol;value.name_=std::move(name);value.signature_=signature;value.literalZero_=false;return value;
}
FormulaExpr FormulaExpr::constantValue(double value) {
    auto result=constantValue(value,Math::signature(0,{}));result.literalZero_=value==0;return result;
}
FormulaExpr FormulaExpr::constantValue(double value,Math::Signature signature) {
    Math::requireKnown(signature);
    if(!std::isfinite(value) || signature.rank)throw std::runtime_error("Mathematical constant requires finite scalar");
    FormulaExpr result;result.constant_=value;result.signature_=signature;result.literalZero_=false;return result;
}
FormulaExpr FormulaExpr::zeroLike(const FormulaExpr& other) {
    FormulaExpr result;result.signature_=other.signature_;result.literalZero_=false;return result;
}
FormulaExpr FormulaExpr::make(Math::Op operation,std::vector<FormulaExpr> args,std::string occurrence) {
    if(operation==Math::Op::Pow) {
        if(args.size()!=2 || args[1].kind()!=Kind::Constant)throw std::runtime_error("Power needs a constant exponent");
        Math::requireCompatible(args[1].signature(),Math::signature(0,{}),"power exponent");
        return power(std::move(args[0]),args[1].constant()).atOccurrence(std::move(occurrence));
    }
    if((operation==Math::Op::Add || operation==Math::Op::Subtract) && args.size()==2) {
        if(args[0].literalZero_)args[0]=zeroLike(args[1]);
        if(args[1].literalZero_)args[1]=zeroLike(args[0]);
    }
    std::vector<Math::Signature> types;for(const auto& a:args)types.push_back(a.signature_);
    FormulaExpr result;result.signature_=Math::resultSignature(operation,types);
    result.operation_=operation;result.literalZero_=false;result.kind_=Kind::Operator;
    switch(operation) {
    case Math::Op::Add:result.kind_=Kind::Add;break;case Math::Op::Subtract:result.kind_=Kind::Subtract;break;
    case Math::Op::Multiply:result.kind_=Kind::Multiply;break;case Math::Op::Divide:result.kind_=Kind::Divide;break;
    case Math::Op::Negate:result.kind_=Kind::Negate;break;default:result.name_=Math::name(operation);break;
    }
    result.arguments_=std::move(args);result.occurrence_=std::move(occurrence);return result;
}
FormulaExpr FormulaExpr::add(FormulaExpr a,FormulaExpr b) {return make(Math::Op::Add,{std::move(a),std::move(b)});}
FormulaExpr FormulaExpr::subtract(FormulaExpr a,FormulaExpr b) {return make(Math::Op::Subtract,{std::move(a),std::move(b)});}
FormulaExpr FormulaExpr::multiply(FormulaExpr a,FormulaExpr b) {return make(Math::Op::Multiply,{std::move(a),std::move(b)});}
FormulaExpr FormulaExpr::divide(FormulaExpr a,FormulaExpr b) {return make(Math::Op::Divide,{std::move(a),std::move(b)});}
FormulaExpr FormulaExpr::negate(FormulaExpr a) {return make(Math::Op::Negate,{std::move(a)});}
FormulaExpr FormulaExpr::op(Math::Op id,std::vector<FormulaExpr> args,std::string occurrence) {return make(id,std::move(args),std::move(occurrence));}
FormulaExpr FormulaExpr::op(std::string name,std::vector<FormulaExpr> args,std::string occurrence) {
    return op(Math::resolveOperator(name),std::move(args),std::move(occurrence));
}
FormulaExpr FormulaExpr::op(const Math::OperatorRule& rule,std::vector<FormulaExpr> args,std::string occurrence) {
    if(!rule.id || !*rule.id || !rule.result)throw std::runtime_error("Invalid mathematical extension rule");
    std::vector<Math::Signature> types;for(const auto& a:args)types.push_back(a.signature_);
    FormulaExpr result;
    try {result.signature_=rule.result(types);Math::requireKnown(result.signature_);}
    catch(const std::runtime_error& error) {
        std::ostringstream message;message<<"Math operator "<<rule.id<<": ";
        for(std::size_t i=0;i<types.size();++i)message<<"operand "<<i<<" {"<<Math::text(types[i])<<"}; ";
        throw std::runtime_error(message.str()+"reason: "+error.what());
    }
    result.literalZero_=false;result.kind_=Kind::Operator;result.operation_=Math::Op::Extension;result.extensionResult_=rule.result;result.extensionIdentity_=Math::extensionIdentity(rule.id);
    result.name_=rule.id;result.arguments_=std::move(args);result.occurrence_=std::move(occurrence);return result;
}
FormulaExpr FormulaExpr::power(FormulaExpr a,double exponent) {
    FormulaExpr result;result.signature_=Math::resultSignature(Math::Op::Pow,{a.signature_,Math::signature(0,{})},exponent);
    result.kind_=Kind::Operator;result.operation_=Math::Op::Pow;result.name_="pow";result.literalZero_=false;
    result.arguments_={std::move(a),constantValue(exponent)};return result;
}
FormulaExpr FormulaExpr::withArguments(std::vector<FormulaExpr> values) const {
    if(kind_==Kind::Symbol || kind_==Kind::Constant)throw std::runtime_error("Leaf AST cannot acquire operands");
    if(extensionResult_)return op(Math::OperatorRule{name_.c_str(),extensionResult_},std::move(values),occurrence_);
    if(operation_==Math::Op::Pow) {
        if(values.size()!=2 || values[1].kind()!=Kind::Constant)throw std::runtime_error("Power needs a constant exponent");
        Math::requireCompatible(values[1].signature(),Math::signature(0,{}),"power exponent");
        return power(std::move(values[0]),values[1].constant()).atOccurrence(occurrence_);
    }
    return make(operation_,std::move(values),occurrence_);
}

Equation::Equation(std::string identity,FormulaExpr left,FormulaExpr right,Provenance source)
    :id(std::move(identity)),lhs(std::move(left)),rhs(std::move(right)),origin(std::move(source)) {
    contextualizeZero();validateTypes();
}
void Equation::contextualizeZero() {
    if(lhs.literalZero_)lhs=FormulaExpr::zeroLike(rhs);
    if(rhs.literalZero_)rhs=FormulaExpr::zeroLike(lhs);
}
void Equation::validateTypes() const {Math::requireCompatible(lhs.signature(),rhs.signature(),"equation "+id);}

void EquationRegistry::add(Equation formula) {
    formula.contextualizeZero();
    if (formula.id.empty()) throw std::runtime_error("Equation ID is empty.");
    validateFormula(formula);
    std::vector<const Equation*> checked{&formula};
    for(const auto& current:formulas_)if(current.id!=formula.id)checked.push_back(&current);
    validateExtensionRules(checked);
    const auto found=std::find_if(formulas_.begin(),formulas_.end(),
        [&](const Equation& current) { return current.id==formula.id; });
    if (found!=formulas_.end())
        throw std::runtime_error("Duplicate Equation ID '"+formula.id+"'.");
    formulas_.push_back(std::move(formula));
}

void EquationRegistry::replace(Equation formula) {
    formula.contextualizeZero();
    validateFormula(formula);
    std::vector<const Equation*> checked{&formula};
    for(const auto& current:formulas_)if(current.id!=formula.id)checked.push_back(&current);
    validateExtensionRules(checked);
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
    if (expression.kind()==FormulaExpr::Kind::Symbol
        && expression.name()==symbol) return true;
    return std::any_of(expression.arguments().begin(),expression.arguments().end(),
        [&](const FormulaExpr& argument) { return dependsOn(argument,symbol); });
}

} // namespace SF::System
