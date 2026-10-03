/// @file SF_termProviderCatalog.cpp
/// @brief Model-independent term capabilities, independent of equation IDs and solve policy.

#include "SF_termProviderCatalog.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {
const FormulaOperatorBinding* selectFormulaBinding(
        const Equation& formula,const FormulaExpr& expression,
        std::string_view occurrence,
        const std::vector<FormulaOperatorBinding>& bindings) {
    const auto priority=[&](const FormulaOperatorBinding& binding) {
        if (binding.formula==formula.id && !occurrence.empty()
            && binding.occurrence==occurrence) return 2;
        if (binding.formula==formula.id && binding.occurrence.empty()) return 1;
        if (binding.formula.empty() && binding.occurrence==expression.name) return 0;
        return -1;
    };
    int highest=-1;
    for (const auto& binding:bindings) highest=std::max(highest,priority(binding));
    const FormulaOperatorBinding* selected=nullptr;
    for (const auto& binding:bindings) {
        if (highest<0 || priority(binding)!=highest) continue;
        if (selected)
            throw std::runtime_error("Equation operator '"+formula.id+"/"
                +std::string(occurrence)+"' has duplicate numerical bindings.");
        selected=&binding;
    }
    return selected;
}


namespace {
// The fused conservative kernels implement these specific mathematical
// operands. Other symbols need their own provider, not the same flux kernel.
bool fusedOperandsSupported(const TermMatchContext& context) {
    if (context.equationMethod!="ConservativeResidual"
        && context.equationMethod!="ConservativePressureMomentum") return true;
    const auto& expression=context.expression;
    const auto symbol=[&](std::size_t n,std::string_view expected) {
        return n<expression.arguments.size()
            && expression.arguments[n].kind==FormulaExpr::Kind::Symbol
            && expression.arguments[n].name==expected;
    };
    if (expression.name=="div") {
        return expression.arguments.size()==1
            && ((context.output=="rho" && symbol(0,"massFlux"))
                || (context.output=="rhoU" && symbol(0,"momentumFlux"))
                || (context.output=="rhoE" && symbol(0,"energyFlux")));
    }
    if (expression.name=="diffusion") {
        return expression.arguments.size()==2
            && ((context.output=="rhoU" && symbol(0,"mu") && symbol(1,"U"))
                || (context.output=="rhoE" && symbol(0,"conductivity") && symbol(1,"T")));
    }
    return true;
}
} // namespace

void TermProviderCatalog::add(TermProviderDescriptor descriptor) {
    if (descriptor.id.empty() || !descriptor.match)
        throw std::invalid_argument("Term provider needs an id and matcher.");
    if (std::any_of(descriptors_.begin(),descriptors_.end(),
            [&](const auto& existing) { return existing.id==descriptor.id; }))
        throw std::invalid_argument("Duplicate term provider: "+descriptor.id);
    descriptors_.push_back(std::move(descriptor));
}

void TermProviderCatalog::addSource(SourceTermProviderDescriptor source) {
    add(TermProviderDescriptor{
        std::move(source.id),std::move(source.match),
        std::move(source.sourceRecipe),
        std::move(source.compilePrimitiveSource),nullptr,
        std::move(source.owner),
        std::move(source.compileConservativeSource)});
}

ResolvedTermProvider TermProviderCatalog::resolve(
        const TermMatchContext& context,
        const std::vector<FormulaOperatorBinding>& bindings) const {
    const auto* selected=selectFormulaBinding(context.formula,context.expression,
                                             context.occurrence,bindings);
    std::vector<const TermProviderDescriptor*> candidates;
    for (const auto& descriptor:descriptors_)
        if ((!selected || descriptor.id==selected->provider)
            && descriptor.match(context)) candidates.push_back(&descriptor);
    std::sort(candidates.begin(),candidates.end(),
              [](const auto* left,const auto* right) {
                  return left->id<right->id;
              });
    ResolvedTermProvider result;
    if (candidates.size()==1) {
        result.id=candidates.front()->id;
        result.status=BindingStatus::Resolved;
        result.owner=candidates.front()->owner;
        if (candidates.front()->compilePrimitiveSource) {
            result.primitiveSource=
                candidates.front()->compilePrimitiveSource();
            if (!result.primitiveSource)
                throw std::runtime_error(
                    "Term provider '"+result.id
                    +"' compiled an empty primitive source kernel.");
            result.compiledDataAvailable=true;
        }
        if (candidates.front()->compileConservativeSource) {
            result.conservativeSource=
                candidates.front()->compileConservativeSource();
            if (!result.conservativeSource)
                throw std::runtime_error(
                    "Term provider '"+result.id
                    +"' compiled an empty conservative source kernel.");
            result.compiledDataAvailable=true;
        }
        result.primitiveSpatial=candidates.front()->primitiveSpatial;
        result.recipe=candidates.front()->sourceRecipe
            ? candidates.front()->sourceRecipe
            : (context.selectedRecipe
                ? std::optional<FDM::TermRecipe>(*context.selectedRecipe)
                : std::nullopt);
    } else if (candidates.empty()) {
        result.reason="no term provider for Equation '"+context.formula.id
            +"' occurrence '"+std::string(context.occurrence)+"'"
            +(selected ? " for binding '"+selected->provider+"'" : "");
    } else {
        result.status=BindingStatus::Invalid;
        result.reason="ambiguous term providers for Equation '"
            +context.formula.id+"' occurrence '"
            +std::string(context.occurrence)+"': ";
        for (std::size_t i=0;i<candidates.size();++i) {
            if (i) result.reason+=", ";
            result.reason+=candidates[i]->id;
        }
    }
    return result;
}

TermProviderCatalog TermProviderCatalog::builtIn() {
    TermProviderCatalog catalog;
    catalog.add({"convection.primitiveUpwind1",
        [](const TermMatchContext& context) {
            return context.expression.kind==FormulaExpr::Kind::Operator
                && context.expression.name=="div"
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::PrimitiveUpwind1
                && context.output=="U";
        },{},nullptr,
        &Discretization::PressureMomentum::primitiveUpwind1});
    catalog.add({"convection.conservativeFlux",
        [](const TermMatchContext& context) {
            return context.expression.kind==FormulaExpr::Kind::Operator
                && context.expression.name=="div"
                && context.selectedRecipe
                && context.selectedRecipe->role()==FDM::TermRole::Convection
                && context.selectedRecipe->id()!=FDM::TermRecipeId::PrimitiveUpwind1
                && fusedOperandsSupported(context);
        },{}});
    catalog.add({"diffusion.central2",
        [](const TermMatchContext& context) {
            return context.expression.kind==FormulaExpr::Kind::Operator
                && context.expression.name=="diffusion"
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::Central2Explicit
                && fusedOperandsSupported(context);
        },{},nullptr,
        &Discretization::PressureMomentum::central2});
    catalog.add({"diffusion.central4",
        [](const TermMatchContext& context) {
            return context.expression.kind==FormulaExpr::Kind::Operator
                && context.expression.name=="diffusion"
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::Central4Explicit
                && fusedOperandsSupported(context);
        },{}});
    return catalog;
}

} // namespace SF::System
