/// @file SF_termProviderCatalog.cpp
/// @brief Model-independent term capabilities, independent of equation IDs and solve policy.

#include "SF_termProviderCatalog.h"

#include <algorithm>
#include <stdexcept>
#include "solver/discretization/source/SF_scalarFeedback.h"

namespace SF::System {
const FormulaOperatorBinding* selectFormulaBinding(
        const Equation& formula,const FormulaExpr& expression,
        std::string_view occurrence,
        const std::vector<FormulaOperatorBinding>& bindings) {
    const auto priority=[&](const FormulaOperatorBinding& binding) {
        if (binding.formula==formula.id && !occurrence.empty()
            && binding.occurrence==occurrence) return 2;
        if (binding.formula==formula.id && binding.occurrence.empty()) return 1;
        if (binding.formula.empty() && binding.operation && *binding.operation==expression.operation()) return 0;
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
        return n<expression.arguments().size()
            && expression.arguments()[n].kind()==FormulaExpr::Kind::Symbol
            && expression.arguments()[n].name()==expected;
    };
    if (expression.operation()==Math::Op::Div) {
        return expression.arguments().size()==1
            && ((context.output=="rho" && symbol(0,"massFlux"))
                || (context.output=="rhoU" && symbol(0,"momentumFlux"))
                || (context.output=="rhoE" && symbol(0,"energyFlux")));
    }
    if (expression.operation()==Math::Op::Diffusion) {
        return expression.arguments().size()==2
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
    TermProviderDescriptor descriptor{
        std::move(source.id),std::move(source.match),
        std::move(source.sourceRecipe),
        std::move(source.compilePrimitiveSource),nullptr,
        std::move(source.owner),
        std::move(source.compileConservativeSource)};
    descriptor.operations={Math::Op::Source};
    add(std::move(descriptor));
}

ResolvedTermProvider TermProviderCatalog::resolve(
        const TermMatchContext& context,
        const std::vector<FormulaOperatorBinding>& bindings) const {
    const auto* selected=selectFormulaBinding(context.formula,context.expression,
                                             context.occurrence,bindings);
    std::vector<const TermProviderDescriptor*> candidates;
    for (const auto& descriptor:descriptors_)
        if ((!selected || descriptor.id==selected->provider)
            && (descriptor.operations.empty() || std::find(descriptor.operations.begin(),descriptor.operations.end(),context.expression.operation())!=descriptor.operations.end())
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
        if (candidates.front()->compileStageSource) {
            result.stageSource=candidates.front()->compileStageSource(context);
            if (!result.stageSource.evaluate) throw std::runtime_error("Empty compiled Stage source.");
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
    const auto add=[&](Math::Op op,TermProviderDescriptor value){value.operations={op};catalog.add(std::move(value));};
    for(bool work:{false,true}) {
        TermProviderDescriptor p;p.id=work?"source.stage.scalar-work":"source.stage.scalar-force";
        p.match=[work](const TermMatchContext& c) {
            const auto& e=c.expression;
            return c.equationMethod=="ConservativeResidual" && c.output==(work?"rhoE":"rhoU")
                && e.operation()==Math::Op::Source && e.arguments().size()==4
                && e.arguments()[0].kind()==FormulaExpr::Kind::Symbol
                && e.arguments()[0].name()==(work?"linearScalarWork":"linearScalarForce");
        };
        p.compileStageSource=[work](const auto& c){
            const auto& a=c.expression.arguments();using K=FormulaExpr::Kind;
            if(a.size()!=4 || a[1].kind()!=K::Symbol || a[2].kind()!=K::Constant || !std::isfinite(a[2].constant())
                || a[3].kind()!=K::Constant || a[3].constant()<0 || a[3].constant()>2 || a[3].constant()!=std::floor(a[3].constant()))
                throw std::runtime_error("Scalar feedback requires Stage symbol, finite coefficient and axis 0/1/2.");
            return Discretization::linearScalarFeedback(a[1].name(),a[2].constant(),static_cast<int>(a[3].constant()),work);
        };
        p.operations={Math::Op::Source};catalog.add(std::move(p));
    }
    add(Math::Op::Div,{ "convection.primitiveUpwind1",
        [](const TermMatchContext& context) {
            return context.expression.kind()==FormulaExpr::Kind::Operator
                && context.expression.operation()==Math::Op::Div
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::PrimitiveUpwind1
                && context.output=="U";
        },{},nullptr,
        &Discretization::PressureMomentum::primitiveUpwind1});
    add(Math::Op::Div,{ "convection.conservativeFlux",
        [](const TermMatchContext& context) {
            return context.expression.kind()==FormulaExpr::Kind::Operator
                && context.expression.operation()==Math::Op::Div
                && context.selectedRecipe
                && context.selectedRecipe->role()==FDM::TermRole::Convection
                && context.selectedRecipe->id()!=FDM::TermRecipeId::PrimitiveUpwind1
                && fusedOperandsSupported(context);
        },{}});
    add(Math::Op::Diffusion,{ "diffusion.central2",
        [](const TermMatchContext& context) {
            return context.expression.kind()==FormulaExpr::Kind::Operator
                && context.expression.operation()==Math::Op::Diffusion
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::Central2Explicit
                && fusedOperandsSupported(context);
        },{},nullptr,
        &Discretization::PressureMomentum::central2});
    add(Math::Op::Diffusion,{ "diffusion.central4",
        [](const TermMatchContext& context) {
            return context.expression.kind()==FormulaExpr::Kind::Operator
                && context.expression.operation()==Math::Op::Diffusion
                && context.selectedRecipe
                && context.selectedRecipe->id()==FDM::TermRecipeId::Central4Explicit
                && fusedOperandsSupported(context);
        },{}});
    return catalog;
}

} // namespace SF::System
