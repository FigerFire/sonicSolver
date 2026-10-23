/// @brief 既有融合后端的操作数限制、term 能力和 lowering；保留原数学及执行顺序。
#include "SF_fusedSpatial.h"
#include "solver/discretization/source/SF_scalarFeedback.h"
#include <algorithm>
#include <stdexcept>
namespace SF::System::FusedSpatial {
namespace {
void appendUnique(std::vector<std::string>& values,std::string value) {
    if (std::find(values.begin(),values.end(),value)==values.end()) values.push_back(std::move(value));
}
struct SelectedFormula {
    std::string formula;
    std::string execution;
    std::string output;
    std::string method;
    bool primitiveSourceRequired = false;
    std::string residualWorkspace;
};

std::vector<SelectedFormula> selectFormulas(
        const CompiledExecutionProgram& program) {
    std::vector<SelectedFormula> selected;
    for (const auto& step:program.steps) {
        if (step.spatialTerms) {
            for (const auto& call:step.calls)
                selected.push_back({call.equation,step.source.equation,
                                    call.target,step.equationMethod,step.primitiveSourceRequired,step.residualWorkspace});
        } else if (!step.sourceMathInputs.empty()) {
            for (const auto& formula:step.sourceMathInputs)
                selected.push_back({formula,step.source.equation,
                    step.source.target.symbol,step.equationMethod,step.primitiveSourceRequired,step.residualWorkspace});
        }
        // Other methods own their local numerical realization. The spatial
        // compiler consumes the method's declared residual/source input set.
    }
    return selected;
}

std::string symbolAt(const FormulaExpr& expression,std::size_t index,
                     const Equation& formula,const std::string& path) {
    if (expression.arguments().size()<=index
        || expression.arguments()[index].kind()!=FormulaExpr::Kind::Symbol)
        throw std::runtime_error("Equation '"+formula.id+"' occurrence '"+path
            +"' has an unsupported non-symbol operator operand.");
    return expression.arguments()[index].name();
}

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


void registerProviders(TermProviderCatalog& catalog) {
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
}
void lower(const ExecutableEquationSystem& equations,const CompiledEquationCall& call,
        const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) {
    if (recipes.time.topology()!=FDM::TimeTopology::ExplicitStages)
        throw std::runtime_error("Current fused numerical recipes require an explicit time topology.");
    CompiledExecutionProgram selection;selection.steps.push_back(call);
    const auto selectedFormulas=selectFormulas(selection);
    const auto bindRecipe=[&](const FDM::TermRecipe& recipe,RecipeConsumerKind kind,std::string consumer) {
        result.recipeBindings.push_back({recipe,kind,std::move(consumer)});
    };
    for (const auto& selected:selectedFormulas) {
        const auto& formula=equations.registry.at(selected.formula);
        int transientCount=0;
        const auto visit=[&](const auto& self,const FormulaExpr& expression,
                             CompiledSpatialBinding::Side side,
                             const std::string& path)->void {
            if (expression.kind()==FormulaExpr::Kind::Subtract
                || expression.kind()==FormulaExpr::Kind::Multiply
                || expression.kind()==FormulaExpr::Kind::Divide
                || expression.kind()==FormulaExpr::Kind::Negate)
                throw std::runtime_error("Equation '"+formula.id
                    +"' at '"+path+"' requires unsupported signed/coefficient "
                    "lowering in the fused numerical backend.");
            if (side==CompiledSpatialBinding::Side::Right && expression.kind()==FormulaExpr::Kind::Symbol)
                throw std::runtime_error("Equation '"+formula.id+"' at '"+path
                    +"' requires an explicit Source AST for RHS symbol '"+expression.name()+"'.");
            const bool source=expression.kind()==FormulaExpr::Kind::Operator
                && expression.operation()==Math::Op::Source;
            const bool convection=expression.kind()==FormulaExpr::Kind::Operator
                && expression.operation()==Math::Op::Div;
            const bool diffusion=expression.kind()==FormulaExpr::Kind::Operator
                && expression.operation()==Math::Op::Diffusion;
            if (source
                    && side!=CompiledSpatialBinding::Side::Right
                || (convection || diffusion)
                    && side!=CompiledSpatialBinding::Side::Left)
                throw std::runtime_error("Equation '"+formula.id+"' at '"+path
                    +"' requires unsupported term-side/sign lowering.");
            if (source || convection || diffusion) {
                const FormulaExpr* numerical=&expression;
                const std::string occurrence=expression.occurrence().empty()
                    ? path : expression.occurrence();
                const auto opName=numerical->operation();
                const std::string primary=symbolAt(*numerical,
                    diffusion?1:0,formula,occurrence);
                const std::string secondary=diffusion
                    ? symbolAt(*numerical,0,formula,occurrence) : "";
                const FDM::TermRecipe* requested=nullptr;
                if (convection) {
                    if (!recipes.convection)
                        throw std::runtime_error("Equation '"+formula.id
                            +"' occurrence '"+occurrence
                            +"' has no convection TermRecipe.");
                    requested=&*recipes.convection;
                } else if (diffusion) {
                    if (!recipes.diffusion)
                        throw std::runtime_error("Equation '"+formula.id
                            +"' occurrence '"+occurrence
                            +"' has no diffusion TermRecipe.");
                    requested=&*recipes.diffusion;
                }
                const auto provider=catalog.resolve({formula,*numerical,
                    occurrence,selected.output,selected.method,requested},bindings);
                const bool isSource=source;
                const bool primitiveSource=isSource && selected.primitiveSourceRequired
                    && provider.primitiveSource;
                if (provider.status!=BindingStatus::Resolved
                    || (!provider.recipe && !primitiveSource && !provider.stageSource.evaluate))
                    throw std::runtime_error("Equation provider "
                        +std::string(provider.status==BindingStatus::Invalid
                            ? "Invalid: " : "Unsupported: ")
                        +(provider.reason.empty()
                            ? "missing recipe for '"+occurrence+"'"
                            : provider.reason));
                if (isSource && !selected.primitiveSourceRequired
                    && !provider.conservativeSource && !provider.stageSource.evaluate)
                    throw std::runtime_error("Conservative source provider '"
                        +provider.id+"' has no compiled execution kernel.");
                result.operators.emplace_back(formula.id,selected.execution,
                    selected.output,selected.method,occurrence,opName,
                    primary,secondary,side,provider.recipe,provider.id,
                    provider.owner,provider.compiledDataAvailable,
                    provider.primitiveSource,provider.conservativeSource,
                    provider.primitiveSpatial);
                result.operators.back().stageSource=provider.stageSource;
                appendUnique(result.providerRequirements,provider.id);
                if (requested) {
                    bindRecipe(*requested,RecipeConsumerKind::EquationTerm,
                        formula.id+"@"+occurrence);
                    if (convection) {
                        result.requiredHaloWidth=std::max(
                            result.requiredHaloWidth,requested->haloWidth());
                        appendUnique(result.workspaceRequirements,selected.residualWorkspace);
                        if (requested->reconstruction()
                            ==FDM::ReconstructionVariable::Characteristic)
                            appendUnique(result.workspaceRequirements,
                                         "characteristicReconstruction");
                    } else {
                        appendUnique(result.workspaceRequirements,
                            selected.residualWorkspace);
                    }
                } else {
                    if (provider.recipe)
                        bindRecipe(*provider.recipe,RecipeConsumerKind::EquationTerm,
                            formula.id+"@"+occurrence);
                    appendUnique(result.workspaceRequirements,
                        selected.residualWorkspace);
                }
                return;
            }
            if (expression.kind()==FormulaExpr::Kind::Operator
                && expression.operation()!=Math::Op::Ddt)
                throw std::runtime_error("Equation '"+formula.id
                    +"' has unsupported numerical operator '"+expression.name()
                    +"' at '"+path+"'.");
            if (expression.kind()==FormulaExpr::Kind::Operator
                && expression.operation()==Math::Op::Ddt) {
                if (++transientCount!=1 || expression.arguments().size()!=1
                    || side!=CompiledSpatialBinding::Side::Left
                    || symbolAt(expression,0,formula,path)!=selected.output)
                    throw std::runtime_error("Equation ddt does not match HOW Output.");
                return;
            }
            if (expression.kind()==FormulaExpr::Kind::Symbol
                || (expression.kind()==FormulaExpr::Kind::Constant
                    && expression.constant()!=0.0))
                throw std::runtime_error("Equation '"+formula.id+"' at '"+path
                    +"' has unsupported unbound residual data.");
            for (std::size_t i=0;i<expression.arguments().size();++i)
                self(self,expression.arguments()[i],side,
                    path+"/"+std::to_string(i));
        };
        visit(visit,formula.lhs,CompiledSpatialBinding::Side::Left,"lhs");
        visit(visit,formula.rhs,CompiledSpatialBinding::Side::Right,"rhs");
        if (transientCount!=1)
            throw std::runtime_error("Fused transport Equation requires exactly one ddt(Output).");
    }
}
} // namespace SF::System::FusedSpatial
