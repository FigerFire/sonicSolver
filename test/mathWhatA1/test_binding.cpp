#include "solver/system/SF_formulaCompiler.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "solver/system/SF_termProviderCatalog.h"
#include "solver/system/SF_builtinState.h"
#include <cassert>
using namespace SF;using namespace SF::System;
int main() {
 BuiltinStateCatalog builtin;assert(builtin.expression("massFlux").signature().rank==1);
 bool noStorage=false;try{builtin.solution("massFlux");}catch(const std::runtime_error&){noStorage=true;}assert(noStorage);
 using E=FormulaExpr;FormulaSymbols symbols;
 symbols.add("mu",Math::signature(0,{1,-1,-1,0,0,0,0}));symbols.add("U",Math::signature(1,{0,1,-1,0,0,0,0}));
 auto diffusion=E::op(Math::Op::Diffusion,{symbols.at("mu"),symbols.at("U")},"momentum.viscous");
 const Equation equation{"momentum",diffusion,E::constantValue(0),{OriginKind::User,"A1 binding"}};
 const auto canonical=canonicalFormula(equation);auto catalog=TermProviderCatalog::builtIn();
 for(auto id:{FDM::TermRecipeId::Central2Explicit,FDM::TermRecipeId::Central4Explicit}) {
  auto recipe=FDM::builtInDiffusionRecipe(id==FDM::TermRecipeId::Central2Explicit?FDM::ViscousScheme::Central2:FDM::ViscousScheme::Central4);
  auto chosen=catalog.resolve({equation,diffusion,"momentum.viscous","rhoU","ConservativeResidual",&recipe});
  assert(chosen.status==BindingStatus::Resolved);assert(chosen.id==(id==FDM::TermRecipeId::Central2Explicit?"diffusion.central2":"diffusion.central4"));
  assert(canonicalFormula(equation)==canonical);
 }
 auto recipe=FDM::builtInDiffusionRecipe(FDM::ViscousScheme::Central2);
 auto incompatible=catalog.resolve({equation,diffusion,"momentum.viscous","rhoU","ConservativeResidual",&recipe},{{"momentum","momentum.viscous","convection.conservativeFlux"}});
 assert(incompatible.status==BindingStatus::Unsupported && !incompatible.reason.empty());
 auto c=E::symbol("c",Math::signature(0,{}));
 auto rhs=E::op(Math::Op::Diffusion,{E::constantValue(.1,Math::signature(0,{0,2,-1,0,0,0,0})),c},"spread");
 const Equation scalar{"scalar",E::op(Math::Op::Ddt,{c}),rhs};FormulaOperatorCatalog values;
 StateRegistry state;const std::vector<std::string> slots={"c"};
 auto first=Discretization::central2ScalarDiffusion(state,slots),second=first;
 second.id="test.central2.independent-registration";values.add(first);values.add(second);
 FormulaValues data;data.spacing=[](int d){return d==0?1.:0.;};data.boundNeighbor=[](std::size_t,int,int,int sign,int){return sign<0?1.:sign>0?9.:4.;};
 for(const auto& id:{first.id,second.id}) {
  const auto& selected=values.resolve(scalar,rhs,{{"scalar","spread",id}});
  assert(selected.mathematicalOperator==rhs.operation());
  auto kernel=selected.compileValue(rhs);assert(kernel(data,0,0)==.2);
 }
 bool rejected=false;try{values.resolve(scalar,E::op(Math::Op::Grad,{c},"spread"),{{"scalar","spread",first.id}});}catch(const std::runtime_error&){rejected=true;}assert(rejected);
}
