#include "core/system/SF_formula.h"
#include "core/system/SF_stateRegistry.h"
#include <cassert>
#include <type_traits>
using namespace SF;using namespace SF::System;
template<class F> void rejects(F f) {bool failed=false;try{f();}catch(const std::runtime_error&){failed=true;}assert(failed);}
int main() {
 static_assert(!std::is_aggregate_v<FormulaExpr>);
 FormulaSymbols symbols;
 symbols.add("rho",Math::signature(0,{1,-3,0,0,0,0,0}));symbols.add("U",Math::signature(1,{0,1,-1,0,0,0,0}));
 symbols.add("p",Math::signature(0,{1,-1,-2,0,0,0,0}));
 auto rho=symbols.at("rho"),u=symbols.at("U"),p=symbols.at("p");using E=FormulaExpr;
 assert(Math::compatible(E::op(Math::Op::Pow,{rho,E::constantValue(2)}).signature(),E::power(rho,2).signature()));
 rejects([&]{E::op(Math::Op::Pow,{rho,rho});});
 const auto momentum=rho*u;const auto transient=E::op(Math::Op::Ddt,{momentum},"momentum.time");
 const auto convection=E::op(Math::Op::Div,{E::op(Math::Op::Outer,{momentum,u})},"momentum.convection");
 const Equation ns{"momentum",transient+convection,E::negate(E::op(Math::Op::Grad,{p},"momentum.pressure"))};
 assert(ns.lhs.signature().rank==1 && ns.lhs.arguments()[1].occurrence()=="momentum.convection");
 symbols.add("c",Math::signature(0,{}));symbols.add("D",Math::signature(0,{0,2,-1,0,0,0,0}));symbols.add("k",Math::signature(0,{0,0,-1,0,0,0,0}));
 auto c=symbols.at("c");Equation nonlinear{"newPDE",E::op(Math::Op::Ddt,{c}),E::op(Math::Op::Diffusion,{symbols.at("D"),c})-symbols.at("k")*E::power(c,2)};
 Equation zero{"mass",E::op(Math::Op::Ddt,{rho}),E::constantValue(0)};
 assert(Math::compatible(zero.lhs.signature(),zero.rhs.signature()));
 rejects([&]{symbols.at("unknown");});rejects([&]{u+p;});rejects([&]{u+E::op(Math::Op::Grad,{p});});
 rejects([&]{E::op(Math::Op::Div,{p});});rejects([&]{E::op(Math::Op::Dot,{p,u});});
 rejects([&]{E::op(Math::Op::Sin,{u});});rejects([&]{E::symbol("unspecified",{});});
 rejects([&]{Equation bad{"wrong",transient,E::constantValue(1)};});
 rejects([&]{E::op(Math::Op::Diffusion,{E::constantValue(.1),c})+E::op(Math::Op::Ddt,{c});});
 const Math::OperatorRule scalarIdentity{"test.scalarIdentity",[](const std::vector<Math::Signature>& a) {
  if(a.size()!=1 || a[0].rank)throw std::runtime_error("scalar argument required");return a[0];
 }};
 bool diagnostic=false;try{E::op(scalarIdentity,{u});}catch(const std::runtime_error& error){const std::string message=error.what();diagnostic=message.find("test.scalarIdentity")!=std::string::npos && message.find("operand 0")!=std::string::npos && message.find("rank=1")!=std::string::npos && message.find("dimensions=")!=std::string::npos;}
 assert(diagnostic);
 rejects([&]{E::op("unknown",{c});});rejects([&]{c/u;});
 static_assert(std::is_const_v<std::remove_reference_t<decltype(c.arguments())>>);
 auto gradient=E::op(Math::Op::Grad,{p});auto rebuilt=gradient.withArguments({u});
 assert(gradient.signature().rank==1 && rebuilt.signature().rank==2);
 rejects([&]{(u+u).withArguments({u,p});});
 assert(canonicalFormula({"a",rho,rho})!=canonicalFormula({"a",E::symbol("rho",Math::signature(0,{})),E::symbol("rho",Math::signature(0,{}))}));
 StateRegistry state;StateSymbol vector;vector.id="twoComponents";vector.shape=ValueShape::Vector;vector.components=2;vector.signature=Math::signature(1,{});
 rejects([&]{state.add(vector);});vector.signature.extent=2;state.add(vector);
 StateSymbol symmetric;symmetric.id="stress";symmetric.shape=ValueShape::Tensor;symmetric.components=6;symmetric.signature=Math::signature(2,{});symmetric.signature.structure=Math::TensorStructure::Symmetric;state.add(symmetric);
 EquationRegistry registry;registry.add(ns);registry.add(nonlinear);assert(registry.at("newPDE").lhs.operation()==Math::Op::Ddt);
}
