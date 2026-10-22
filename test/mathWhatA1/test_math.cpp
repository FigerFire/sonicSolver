#include "methods/math/SF_signature.h"
#include "methods/math/SF_tensor.h"
#include <cassert>
using namespace SF::Math;
template<class F> void rejects(F f) {bool failed=false;try{f();}catch(const std::runtime_error&){failed=true;}assert(failed);}
int main() {
 static_assert(sizeof(Signature)<=128,"Ordinary signatures must not inline heterogeneous products");
 const auto rho=signature(0,{1,-3,0,0,0,0,0}),u=signature(1,{0,1,-1,0,0,0,0}),p=signature(0,{1,-1,-2,0,0,0,0});
 auto momentum=resultSignature(Op::Multiply,{rho,u});
 assert(compatible(resultSignature(Op::Ddt,{momentum}),resultSignature(Op::Grad,{p})));
 assert(resultSignature(Op::Div,{u}).rank==0);assert(resultSignature(Op::Grad,{u}).rank==2);
 assert(resultSignature(Op::Laplacian,{u}).rank==1);
 assert(equivalent(resultSignature(Op::Laplacian,{u}).dimensions,Dimensions::specified({0,-1,-1,0,0,0,0})));
 assert(resultSignature(Op::Partial,{u}).rank==1);
 assert(equivalent(resultSignature(Op::Partial,{u}).dimensions,Dimensions::specified({0,0,-1,0,0,0,0})));
 assert(resultSignature(Op::Trace,{signature(2,{})}).rank==0);
 assert(resultSignature(Op::Outer,{signature(3,{}),signature(2,{})}).rank==5);
 assert(resultSignature(Op::Outer,{u,u}).rank==2);assert(resultSignature(Op::Dot,{u,u}).rank==0);
 assert(resultSignature(Op::Cross,{u,u}).rank==1);assert(resultSignature(Op::Divide,{momentum,rho}).rank==1);
 assert(equivalent(resultSignature(Op::Pow,{p},.5).dimensions,scaled(p.dimensions,.5)));
 assert(equivalent(Dimensions::specified({0,1+1e-13,0,0,0,0,0}),Dimensions::specified({0,1,0,0,0,0,0})));
 rejects([&]{resultSignature(Op::Add,{u,p});});rejects([&]{resultSignature(Op::Add,{u,resultSignature(Op::Grad,{p})});});
 rejects([&]{resultSignature(Op::Div,{p});});rejects([&]{resultSignature(Op::Multiply,{u,u});});
 rejects([&]{resultSignature(Op::Divide,{p,u});});rejects([&]{resultSignature(Op::Dot,{p,u});});
 rejects([&]{resultSignature(Op::Sin,{u});});rejects([&]{resultSignature(Op::Ddt,{{}});});
 rejects([&]{resolveOperator("unknown");});rejects([&]{resultSignature(Op::Grad,{u,u});});
 auto other=u;other.extent=2;rejects([&]{resultSignature(Op::Dot,{u,other});});
 rejects([&]{resultSignature(static_cast<Op>(10000),{u});});
 assert(resultSignature(Op::Curl,{u}).rank==1);
 assert(resultSignature(Op::DoubleDot,{signature(2,{}),signature(2,{})}).rank==0);
 const Vector3 a{1,2,3},b{4,5,6};assert(dot(a,b)==32);assert((a+b).x==5);
}
