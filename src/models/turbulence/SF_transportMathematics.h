#pragma once

/// @brief Single-fluid RAS mathematical contributions.
/// The stage recipe declares full conservative transport. The explicit legacy
/// selection retains its real no-convection, once-step arithmetic contract.
#include "core/system/SF_formula.h"
#include <stdexcept>

namespace SF::Turbulence {
inline Math::Signature mathematicalSignature(std::string name) {
    const auto scalar=[](std::array<double,7> d){return Math::signature(0,d);};
    const auto prefix=name.find('.');
    if(prefix!=std::string::npos && (name.substr(0,prefix)=="SST" || name.substr(0,prefix)=="kOmegaSST" || name.substr(0,prefix)=="kEpsilon")) {
        name=name.substr(prefix+1);name=name.substr(0,name.find('.',name.rfind("nu.",0)==0 || name.rfind("production.",0)==0 || name.rfind("sink.",0)==0 ? name.find('.')+1:0));
        if(name=="nu.k" || name=="nu.omega")return scalar({0,2,-1,0,0,0,0});
        if(name=="production.k" || name=="limitedProduction")return scalar({1,-1,-3,0,0,0,0});
        if(name=="production.omega" || name=="cross" || name=="crossSource" || name=="crossDiffusion")return scalar({1,-3,-2,0,0,0,0});
        if(name=="sink.k" || name=="sink.omega")return scalar({1,-3,-1,0,0,0,0});
        if(name=="maxRhoK")return scalar({1,-1,-2,0,0,0,0});
        if(name=="maxK")return scalar({0,2,-2,0,0,0,0});
        if(name=="maxMuT" || name=="laminarMu" || name=="eddyDynamicViscosity")return scalar({1,-1,-1,0,0,0,0});
        for(const auto* key:{"betaStar","beta","gamma","C1","C2","sigmaK","sigmaOmega","sigmaEpsilon"})if(name==key)return scalar({});
    } else {
        name=name.substr(0,prefix);
        if(name=="rho" || name=="phaseMass")return scalar({1,-3,0,0,0,0,0});
        if(name=="mu" || name=="mu_t")return scalar({1,-1,-1,0,0,0,0});
        if(name=="k")return scalar({0,2,-2,0,0,0,0});
        if(name=="omega")return scalar({0,0,-1,0,0,0,0});
        if(name=="epsilon")return scalar({0,2,-3,0,0,0,0});
        if(name=="alpha")return scalar({});
        if(name=="massFaceFlux")return scalar({1,0,-1,0,0,0,0});
        if(name=="U")return Math::signature(1,{0,1,-1,0,0,0,0});
        if(name=="rhoU")return Math::signature(1,{1,-2,-1,0,0,0,0});
        if(name=="Smagorinsky")return scalar({1,-1,-1,0,0,0,0});
    }
    throw std::runtime_error("Unknown declared turbulence mathematical symbol: "+name);
}
inline System::FormulaExpr mathematicalSymbol(const std::string& name) {return System::FormulaExpr::symbol(name,mathematicalSignature(name));}
/// Integrated mass rate divided by cell volume is distinct from continuous div.
inline System::FormulaExpr faceMassDivergence(System::FormulaExpr value) {
    static const Math::OperatorRule rule{"faceMassDivergence",[](const std::vector<Math::Signature>& a) {
        if(a.size()!=1 || a[0].rank)throw std::runtime_error("Face mass divergence requires scalar face rate");
        auto result=a[0];result.dimensions=Math::combine(result.dimensions,Math::Dimensions::specified({0,3,0,0,0,0,0}),-1);return result;
    }};return System::FormulaExpr::op(rule,{std::move(value)});
}
}

namespace SF::Turbulence {
inline System::Equation stageTransportMathematics(const std::string& target) {
    using E=System::FormulaExpr;auto s=[](const std::string& n){return mathematicalSymbol(n);};
    if(target!="k" && target!="omega")throw std::runtime_error("Stage SST target unsupported.");
    auto q=E::multiply(s("rho"),s(target));
    auto rhs=E::add(E::op(Math::Op::Diffusion,{E::multiply(s("rho"),s("SST.nu."+target)),s(target)}),
        E::subtract(s("SST.production."+target),E::multiply(s("SST.sink."+target),s(target))));
    if(target=="omega")rhs=E::add(rhs,s("SST.cross"));
    return {target,E::add(E::op(Math::Op::Ddt,{q}),E::op(Math::Op::Div,{s("rhoU"),s(target)})),rhs,{System::OriginKind::Model,"kOmegaSST.stageConservative"}};
}
inline System::Equation transportMathematics(const std::string& second,
                                            const std::string& target) {
    using E=System::FormulaExpr;
    const auto s=[](const std::string& name) { return mathematicalSymbol(name); };
    const auto mul=[](E a,E b) { return E::multiply(std::move(a),std::move(b)); };
    const auto div=[](E a,E b) { return E::divide(std::move(a),std::move(b)); };
    if ((second!="epsilon" && second!="omega") || (target!="k" && target!=second))
        throw std::runtime_error("Unsupported single-fluid RAS mathematical contract.");
    const bool epsilon=second=="epsilon";
    const std::string model=epsilon ? "kEpsilon" : "kOmegaSST";
    const auto coefficient=[&](const std::string& name) { return s(model+"."+name); };
    const auto production=coefficient("limitedProduction");
    E source;
    E diffusivity;
    if (target=="k") {
        const auto sink=epsilon ? s("epsilon")
            : mul(coefficient("betaStar"),mul(s("k"),s("omega")));
        source=E::subtract(div(production,s("rho")),sink);
    } else if (epsilon) {
        source=E::subtract(div(mul(coefficient("C1"),mul(s("epsilon"),production)),
                                  coefficient("maxRhoK")),
            div(mul(coefficient("C2"),mul(s("epsilon"),s("epsilon"))),coefficient("maxK")));
    } else {
        source=E::subtract(div(mul(coefficient("gamma"),production),coefficient("maxMuT")),
            mul(coefficient("beta"),mul(s("omega"),s("omega"))));
    }
    diffusivity=epsilon ? div(s("mu_t"),coefficient(target=="k" ? "sigmaK" : "sigmaEpsilon"))
        : E::add(coefficient("laminarMu"),mul(coefficient(target=="k" ? "sigmaK" : "sigmaOmega"),s("mu_t")));
    auto diffusion=E::op(Math::Op::Diffusion,{diffusivity,s(target)});
    if (!epsilon && target=="omega")
        diffusion=E::add(std::move(diffusion),coefficient("crossDiffusion"));
    return {target,E::op(Math::Op::Ddt,{s(target)}),
        E::add(std::move(source),div(std::move(diffusion),s("rho"))),
        {System::OriginKind::Model,model}};
}
} // namespace SF::Turbulence
