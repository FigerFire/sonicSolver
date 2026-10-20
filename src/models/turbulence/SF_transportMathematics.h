#pragma once

/// @brief Single-fluid RAS mathematical contributions.
/// The stage recipe declares full conservative transport. The explicit legacy
/// selection retains its real no-convection, once-step arithmetic contract.
#include "core/system/SF_formula.h"
#include <stdexcept>

namespace SF::Turbulence {
inline System::Equation stageTransportMathematics(const std::string& target) {
    using E=System::FormulaExpr;auto s=[](const std::string& n){return E::symbol(n);};
    if(target!="k" && target!="omega")throw std::runtime_error("Stage SST target unsupported.");
    auto q=E::multiply(s("rho"),s(target));
    auto rhs=E::add(E::op("diffusion",{E::multiply(s("rho"),s("SST.nu."+target)),s(target)}),
        E::subtract(s("SST.production."+target),E::multiply(s("SST.sink."+target),s(target))));
    if(target=="omega")rhs=E::add(rhs,s("SST.cross"));
    return {target,E::add(E::op("ddt",{q}),E::op("div",{s("rhoU"),s(target)})),rhs,{System::OriginKind::Model,"kOmegaSST.stageConservative"},true};
}
inline System::Equation transportMathematics(const std::string& second,
                                            const std::string& target) {
    using E=System::FormulaExpr;
    const auto s=[](const std::string& name) { return E::symbol(name); };
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
    auto diffusion=E::op("diffusion",{diffusivity,s(target)});
    if (!epsilon && target=="omega")
        diffusion=E::add(std::move(diffusion),coefficient("crossDiffusion"));
    return {target,E::op("ddt",{s(target)}),
        E::add(std::move(source),div(std::move(diffusion),s("rho"))),
        {System::OriginKind::Model,model},true};
}
} // namespace SF::Turbulence
