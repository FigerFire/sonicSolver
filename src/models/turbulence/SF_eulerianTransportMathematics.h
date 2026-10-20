#pragma once
/// @file SF_eulerianTransportMathematics.h
/// @brief 相质量加权湍流 WHAT；对应既有 Eulerian 半隐式输运及模型闭合。
#include "core/system/SF_formula.h"
#include <stdexcept>

namespace SF::Turbulence {
inline System::Equation eulerianClosureMathematics(const std::string& model,
        const std::string& phase) {
    using E=System::FormulaExpr;
    return {"mu_t."+phase,E::symbol("mu_t."+phase),
        E::symbol(model+".eddyDynamicViscosity."+phase),
        {System::OriginKind::Model,model},true};
}
inline System::Equation eulerianTransportMathematics(const std::string& second,
        const std::string& phase,const std::string& variable) {
    using E=System::FormulaExpr;
    if ((second!="epsilon" && second!="omega") || (variable!="k" && variable!=second))
        throw std::runtime_error("Unsupported Eulerian RAS mathematical contract.");
    const bool epsilon=second=="epsilon";
    const std::string model=epsilon ? "kEpsilon" : "kOmegaSST";
    const auto s=[&](const std::string& name) { return E::symbol(name+"."+phase); };
    const auto c=[&](const std::string& name) { return s(model+"."+name); };
    const auto mul=[](E a,E b) { return E::multiply(std::move(a),std::move(b)); };
    const auto div=[](E a,E b) { return E::divide(std::move(a),std::move(b)); };
    const auto production=c("limitedProduction");
    E source,sink,diffusivity;
    if (epsilon) {
        diffusivity=mul(s("alpha"),E::add(s("mu"),div(s("mu_t"),c(variable=="k" ? "sigmaK" : "sigmaEpsilon"))));
        sink=div(mul(s("phaseMass"),s("epsilon")),s("k"));
        source=production;
        if (variable=="epsilon") {
            source=mul(div(mul(c("C1"),s("epsilon")),s("k")),production);
            sink=mul(c("C2"),std::move(sink));
        }
    } else {
        diffusivity=mul(s("alpha"),E::add(s("mu"),mul(c(variable=="k" ? "sigmaK" : "sigmaOmega"),s("mu_t"))));
        sink=mul(c(variable=="k" ? "betaStar" : "beta"),mul(s("phaseMass"),s("omega")));
        source=variable=="k" ? production : E::add(mul(div(mul(c("gamma"),s("rho")),s("mu_t")),production),c("crossSource"));
    }
    const auto transient=E::op("ddt",{mul(s("phaseMass"),s(variable))});
    const auto convection=E::op("div",{mul(s("massFaceFlux"),s(variable))});
    return {variable+"."+phase,E::add(E::add(transient,convection),mul(std::move(sink),s(variable))),
        E::add(E::op("diffusion",{std::move(diffusivity),s(variable)}),std::move(source)),
        {System::OriginKind::Model,model},true};
}
}
