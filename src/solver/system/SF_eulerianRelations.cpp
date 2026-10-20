#include "SF_eulerianRelations.h"
namespace SF::System {
namespace {
using E=FormulaExpr;
E s(const std::string& name) { return E::symbol(name); }
E op(const std::string& name,std::vector<E> args) { return E::op(name,std::move(args)); }
}
std::vector<Equation> eulerianPhaseRelations(const std::string& phase,bool reference) {
    const auto q=[&](const char* name) { return s(std::string(name)+"."+phase); };
    const std::string suffix="."+phase;
    Equation mass{"E_CONTINUITY"+suffix,
        E::add(op("ddt",{q("phaseMass")}),op("div",{q("phaseMassFlux")})),q("phaseMassSources"),{}};
    if (reference) {
        // The reference continuum realizes sum(alpha)=1, not a second explicit mass update.
        mass.lhs=q("phaseMass");
        mass.rhs=E::multiply(q("rho"),s("referencePhaseRemainder"));
    }
    return {std::move(mass),
        {"momentum"+suffix,E::add(op("ddt",{q("momentum")}),
            E::add(op("div",{q("phaseMomentumFlux")}),E::multiply(q("alpha"),op("gradient",{s("p")})))),
            E::add(op("diffusion",{q("effectivePhaseViscosity"),q("U")}),q("phaseMomentumSources")),{}},
        {"E_ENTHALPY"+suffix,E::add(op("ddt",{q("enthalpy")}),op("div",{q("phaseEnthalpyFlux")})),
            E::add(op("diffusion",{q("effectivePhaseConductivity"),q("h")}),
                E::add(q("phaseEnthalpySources"),E::multiply(q("alpha"),q("materialPressureRate")))),{}},
        {"correctPhase"+suffix,q("momentum"),E::multiply(q("phaseMass"),
            E::subtract(q("predictedVelocity"),E::multiply(E::divide(q("alpha"),q("momentumDiagonal")),op("gradient",{s("pPrime")})))),{}},
        {"correctPhaseFlux"+suffix,op("pairedFaceFlux",{q("volumeFaceFlux"),q("massFaceFlux")}),
            op("pairedFaceFlux",{
                E::subtract(q("volumeFaceFlux"),E::multiply(q("facePressureResponse"),op("faceGradient",{s("pPrime")}))),
                E::subtract(q("massFaceFlux"),E::multiply(q("faceDensity"),
                    E::multiply(q("facePressureResponse"),op("faceGradient",{s("pPrime")}))))}),{}}

    };
}
std::vector<Equation> eulerianPressureRelations() {
    return {
        {"E_SHARED_PRESSURE",E::subtract(E::multiply(s("mixtureCompressibilityOverDt"),s("pPrime")),
            op("diffusion",{s("sumFaceAlphaSquaredOverDiagonal"),s("pPrime")})),
            E::subtract(E::subtract(s("sumMassSourceOverDensity"),op("div",{s("mixtureVolumeFlux")})),s("pressureHistory")),{}},
        {"correctSharedP",s("p"),E::add(s("p"),E::multiply(s("pressureRelaxation"),s("pPrime"))),{}}
    };
}
}
