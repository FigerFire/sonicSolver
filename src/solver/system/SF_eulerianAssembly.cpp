#include "SF_eulerianAssembly.h"
#include "SF_eulerianTurbulence.h"
#include <stdexcept>
#include <array>
namespace SF::System {
const char* toString(EulerianAssemblyRelation value) {
    switch (value) {
        case EulerianAssemblyRelation::PhaseMassBalance:return "phase mass balance";
        case EulerianAssemblyRelation::ReferencePhaseClosure:return "reference phase remainder closure";
        case EulerianAssemblyRelation::PhaseMomentumBalance:return "phase momentum predictor";
        case EulerianAssemblyRelation::SharedPressureCorrection:return "shared pressure correction";
        case EulerianAssemblyRelation::PhaseMomentumCorrection:return "phase momentum correction";
        case EulerianAssemblyRelation::SharedPressureUpdate:return "shared pressure update";
        case EulerianAssemblyRelation::PairedFaceFluxCorrection:return "paired volume/mass face-flux correction";
        case EulerianAssemblyRelation::PhaseEnthalpyBalance:return "phase enthalpy with material pressure work";
    }
    throw std::runtime_error("Invalid compiled Eulerian relation.");
}
const char* toString(EulerianSourceExtension value) {
    switch (value) {
        case EulerianSourceExtension::Gravity:return "gravity";
        case EulerianSourceExtension::MRF:return "MRF";
        case EulerianSourceExtension::WallHeat:return "wallHeat";
    }
    throw std::runtime_error("Invalid compiled Eulerian source extension.");
}
void validateEulerianAssemblyBindings(const CompiledExecutionProgram& program,
        const std::vector<std::string>& names,std::size_t reference) {
    if (names.size()<2 || reference>=names.size())
        throw std::runtime_error("Invalid Eulerian runtime phase layout.");
    std::vector<std::array<int,8>> phaseCounts(names.size());
    int pressure=0,publication=0;
    for (const auto& call:program.steps) {
        if (eulerianTurbulenceContract(call)) continue;
        const auto* c=eulerianAssemblyContract(call);
        if (!c || c->provider!="flow.eulerian-pressure" || c->provider!=call.backendProvider
            || c->method!=call.equationMethod || c->equation!=call.source.equation
            || c->target.symbol!=call.target.symbol || c->target.kind!=call.target.kind)
            throw std::runtime_error("Eulerian runtime requires intact frozen provider contracts.");
        using R=EulerianAssemblyRelation;
        if (c->phaseSlot) {
            const auto slot=*c->phaseSlot;
            if (slot>=names.size() || c->phaseName!=names[slot]
                || c->relation==R::SharedPressureCorrection || c->relation==R::SharedPressureUpdate
                || (c->relation==R::ReferencePhaseClosure && slot!=reference)
                || (c->relation==R::PhaseMassBalance && slot==reference))
                throw std::runtime_error("Frozen Eulerian phase slot/reference differs from PhaseSystem backing.");
            ++phaseCounts[slot].at(static_cast<std::size_t>(c->relation));
        } else if (c->relation==R::SharedPressureCorrection) ++pressure;
        else if (c->relation==R::SharedPressureUpdate) ++publication;
        else throw std::runtime_error("Eulerian phase relation lacks a frozen phase slot.");
    }
    if (pressure!=1 || publication!=1)
        throw std::runtime_error("Eulerian compiled contracts lack unique shared pressure solve/update.");
    for (std::size_t slot=0;slot<names.size();++slot) {
        using R=EulerianAssemblyRelation;
        auto expected=std::array<int,8>{};
        for (auto relation:{slot==reference ? R::ReferencePhaseClosure : R::PhaseMassBalance,
                R::PhaseMomentumBalance,R::PhaseMomentumCorrection,R::PairedFaceFluxCorrection,R::PhaseEnthalpyBalance})
            expected[static_cast<std::size_t>(relation)]=1;
        if (phaseCounts[slot]!=expected)
            throw std::runtime_error("Eulerian compiled contracts do not cover each declared phase exactly once.");
    }
}
}
