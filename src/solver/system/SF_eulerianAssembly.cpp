#include "SF_eulerianAssembly.h"
#include "SF_eulerianTurbulence.h"
#include <stdexcept>
#include <array>
#include <map>
#include <tuple>
#include <cmath>
#include <algorithm>
#include "SF_runtimeRequirements.h"
#include "core/system/SF_operationIds.h"
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
const PhaseNumericalConfig* requiredPhaseInputs(const Program& program) {
    std::map<std::size_t,std::string> slots;
    std::optional<std::size_t> reference;
    bool active=false;
    for(const auto& call:program.solve.compiledProgram.steps)if(const auto* contract=eulerianAssemblyContract(call)) {
        active=true;
        if(contract->phaseSlot) {
            const auto [it,added]=slots.emplace(*contract->phaseSlot,contract->phaseName);
            if(!added && it->second!=contract->phaseName)throw std::runtime_error("Conflicting phase slot resources: "+call.source.occurrence);
            if(contract->relation==EulerianAssemblyRelation::ReferencePhaseClosure)reference=contract->phaseSlot;
        }
    }
    if(!active)return nullptr;
    if(!program.numerics.phaseOperator || !reference || slots.size()<2)
        throw std::runtime_error("Phase resources require frozen numerical inputs, at least two phases and a reference-phase closure.");
    std::vector<std::string> names;
    for(const auto& [slot,name]:slots) {
        if(slot!=names.size())throw std::runtime_error("Non-contiguous compiled phase storage slots.");
        names.push_back(name);
        const auto prefix="phase"+std::to_string(slot)+".";
        for(const auto& [symbol,key,components]:std::vector<std::tuple<std::string,std::string,int>>{
                {"phaseMass."+name,"mass",1},{"momentum."+name,"momentum",3},{"enthalpy."+name,"enthalpy",1},
                {"alpha."+name,"alpha",1},{"rho."+name,"density",1},{"U."+name,"velocity",3},
                {"h."+name,"primitiveEnthalpy",1},{"T."+name,"temperature",1}}) {
            const auto& state=program.state.at(symbol);
            if(state.storageKey!=prefix+key || state.storageBinding!=StorageBinding::ProviderDistributed
                || state.components!=components || state.componentOffset!=0 || state.readOnly())
                throw std::runtime_error("Unsupported phase backing: "+symbol+" storage="+state.storageKey);
        }
    }
    for(const auto& call:program.solve.compiledProgram.steps)
        if(!eulerianAssemblyContract(call) && !eulerianTurbulenceContract(call))
            throw std::runtime_error("Unsupported phase resource/time-port composition: step="+call.source.occurrence+" method="+call.equationMethod);
    validateEulerianAssemblyBindings(program.solve.compiledProgram,names,*reference);
    const auto view=[&](const std::string& storage,StateViewKind kind,VariableLocation location,OwnershipKind ownership) {
        return std::any_of(program.solve.compiledProgram.stateViews.begin(),program.solve.compiledProgram.stateViews.end(),
            [&](const auto& v){return v.storage==storage && v.kind==kind && v.location==location
                && v.ownership==ownership && v.owner==StateViewOwner::NumericalProvider;});
    };
    if(!view("pressureCorrection",StateViewKind::Correction,VariableLocation::EulerianCell,OwnershipKind::EulerianGlobalDof))
        throw std::runtime_error("Phase resources lack the shared pressure Correction view.");
    for(std::size_t slot=0;slot<names.size();++slot) {
        const auto prefix="phase"+std::to_string(slot)+".";
        if(!view(prefix+"volumeFaceFlux",StateViewKind::Workspace,VariableLocation::EulerianFace,OwnershipKind::CanonicalFace))
            throw std::runtime_error("Phase resources lack canonical paired volume/mass flux workspace: "+prefix);
        for(const auto* suffix:{"volumeFaceFlux","massFaceFlux"})
            if(std::none_of(program.runtime.workspaceRequirements.begin(),program.runtime.workspaceRequirements.end(),
                [&](const auto& w){return w.id==prefix+suffix && w.components==3 && w.location==VariableLocation::EulerianFace && w.ownership==OwnershipKind::CanonicalFace;}))
                throw std::runtime_error("Missing phase canonical workspace requirement: "+prefix+suffix);
    }
    const auto& p=program.state.at("p");
    if(p.storageBinding!=StorageBinding::ProviderDistributed || p.storageKey!="pressure" || p.components!=1 || p.componentOffset!=0 || p.readOnly())
        throw std::runtime_error("Phase resources require the unique shared-pressure array.");
    for(const char* id:{OpIds::EeDtCompute,OpIds::EeStepBegin,OpIds::EeMomentumSolve,
            OpIds::EePressureSolve,OpIds::EePressurePublish,OpIds::EePressureSync,
            OpIds::EeFaceFluxCorrect,OpIds::EeFaceFluxCanonicalPressure,OpIds::EeEnergySolve,
            OpIds::EeStepCommit,OpIds::TimeCommit}) {
        const auto binding=std::find_if(program.runtime.operationBindings.begin(),program.runtime.operationBindings.end(),
            [&](const auto& item){return item.operation==id;});
        if(!requiresOperation(program.runtime,id) || binding==program.runtime.operationBindings.end()
            || (binding->status==BindingStatus::Resolved && binding->provider!="flow.eulerian-pressure"))
            throw std::runtime_error("Missing/conflicting phase operation resource: "+std::string(id));
    }
    if(!std::isfinite(program.numerics.dt.cfl) || program.numerics.dt.cfl<=0
        || !std::isfinite(program.numerics.dt.maxDeltaT) || program.numerics.dt.maxDeltaT<=0)
        throw std::runtime_error("Phase operations require a finite positive compiled dt policy.");
    const auto& inputs=*program.numerics.phaseOperator;
    if(!std::isfinite(inputs.momentumRelaxation) || inputs.momentumRelaxation<=0 || inputs.momentumRelaxation>1
        || !std::isfinite(inputs.pressureRelaxation) || inputs.pressureRelaxation<=0 || inputs.pressureRelaxation>1
        || inputs.reference.referenceCell<0 || !std::isfinite(inputs.reference.referencePressure))
        throw std::runtime_error("Invalid frozen phase relaxation/reference inputs.");
    FDM::validateLinearSolverConfig(inputs.linear.pressure,"phase pressure");
    FDM::validateLinearSolverConfig(inputs.linear.momentum,"phase momentum");
    FDM::validateLinearSolverConfig(inputs.linear.energy,"phase energy");
    FDM::validateLinearSolverConfig(inputs.linear.turbulence,"phase turbulence");
    return &*program.numerics.phaseOperator;
}

void validateEulerianHow(const CompiledExecutionProgram& program) {
    const auto& source=program.root;const auto& calls=program.steps;
    // This backend supports the frozen segregated micro-lifecycle only.
    // Validate authored scopes; their repetition counts remain source HOW values.
    const auto methodOf=[&](const ExecutionScope& node,const char* method) {
        const auto belongs=[&](const ExecutionScope& member) {
            return member.kind==ExecutionKind::EquationCall && std::any_of(calls.begin(),calls.end(),[&](const auto& call) {
                return call.source.occurrence==member.step.occurrence && call.equationMethod==method;
            });
        };
        return node.kind==ExecutionKind::EquationCall ? belongs(node)
            : node.kind==ExecutionKind::Sequence && !node.children.empty()
                && std::all_of(node.children.begin(),node.children.end(),belongs);
    };
    bool compatible=source.kind==ExecutionKind::Sequence && source.children.size()==2
        && source.children[0].kind==ExecutionKind::Loop && source.children[1].kind==ExecutionKind::Commit;
    if (compatible) {
        const auto& outer=source.children[0];
        const bool closure=!outer.children.empty() && methodOf(outer.children[0],"EulerianTurbulenceClosure");
        const std::size_t start=closure ? 1 : 0;
        const bool transport=outer.children.size()>start+4
            && methodOf(outer.children[start+4],"EulerianTurbulenceTransport");
        compatible=outer.children.size()==start+4+(transport ? 1 : 0)
            && methodOf(outer.children[start],"EulerianPhaseContinuity")
            && methodOf(outer.children[start+1],"EulerianPhaseMomentum")
            && outer.children[start+2].kind==ExecutionKind::Loop
            && methodOf(outer.children[start+3],"EulerianPhaseEnthalpy") && (!transport || closure);
        if (compatible) {
            const auto& pressure=outer.children[start+2];
            compatible=pressure.children.size()==1 && pressure.children[0].kind==ExecutionKind::Loop;
            if (compatible) {
                const auto& pass=pressure.children[0];
                compatible=pass.children.size()==3 && methodOf(pass.children[0],"EulerianSharedPressureCorrection")
                    && methodOf(pass.children[1],"EulerianPhaseCorrection") && methodOf(pass.children[2],"EulerianPhaseFluxCorrection");
            }
        }
    }
    if (!compatible) throw std::runtime_error("Unsupported Eulerian HOW topology for the selected segregated backend; source scopes were not reordered.");

    validateEulerianTurbulenceProcess(program);
}
void validateEulerianAssemblyBindings(const CompiledExecutionProgram& program,
        const std::vector<std::string>& names,std::size_t reference) {
    validateEulerianHow(program);
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
