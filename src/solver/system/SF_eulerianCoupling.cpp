#include "SF_eulerianCoupling.h"
#include <stdexcept>
#include <algorithm>
namespace SF::System {
void applyEulerianExecution(ExecutionProgram& program,std::vector<NumericalBinding>& bindings,
        const CouplingPresetRequest& request) {
    if (!request.explicitlyRegistered || request.presetKind!=FDM::PressureCouplingPreset::PIMPLE)
        throw std::runtime_error("Unsupported Eulerian shared-pressure coupling: explicitly selected PIMPLE is required; HOW was not replaced.");
    if (request.outerCorrectors<1 || request.pressureCorrectors<1 || request.nonOrthogonalCorrectors<0)
        throw std::runtime_error("Invalid Eulerian HOW loop counts.");
    const Provenance origin=request.origin;
    const auto scope=[&](ExecutionKind kind,const std::string& id,int count=1) {
        ExecutionScope node;node.kind=kind;node.id=id;node.origin=origin;node.repetitions=count;return node;
    };
    const auto call=[&](const std::string& equation,Target target,const std::string& method) {
        auto node=scope(ExecutionKind::EquationCall,equation+".call");
        node.step={equation,std::move(target),node.id};
        bindings.push_back({equation,method,{},node.id});return node;
    };
    std::vector<ExecutionScope> mass,momentum,energy,closures,turbulence;
    for (auto node:program.root.children) {
        if (node.kind!=ExecutionKind::EquationCall)
            throw std::runtime_error("Unsupported Eulerian composition: non-equation default scope.");
        node.id=node.step.equation+".call";node.step.occurrence=node.id;node.order=0;
        const auto selected=std::find_if(bindings.begin(),bindings.end(),[&](const auto& b) {
            return b.equation==node.step.equation && b.occurrence.empty();
        });
        if (selected!=bindings.end() && selected->method=="EulerianTurbulenceClosure") closures.push_back(node);
        else if (selected!=bindings.end() && selected->method=="EulerianTurbulenceTransport") turbulence.push_back(node);
        else if (node.step.equation.rfind("E_CONTINUITY.",0)==0) mass.push_back(node);
        else if (node.step.equation.rfind("momentum.",0)==0) momentum.push_back(node);
        else if (node.step.equation.rfind("E_ENTHALPY.",0)==0) energy.push_back(node);
        else throw std::runtime_error("Unsupported Eulerian equation default: "+node.step.equation);
    }
    if (mass.size()<2 || mass.size()!=momentum.size() || mass.size()!=energy.size())
        throw std::runtime_error("Eulerian HOW requires complete phase transport groups.");
    const auto group=[&](const std::string& id,const std::vector<ExecutionScope>& members) {
        auto node=scope(ExecutionKind::Sequence,id);node.children=members;
        std::vector<std::string> targets;
        for (const auto& member:members) targets.push_back(member.step.target.symbol);
        for (auto& binding:bindings)
            for (const auto& member:members)
                if (binding.equation==member.step.equation && (binding.occurrence.empty() || binding.occurrence==member.step.occurrence))
                    binding.inputs=targets;
        return node;
    };
    auto outer=scope(ExecutionKind::Loop,"EE.outer",request.outerCorrectors);
    if (!closures.empty()) outer.children.push_back(group("phase.turbulence.closure",closures));
    outer.children.push_back(group("phase.continuity",mass));
    outer.children.push_back(group("phase.momentum",momentum));
    auto pressure=scope(ExecutionKind::Loop,"EE.pressure",request.pressureCorrectors);
    auto nonOrth=scope(ExecutionKind::Loop,"EE.nonOrthogonal",request.nonOrthogonalCorrectors+1);
    nonOrth.children.push_back(call("E_SHARED_PRESSURE",{"p",TargetKind::Correction},"EulerianSharedPressureCorrection"));
    std::vector<ExecutionScope> corrections,fluxes;
    for (const auto& item:momentum) {
        const auto suffix=item.step.equation.substr(std::string("momentum").size());
        corrections.push_back(call("correctPhase"+suffix,item.step.target,"EulerianPhaseCorrection"));
        fluxes.push_back(call("correctPhaseFlux"+suffix,{"volumeFaceFlux"+suffix,TargetKind::Workspace},"EulerianPhaseFluxCorrection"));
    }
    corrections.push_back(call("correctSharedP",{"p"},"EulerianPhaseCorrection"));
    nonOrth.children.push_back(group("phase.correction",corrections));
    nonOrth.children.push_back(group("phase.fluxCorrection",fluxes));
    pressure.children.push_back(std::move(nonOrth));outer.children.push_back(std::move(pressure));
    outer.children.push_back(group("phase.enthalpy",energy));
    if (!turbulence.empty()) outer.children.push_back(group("phase.turbulence.transport",turbulence));
    if (!closures.empty()) program.requirements.push_back({"EE.outer","phase.turbulence.closure",{}, {"phase.continuity","phase.momentum"}});
    if (!turbulence.empty()) program.requirements.push_back({"EE.outer","phase.turbulence.transport",{"phase.enthalpy"},{}});
    program.root=scope(ExecutionKind::Sequence,"EE.step");
    program.root.children.push_back(std::move(outer));
    program.root.children.push_back(scope(ExecutionKind::Commit,"EE.commit"));
}
}
