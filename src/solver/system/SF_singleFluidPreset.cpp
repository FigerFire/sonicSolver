/// @file SF_singleFluidPreset.cpp
/// @brief 可读的内置守恒 NS 声明及显式默认调用；不执行数值操作。
#include "SF_singleFluidPreset.h"
#include "SF_builtinState.h"
#include "SF_presets.h"
#include "SF_pressureCoupling.h"
#include <algorithm>

namespace SF::System::Preset {
EquationCompositionConfig singleFluidDefaults() {
    EquationCompositionConfig defaults;
    defaults.stateDeclared=true;
    defaults.solutionVariables={"rho","rhoU","rhoE"};
    defaults.equations={"Continuity","Momentum","Energy"};
    defaults.algorithm="Explicit";
    return defaults;
}

SystemContribution singleFluidEquations(const SingleFluidPresetSpec& spec) {
    const BuiltinStateCatalog symbols;
    const auto transient=[&](const char* symbol,const char* occurrence) {
        return FormulaExpr::op(Math::Op::Ddt,{symbols.expression(symbol)},occurrence);
    };
    const auto divergence=[&](const char* symbol,const char* occurrence) {
        return FormulaExpr::op(Math::Op::Div,{symbols.expression(symbol)},occurrence);
    };
    const auto diffusion=[&](const char* coefficient,const char* field,const char* occurrence) {
        return FormulaExpr::op(Math::Op::Diffusion,
            {symbols.expression(coefficient),symbols.expression(field)},occurrence);
    };
    const auto mass=transient("rho","mass.time")+divergence("massFlux","mass.convection");
    auto momentum=transient("rhoU","momentum.time")+divergence("momentumFlux","momentum.convection");
    auto energy=transient("rhoE","energy.time")+divergence("energyFlux","energy.convection");
    if(spec.diffusion) {
        // Preserve (ddt + div) + diffusion, including its original residual sign.
        momentum=momentum+diffusion("mu","U","momentum.diffusion");
        energy=energy+diffusion("conductivity","T","energy.diffusion");
    }
    const Provenance origin{OriginKind::BuiltinPreset,"builtin conservativeNS"};
    SystemContribution result;
    result.recordContribution("builtin.conservativeNS","single-fluid conservative balances");
    for(const auto* id:{"U","p","T"})result.requireState(id);
    if(spec.diffusion)for(const auto* id:{"mu","conductivity"})result.requireState(id);
    result.addEquation({"continuity",mass,FormulaExpr::constantValue(0),origin});
    result.addEquation({"momentum",momentum,FormulaExpr::constantValue(0),origin});
    result.addEquation({"energy",energy,FormulaExpr::constantValue(0),origin});
    return result;
}

void installSingleFluid(SystemCompositionBuilder& system,const SingleFluidPresetSpec& spec,
                       const EquationCompositionConfig& composition) {
    auto balances=singleFluidEquations(spec);
    const auto selected=[&](const char* id) {
        return composition.equations.empty() || std::find(composition.equations.begin(),
            composition.equations.end(),id)!=composition.equations.end();
    };
    const char* choices[]={"Continuity","Momentum","Energy"};
    for(std::size_t i=balances.registeredEquations.size();i-- > 0;) {
        const auto& id=balances.registeredEquations[i].id;
        if(std::count(composition.equations.begin(),composition.equations.end(),choices[i])
            +std::count(composition.equations.begin(),composition.equations.end(),id)>1)
            throw std::runtime_error("Duplicate WHAT reference to built-in equation: "+id);
        if(!selected(choices[i]) && !selected(id.c_str()))
            balances.registeredEquations.erase(balances.registeredEquations.begin()+i);
    }
    if(balances.registeredEquations.empty())
        throw std::runtime_error("Selected WHAT has no equation from the built-in conservative NS package.");
    system.applyContribution(std::move(balances));
    system.recordContribution("builtin.conservativeExplicit","conservative equations + explicit execution");
    std::vector<ExecutionScope> calls;
    const auto addCall=[&](const char* equation,const char* target,Order order)->std::optional<EquationCall> {
        if(!system.rawSystem().registry.contains(equation))return std::nullopt;
        if(!system.rawSystem().state.isSolution(target))
            throw std::runtime_error(std::string("No selected solution STATE target for equation ")+equation+": "+target);
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;
        call.step={equation,{target}};call.order=order;
        call.origin={OriginKind::BuiltinPreset,"builtin conservativeExplicit"};
        const auto reference=call.step;
        calls.push_back(std::move(call));
        NumericalBinding method{equation,"ConservativeResidual"};
        method.origin={OriginKind::BuiltinPreset,"builtin conservativeExplicit"};
        system.bindNumerics(std::move(method));
        return reference;
    };
    addCall("continuity","rho",10);
    if(auto predictor=addCall("momentum","rhoU",20))
        system.pressureInputs=pressureCouplingInputs(*predictor,PressurePredictorKind::ConservativeMomentum);
    addCall("energy","rhoE",60);
    for(const auto& equation:system.rawSystem().registry.entries())
        if(std::none_of(calls.begin(),calls.end(),[&](const auto& call){return call.step.equation==equation.id;}))
            throw std::runtime_error("Selected equation has no built-in default call: "+equation.id);
    auto execution=Compose::explicitProgram(std::move(calls));
    for(auto& node:execution.root.children)system.addExecution(std::move(node));
}
} // namespace SF::System::Preset
