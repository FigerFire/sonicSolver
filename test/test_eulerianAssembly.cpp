#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_eulerianAssembly.h"
#include "solver/system/SF_eulerianRelations.h"
#include "solver/system/SF_pressureCoupling.h"
#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
ResolvedSimulationSystem fixture(int count) {
    SF::FDM::SolverConfig config;config.numerics.maxDeltaT=0.001;
    config.pressure.coupling.preset=SF::FDM::PressureCouplingPreset::PIMPLE;
    BuildRequest request;request.templateOrigin=PhysicsTemplateKind::EulerianEulerian;
    request.phaseNames={"liquid","gas"};if (count==3) request.phaseNames.push_back("third");
    request.referencePhase="gas";request.coupling=couplingRequestFrom(config.pressure.coupling,true);
    return build(config,request);
}
std::string signature(const CompiledExecutionProgram& program) {
    std::ostringstream out;
    for (const auto& call:program.steps) {
        const auto* c=eulerianAssemblyContract(call);
        require(c,"A native Eulerian occurrence lacks its typed assembly contract.");
        out<<c->equation<<':'<<c->method<<':'<<c->provider<<':'<<c->target.symbol
            <<':'<<static_cast<int>(c->target.kind)<<':'<<static_cast<int>(c->relation)<<':'<<c->phaseName;
        if (c->phaseSlot) out<<*c->phaseSlot;
        for (auto extension:c->extensions) out<<static_cast<int>(extension);
        out<<';';
    }
    return out.str();
}
void editExpr(FormulaExpr& expr,const std::function<bool(FormulaExpr&)>& edit) {
    if (edit(expr)) return;
    for (auto& child:expr.arguments) editExpr(child,edit);
}
std::vector<std::string> categories(const Equation& equation) {
    std::vector<std::string> result;
    const auto inventory=[&](const auto& self,const FormulaExpr& e)->void {
        if (e.kind==FormulaExpr::Kind::Operator) {
            result.push_back(e.name);
            return;
        }
        for (const auto& child:e.arguments) self(self,child);
    };
    inventory(inventory,equation.lhs);result.push_back("=");inventory(inventory,equation.rhs);
    return result;
}
}
int main() {
    using R=EulerianAssemblyRelation;using E=FormulaExpr;
    for (const int count:{2,3}) {
        const auto resolved=fixture(count);
        const auto& program=resolved.solvePlan.compiledProgram;
        require(resolved.rawSystem.legacyEquations.empty() && resolved.rawSystem.legacyDefinitions.equations().empty()
            && resolved.executableSystem.legacyEquations.empty() && resolved.executableSystem.legacyDefinitions.equations().empty(),
            "Native Eulerian retained duplicate legacy equations/definitions.");
        require(program.steps.size()==static_cast<std::size_t>(5*count+2),"Missing phase/shared contract.");
        int references=0;
        for (const auto& call:program.steps) {
            const auto* c=eulerianAssemblyContract(call);require(c,"Missing native contract.");
            require(c->provider==call.backendProvider && c->method==call.equationMethod
                && c->target.symbol==call.target.symbol,"Contract lost frozen implementation/target provenance.");
            if (c->phaseSlot) {
                const std::vector<std::string> names={"liquid","gas","third"};
                require(*c->phaseSlot<static_cast<std::size_t>(count) && c->phaseName==names[*c->phaseSlot],"Phase slot is not the declaration order.");
            }
            if (c->relation==R::ReferencePhaseClosure) {
                ++references;require(c->phaseSlot==1,"Reference closure lost its actual phase slot.");
            }
            if (c->relation==R::SharedPressureCorrection) require(!c->phaseSlot,"Shared pressure gained a phase slot.");
        }
        require(references==1,"Reference phase contract must be unique.");
        auto names=std::vector<std::string>{"liquid","gas"};if (count==3) names.push_back("third");
        validateEulerianAssemblyBindings(program,names,1);
        bool wrongRuntimeRejected=false;
        try { validateEulerianAssemblyBindings(program,names,0); } catch (const std::runtime_error&) {wrongRuntimeRejected=true;}
        require(wrongRuntimeRejected,"Runtime accepted a reference closure differing from its frozen phase slot.");
        auto absent=program;absent.steps.front().providerContract.reset();wrongRuntimeRejected=false;
        try {validateEulerianAssemblyBindings(absent,names,1);} catch (const std::runtime_error&) {wrongRuntimeRejected=true;}
        require(wrongRuntimeRejected,"Runtime rediscovered an implementation without its compiled contract.");
        auto changed=resolved.executableSystem;changed.state=StateRegistry{};
        StateSymbol unrelated;unrelated.id="unrelated";unrelated.storageKey="unrelated";changed.state.add(unrelated);
        const auto& symbols=resolved.executableSystem.state.symbols();
        for (auto it=symbols.rbegin();it!=symbols.rend();++it) changed.state.add(*it);
        auto compile=[&](const auto& equations) {return compileExecutionProgram(equations,resolved.solvePlan.sourceProgram,resolved.numericalSelection.bindings,builtinProviders());};
        require(signature(compile(changed))==signature(program),"Registry insertion changed frozen phase slots.");
        changed=resolved.executableSystem;changed.legacyDefinitions=SF::Equation::System{};
        require(signature(compile(changed))==signature(program),"Native compiled contract reads legacy Definition.");
        // The shadow legacy system is deliberately poisoned; AST-derived contracts stay identical.
        changed.legacyDefinitions.add(SF::Equation::named("momentum.liquid",SF::Equation::algebraic({"fake"})==SF::Equation::Symbol{"zero"}));
        require(signature(compile(changed))==signature(program),"Legacy mutation altered native compiled authority.");
        const auto mutation=[&](const std::string& id,const std::function<void(Equation&)>& edit,bool sameCategories) {
            auto equations=resolved.executableSystem;auto formula=equations.registry.at(id);
            const auto before=categories(formula);edit(formula);
            if (sameCategories) require(before==categories(formula),"Mutation unexpectedly changed flat term categories.");
            equations.registry.replace(formula);bool rejected=false;
            try { (void)compile(equations); } catch (const std::runtime_error& e) {rejected=std::string(e.what()).find("mathematical contract")!=std::string::npos;}
            require(rejected,"Malformed AST inherited a supported Eulerian implementation.");
        };
        mutation("momentum.liquid",[](auto& eq) {
            editExpr(eq.lhs,[](auto& e) {if (e.kind!=E::Kind::Multiply) return false;e=E::op("gradient",{E::multiply(E::symbol("alpha.liquid"),E::symbol("p"))});return true;});
        },true);
        mutation("momentum.liquid",[](auto& eq) {eq.rhs.arguments[0]=E::negate(eq.rhs.arguments[0]);},true);
        mutation("E_ENTHALPY.liquid",[](auto& eq) {auto source=eq.rhs.arguments[1].arguments[0];eq.rhs.arguments[1]=std::move(source);},true);
        mutation("E_SHARED_PRESSURE",[](auto& eq) {eq.lhs.arguments[0].arguments[0]=E::symbol("wrongCompressibility");},true);
        mutation("correctPhase.liquid",[](auto& eq) {editExpr(eq.rhs,[](auto& e) {if (e.kind!=E::Kind::Divide) return false;e.arguments[0]=E::symbol("wrongAlpha");return true;});},true);
        mutation("correctPhaseFlux.liquid",[](auto& eq) {eq.rhs.arguments.pop_back();},true);
        mutation("correctSharedP",[](auto& eq) {eq.rhs.arguments[1].arguments[0]=E::symbol("wrongRelaxation");},true);
        auto extensions=resolved.executableSystem;
        const auto append=[&](const std::string& id,const std::string& source) {
            auto formula=extensions.registry.at(id);formula.rhs=E::add(formula.rhs,E::op("source",{E::symbol(source)}));extensions.registry.replace(formula);
        };
        append("momentum.liquid","gravity");append("momentum.liquid","MRF");append("E_ENTHALPY.liquid","wallHeat");
        auto supported=compile(extensions);
        const auto momentum=std::find_if(supported.steps.begin(),supported.steps.end(),[](const auto& c) {return c.source.equation=="momentum.liquid";});
        require(eulerianAssemblyContract(*momentum)->extensions==std::vector<EulerianSourceExtension>{EulerianSourceExtension::Gravity,EulerianSourceExtension::MRF},"Supported source extensions were not frozen from AST.");
        append("momentum.liquid","unknownSource");bool rejected=false;
        try {(void)compile(extensions);} catch (const std::runtime_error&) {rejected=true;}
        require(rejected,"Unknown extension silently entered the source ledger.");
        extensions=resolved.executableSystem;append("momentum.liquid","MRF");append("momentum.liquid","MRF");rejected=false;
        try {(void)compile(extensions);} catch (const std::runtime_error&) {rejected=true;}
        require(rejected,"Duplicate source AST exceeded the implemented ledger contract.");
    }
}
