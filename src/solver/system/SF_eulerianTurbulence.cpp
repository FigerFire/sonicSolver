#include "SF_eulerianTurbulence.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>

namespace SF::System {
namespace {
// Capability signatures of this frozen numerical implementation, never registered WHAT.
Equation supportedClosureMathematics(const std::string& model,
        const std::string& phase) {
    using E=FormulaExpr;
    return {"mu_t."+phase,E::symbol("mu_t."+phase),
        E::symbol(model+".eddyDynamicViscosity."+phase),
        {OriginKind::Model,model},true};
}
Equation supportedTransportMathematics(const std::string& second,
        const std::string& phase,const std::string& variable) {
    using E=FormulaExpr;
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
        {OriginKind::Model,model},true};
}
SolvePlanNode op(const char* id,const char* owner) {
    SolvePlanNode node;node.kind=PlanNodeKind::Update;node.id=id;node.operation=id;node.provider=owner;
    return node;
}
class EulerianTurbulenceMethod final:public IProvider {
public:
    explicit EulerianTurbulenceMethod(bool transport):transport_(transport) {}
    std::string_view id() const override {return transport_ ? "EulerianTurbulenceTransport" : "EulerianTurbulenceClosure";}
    std::string_view runtimeProvider() const override {return "flow.eulerian-turbulence";}
    std::string_view operationProvider(const OpId& operation) const override {
        return operation==OpIds::EeInterphaseCompute || operation==OpIds::EeSourcesAssemble
            ? "flow.eulerian-pressure" : runtimeProvider();
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& step,
            const NumericalBinding& binding) const override {
        if (!hasConstraint(system,"C_SHARED_PRESSURE") || !hasConstraint(system,"C_VOLUME_FRACTION")
            || binding.inputs.empty()) throw std::runtime_error("Eulerian turbulence requires explicit phase/model grouping and shared-pressure state.");
        CompiledEulerianTurbulenceContract contract;
        contract.equation=step.equation;contract.target=step.target;contract.method=std::string(id());
        contract.transport=transport_;
        std::string model,second;
        std::vector<CompiledMathRef> members;
        for (std::size_t i=0;i<binding.inputs.size();) {
            const auto& state=system.state.at(binding.inputs[i]);
            const auto& phase=state.nameSpace;
            const auto& closure=system.registry.at("mu_t."+phase);
            if (closure.rhs.kind!=FormulaExpr::Kind::Symbol)
                throw std::runtime_error("Unsupported Eulerian turbulence closure AST.");
            std::string selected;
            for (const auto* candidate:{"kEpsilon","kOmegaSST","Smagorinsky"})
                if (canonicalFormula(closure)==canonicalFormula(supportedClosureMathematics(candidate,phase))) selected=candidate;
            if (selected.empty() || (!model.empty() && selected!=model) || (transport_ && selected=="Smagorinsky"))
                throw std::runtime_error("Eulerian turbulence requires a consistent supported mathematical model.");
            model=selected;second=model=="kEpsilon" ? "epsilon" : "omega";
            contract.model=model=="kEpsilon" ? FDM::TurbulenceModelKind::kEpsilon
                : model=="kOmegaSST" ? FDM::TurbulenceModelKind::kOmegaSST : FDM::TurbulenceModelKind::Smagorinsky;
            const auto& mass=system.state.at("phaseMass."+phase);
            const auto end=mass.storageKey.find('.');
            if (mass.storageKey.rfind("phase",0)!=0 || end==std::string::npos || mass.storageKey.substr(end)!=".mass")
                throw std::runtime_error("Invalid Eulerian turbulence phase backing.");
            const auto slot=std::stoul(mass.storageKey.substr(5,end-5));
            const auto prefix="phase"+std::to_string(slot)+".";
            if (std::any_of(contract.phases.begin(),contract.phases.end(),[&](const auto& item) {return item.name==phase || item.slot==slot;}))
                throw std::runtime_error("Duplicate Eulerian turbulence phase group.");
            contract.phases.push_back({phase,slot});
            const auto checkState=[&](const std::string& variable,bool derived) {
                const auto& value=system.state.at(variable+"."+phase);
                if (value.components!=1 || value.componentOffset!=0 || value.nameSpace!=phase
                    || value.storageBinding!=StorageBinding::ProviderDistributed || value.storageKey!=prefix+variable
                    || value.role!=(derived ? StateRole::Derived : StateRole::Transported))
                    throw std::runtime_error("Eulerian turbulence STATE must alias its declared phase/model backing.");
            };
            checkState("mu_t",true);
            const auto variables=transport_ ? std::vector<std::string>{"k",second} : std::vector<std::string>{"mu_t"};
            for (const auto& variable:variables) {
                const auto name=variable+"."+phase;
                if (i>=binding.inputs.size() || binding.inputs[i++]!=name)
                    throw std::runtime_error("Eulerian turbulence requires ordered complete phase equation pairs.");
                checkState(variable,!transport_);
                if (transport_ && canonicalFormula(system.registry.at(name))
                        !=canonicalFormula(supportedTransportMathematics(second,phase,variable)))
                    throw std::runtime_error("Unsupported Eulerian turbulence mathematical contract: "+name);
                members.push_back({name,name});
            }
        }
        if (step.target.kind!=TargetKind::Physical || step.target.symbol!=step.equation
            || std::none_of(members.begin(),members.end(),[&](const auto& item) {return item.equation==step.equation;}))
            throw std::runtime_error("Eulerian turbulence requires each grouped equation's own physical target.");
        CompiledEquationCall result;result.source=step;result.target.symbol=step.target.symbol;
        result.target.kind=step.target.kind;result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.equationMethod=std::string(id());result.calls={{step.equation,step.target.symbol}};
        result.fusionKey=std::string(id())+"/"+model;result.fusionMembers=members;
        result.reads={"p"};
        for (const auto& phase:contract.phases) {
            for (const auto* name:{"phaseMass","rho","alpha","U","T","mu_t"}) result.reads.push_back(std::string(name)+"."+phase.name);
            if (model!="Smagorinsky") {result.reads.push_back("k."+phase.name);result.reads.push_back(second+"."+phase.name);}
            result.writes.push_back("mu_t."+phase.name);
        }
        if (transport_) for (const auto& member:members) {
            result.writes.push_back(member.target);result.stateEffects.push_back({member.target,false,false});
        }
        for (const auto& phase:contract.phases) {
            const auto coefficient="turbulence.coefficients."+phase.name;
            result.stateEffects.push_back({"mu_t."+phase.name,true,false});
            if (transport_) {
                result.stateUses.push_back({coefficient,StateVersion::Frozen,false,{},true});
                result.stateUses.push_back({"phaseMass."+phase.name});
            } else {
                for (const auto* name:{"rho","U","T"}) result.stateUses.push_back({std::string(name)+"."+phase.name});
                if (model!="Smagorinsky") {
                    result.stateUses.push_back({"k."+phase.name});result.stateUses.push_back({second+"."+phase.name});
                }
                result.stateEffects.push_back({coefficient,true,false});
            }
        }
        result.fragment.kind=PlanNodeKind::Sequence;result.fragment.id=std::string(id());
        const auto add=[&](const char* operation,const char* owner,OperationCapability capability) {
            auto leaf=op(operation,owner);leaf.equationCalls=result.calls;result.fragment.children.push_back(std::move(leaf));
            result.operations.push_back({operation,operation,OperationStage::Prepare,{capability},{OriginKind::Generated,std::string(id())}});
        };
        if (!transport_) {
            add(OpIds::EeInterphaseCompute,"flow.eulerian-pressure",OperationCapability::EulerianPhaseExecution);
            add(OpIds::EeSourcesAssemble,"flow.eulerian-pressure",OperationCapability::EulerianPhaseExecution);
        }
        add(transport_ ? OpIds::EeTurbulenceSolve : OpIds::EeTurbulencePrepare,
            "flow.eulerian-turbulence",OperationCapability::EulerianTurbulenceExecution);
        result.temporalMethod=transport_ ? "phase-mass weighted implicit transport; fixed physical old time" : "pre-predictor algebraic closure refresh";
        result.requirements={"complete ordered selected phase/model group","existing Upwind/diffusion/matrix and canonical mass flux","pre-predictor model coefficients stay lagged through transport; post-solve prepare refreshes caches",
            "existing boundary/prepare/commit semantics"};
        result.providerContract=std::make_shared<const CompiledEulerianTurbulenceContract>(std::move(contract));
        return result;
    }
    void lowerLifecycle(const ExecutionScope& source,bool root,const std::vector<CompiledEquationCall>& calls,SolvePlanNode& result) const override {
        if (transport_) return;
        if (root) {
            // A transported model cannot silently become only a cached viscosity closure.
            std::vector<const CompiledEquationCall*> closures,transport;
            for (const auto& call:calls) if (const auto* c=eulerianTurbulenceContract(call))
                (c->transport ? transport : closures).push_back(&call);
            if (closures.empty()) throw std::runtime_error("Missing Eulerian turbulence closure group.");
            const auto* c=eulerianTurbulenceContract(*closures.front());
            if (c->model==FDM::TurbulenceModelKind::Smagorinsky) {
                if (!transport.empty()) throw std::runtime_error("Algebraic LES cannot inherit a transported RAS kernel.");
            } else {
                if (transport.size()!=2*closures.size()) throw std::runtime_error("Eulerian RAS requires complete native transport occurrences.");
                const auto* t=eulerianTurbulenceContract(*transport.front());
                if (c->model!=t->model || c->phases.size()!=t->phases.size()) throw std::runtime_error("Inconsistent Eulerian closure/transport group.");
                for (std::size_t i=0;i<c->phases.size();++i)
                    if (c->phases[i].slot!=t->phases[i].slot || c->phases[i].name!=t->phases[i].name)
                        throw std::runtime_error("Eulerian closure/transport phase order differs.");
            }
        }
        if (source.kind!=ExecutionKind::EquationCall) return;
        const auto member=std::find_if(calls.begin(),calls.end(),[&](const auto& call) {return call.source.occurrence==source.step.occurrence;});
        if (member==calls.end() || member->equationMethod!="EulerianPhaseContinuity") return;
        // The explicit closure HOW group now owns the original source prerequisites.
        // Remove their duplicate realization from continuity, without moving a source sibling.
        if (result.children.size()<2 || result.children[0].operation!=OpIds::EeInterphaseCompute
            || result.children[1].operation!=OpIds::EeSourcesAssemble)
            throw std::runtime_error("Unsupported Eulerian pre-closure prerequisite topology.");
        result.children.erase(result.children.begin(),result.children.begin()+2);
    }
private:bool transport_;
};
}
void addEulerianTurbulenceMethods(ProviderRegistry& registry) {
    static const EulerianTurbulenceMethod closure(false),transport(true);
    registry.add(closure);registry.add(transport);
}
void validateEulerianTurbulenceBindings(const CompiledExecutionProgram& program,
        const std::vector<std::string>& names,const FDM::TurbulenceConfig& config) {
    const bool active=config.enabled && config.family!=FDM::TurbulenceFamily::None && config.family!=FDM::TurbulenceFamily::DNS;
    const bool ras=active && config.family==FDM::TurbulenceFamily::RAS;
    std::size_t closures=0,transport=0;
    for (const auto& call:program.steps) if (const auto* c=eulerianTurbulenceContract(call)) {
        if (!active || c->model!=config.model || c->phases.size()!=config.phaseNames.size()
            || call.backendProvider!="flow.eulerian-turbulence"
            || c->equation!=call.source.equation || c->target.symbol!=call.target.symbol
            || c->target.kind!=call.target.kind || c->method!=call.equationMethod
            || c->method!=(c->transport ? "EulerianTurbulenceTransport" : "EulerianTurbulenceClosure"))
            throw std::runtime_error("Frozen Eulerian turbulence implementation differs from active model.");
        for (std::size_t i=0;i<c->phases.size();++i) {
            const auto& phase=c->phases[i];
            if (phase.slot>=names.size() || phase.name!=names[phase.slot] || phase.name!=config.phaseNames[i])
                throw std::runtime_error("Frozen Eulerian turbulence phase order differs from actual storage.");
        }
        (c->transport ? transport : closures)++;
    }
    if (closures!=(active ? config.phaseNames.size() : 0)
        || transport!=(ras ? 2*config.phaseNames.size() : 0))
        throw std::runtime_error("Eulerian turbulence runtime lacks complete compiled closure/transport contracts.");
}
}
