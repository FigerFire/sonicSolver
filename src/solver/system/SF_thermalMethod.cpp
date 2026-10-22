#include "SF_thermalMethod.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <set>
#include <stdexcept>
namespace SF::System {
namespace {
class ThermalMethod final:public IProvider {
public:
    std::string_view id() const override {return "ThermalEulerTransport";}
    std::string_view runtimeProvider() const override {return ThermalProvider;}
    void lowerLifecycle(const ExecutionScope& scope,bool,const std::vector<CompiledEquationCall>&,SolvePlanNode& node) const override {
        if(scope.kind!=ExecutionKind::Commit)return;
        for(const auto* op:{"thermal.publish",OpIds::TimeCommit}) {
            SolvePlanNode leaf;leaf.kind=PlanNodeKind::Commit;leaf.id=leaf.operation=op;leaf.provider=ThermalProvider;
            node.children.push_back(std::move(leaf));
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding& binding) const override {
        const auto& state=system.state.at(call.target.symbol);
        if (call.target.symbol!="thermal.T" || call.target.kind!=TargetKind::Physical
            || state.components!=1 || state.componentOffset!=0 || state.readOnly()
            || state.storageBinding!=StorageBinding::NamedDistributed || state.storageKey.empty()
            || state.storageKey=="conservative" || !binding.inputs.empty())
            throw std::runtime_error("Thermal transport requires an independent named scalar STATE, not the EOS temperature/energy authority.");
        if(system.state.contains("T") && producedBy(system.state.at("T"),"closure.thermodynamic.T") && state.storageKey=="T")
            throw std::runtime_error("Independent thermal target conflicts with EOS T identity: choose an explicit distinct fieldName, e.g. thermalTemperature.");
        const std::set<std::string> names={"convection","diffusion","D","bounded","lower","upper","boundMode","ilw","ilwOrder"};
        if (binding.parameters.size()!=names.size()) throw std::runtime_error("Incomplete thermal numerical recipe.");
        for(const auto& p:binding.parameters)if(!names.count(p.first) || !std::isfinite(p.second))
            throw std::runtime_error("Invalid thermal numerical parameter: "+p.first);
        const auto& p=binding.parameters;
        for(const auto* key:{"convection","diffusion","bounded","ilw"})if(p.at(key)!=0 && p.at(key)!=1)
            throw std::runtime_error("Invalid thermal numerical flag: "+std::string(key));
        if (p.at("D")<0 || (p.at("diffusion") && p.at("D")<=0) || p.at("lower")>p.at("upper")
            || (p.at("boundMode")!=0 && p.at("boundMode")!=1 && p.at("boundMode")!=2))
            throw std::runtime_error("Invalid thermal diffusivity/bound recipe.");
        if (p.at("ilw") && p.at("ilwOrder")!=3 && p.at("ilwOrder")!=5 && p.at("ilwOrder")!=7 && p.at("ilwOrder")!=9)
            throw std::runtime_error("Thermal ILW requires order 3/5/7/9.");
        if(p.at("convection")) {
            const auto& u=system.state.at("U");const auto& rho=system.state.at("rho");const auto& momentum=system.state.at("rhoU");
            if(!producedBy(u,"closure.thermodynamicVelocity") || u.components!=3
                || rho.storageBinding!=StorageBinding::PackedDistributed || rho.storageKey!="conservative" || rho.componentOffset!=0
                || momentum.storageBinding!=StorageBinding::PackedDistributed || momentum.storageKey!="conservative"
                || momentum.componentOffset!=1 || momentum.components!=3)
                throw std::runtime_error("Unsupported thermal velocity port: requires U closure from five-variable packed rho/rhoU, not an unrelated STATE alias.");
        }
        using E=FormulaExpr;auto lhs=E::op(Math::Op::Ddt,{system.state.expression("thermal.T")});
        if(p.at("convection"))lhs=E::add(std::move(lhs),E::op(Math::Op::Div,{system.state.expression("U"),system.state.expression("thermal.T")}));
        auto rhs=p.at("diffusion")?E::op(Math::Op::Diffusion,{E::constantValue(p.at("D"),Math::signature(0,{0,2,-1,0,0,0,0})),system.state.expression("thermal.T")}):E::constantValue(0);
        if(canonicalFormula(system.registry.at(call.equation))!=canonicalFormula(Equation{call.equation,std::move(lhs),std::move(rhs)}))
            throw std::runtime_error("Unsupported thermal mathematics: numerical recipe does not match WHAT.");
        CompiledEquationCall r;r.source=call;r.equationMethod=id();r.target.symbol=call.target.symbol;r.target.kind=call.target.kind;
        r.providerContract=p;r.temporalResidual=r.synchronousStages=true;r.publishesStageToPhysicalTarget=false;
        r.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        r.backendOperation="thermal.stage";r.temporalStepSize="thermal.dt";
        r.temporalSnapshot="thermal.snapshot";r.temporalStagePrepare="thermal.boundary";
        r.temporalRhs="thermal.rhs";r.temporalAdvance="thermal.advance";r.temporalPublish="thermal.publish";
        r.oldTimeWorkspace="thermal.old";r.stageWorkspace="thermal.stage";r.residualWorkspace="thermal.rhs";
        r.requiredHaloWidth=1;r.numericalWorkspaces={r.oldTimeWorkspace,r.stageWorkspace,r.residualWorkspace};
        r.calls={{call.equation,call.target.symbol}};r.reads={call.target.symbol};r.writes={call.target.symbol};
        r.stateUses={{call.target.symbol,StateVersion::Stage}};
        if(p.at("convection")){r.reads.push_back("U");r.stateUses.push_back({"U",StateVersion::Current});}
        r.stateEffects={{r.residualWorkspace,true,false,false}};
        r.requirements={"serial single patch; Forward Euler only; independent temperature never changes EOS T or rhoE",
            "same Upwind1/Central metric kernel, BC and configured bounds as the former thermal update"};
        for(const auto& op:{r.temporalSnapshot,r.temporalStagePrepare,r.temporalRhs,r.temporalAdvance,r.temporalPublish})
            r.operations.push_back({op,op,op==r.temporalPublish?OperationStage::StepCommit:OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        return r;
    }
};
}
void addThermalMethod(ProviderRegistry& registry){static const ThermalMethod method;registry.add(method);}
}
