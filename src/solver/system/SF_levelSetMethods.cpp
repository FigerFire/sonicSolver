#include "SF_levelSetMethods.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <set>

namespace SF::System {
namespace {
using E=FormulaExpr;
Equation supportedMathematics(std::string_view method) {
    if (method=="LevelSetAdvection") return {"",E::add(E::op("ddt",{E::symbol("phi")}),
        E::op("advect",{E::symbol("U"),E::symbol("phi")})),E::constantValue(0.0)};
    if (method=="LevelSetReference") return {"",E::symbol("phi0"),E::symbol("phi")};
    if (method=="LevelSetReinitialization") return {"",E::op("ddtau",{E::symbol("phi")}),
        E::negate(E::multiply(E::op("smoothedSign",{E::symbol("phi0")}),
            E::subtract(E::op("gradientNorm",{E::symbol("phi")}),E::constantValue(1.0))))};
    return {"",E::symbol("interfaceCurvature"),E::op("meanCurvature",{E::symbol("phi")})};
}
class LevelSetMethod final : public IProvider {
    std::string method_,target_,operation_;
public:
    LevelSetMethod(const char* method,const char* target,const char* operation)
        :method_(method),target_(target),operation_(operation) {}
    std::string_view id() const override { return method_; }
    std::string_view runtimeProvider() const override {
        return method_=="LevelSetAdvection" ? "flow.conservative" : "interface.level-set";
    }
    void lowerLifecycle(const ExecutionScope& scope,bool root,
            const std::vector<CompiledEquationCall>& calls,SolvePlanNode& lowered) const override {
        if (method_!="LevelSetGeometry" || root || scope.kind!=ExecutionKind::Sequence) return;
        const auto owns=[&](const auto& self,const ExecutionScope& node)->bool {
            if (node.kind==ExecutionKind::EquationCall)
                return std::any_of(calls.begin(),calls.end(),[&](const auto& call) {
                    return call.equationMethod==method_ && call.source.occurrence==node.step.occurrence;
                });
            return false;
        };
        if (std::none_of(scope.children.begin(),scope.children.end(),[&](const auto& child) {return owns(owns,child);})) return;
        SolvePlanNode pre;pre.kind=PlanNodeKind::Update;pre.id="flow.finalize.beforeInterface";
        pre.operation=OpIds::FlowFinalizeBegin;pre.provider="flow.conservative";
        lowered.children.insert(lowered.children.begin(),std::move(pre));
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding& binding) const override {
        auto expected=supportedMathematics(method_);expected.id=call.equation;
        if (canonicalFormula(system.registry.at(call.equation))!=canonicalFormula(expected)
            || call.target.symbol!=target_ || call.target.kind!=(method_=="LevelSetReference" ? TargetKind::Workspace : TargetKind::Physical))
            throw std::runtime_error("Unsupported level-set mathematical/output contract: "+method_);
        const std::set<std::string> parameters=method_=="LevelSetAdvection"
            ? std::set<std::string>{"order","epsilon","power","csf","ghostFluid","sigma","width"}
            : method_=="LevelSetReference" ? std::set<std::string>{"order","epsilon","power","pseudoDt","signFactor"}
            : std::set<std::string>{};
        for (const auto& value:binding.parameters) if (!parameters.count(value.first))
            throw std::runtime_error("Unknown "+method_+" numerical parameter: "+value.first);
        if (method_!="LevelSetReinitialization" && !binding.inputs.empty())
            throw std::runtime_error("Unsupported "+method_+" external numerical inputs.");
        for (const auto& value:binding.parameters)
            if (!std::isfinite(value.second) || value.second<0.0)
                throw std::runtime_error("Invalid level-set numerical parameter: "+value.first);
        if (method_=="LevelSetAdvection" || method_=="LevelSetReference") {
            const auto order=binding.parameters.at("order");
            if (order!=3 && order!=5 && order!=7) throw std::runtime_error("Level-set HJ-WENO supports only order 3/5/7.");
            if (binding.parameters.at("epsilon")<=0 || binding.parameters.at("power")<=0)
                throw std::runtime_error("Level-set WENO epsilon/power must be positive.");
            if (method_!="LevelSetAdvection") {
                if (binding.parameters.at("pseudoDt")<=0 || binding.parameters.at("signFactor")<=0)
                    throw std::runtime_error("Level-set pseudo-time step/sign smoothing must be positive.");
            }
        }
        if (method_=="LevelSetReinitialization"
            && (!binding.parameters.empty() || binding.inputs!=std::vector<std::string>{"levelSet.phi0"}))
            throw std::runtime_error("Pseudo-time stage must consume the prepared levelSet.phi0 workspace and its single frozen numerical recipe.");
        const auto& phi=system.state.at("phi");
        if (phi.components!=1 || phi.storageKey!="phi" || phi.storageBinding!=StorageBinding::NamedDistributed)
            throw std::runtime_error("Level-set requires the registered scalar phi authority.");
        CompiledEquationCall result;result.source=call;result.equationMethod=method_;
        result.temporalPreparation={OpIds::FlowStepPrepare};
        result.temporalStepSize=OpIds::FlowDtCompute;
        result.temporalSnapshot=OpIds::FlowStepBegin;
        result.target.symbol=target_;result.target.kind=call.target.kind;
        result.calls={{call.equation,target_}};result.writes={target_};result.reads={"phi"};
        result.stateUses={{"phi"}};
        if (method_=="LevelSetAdvection") {result.stateUses.push_back({"U"});result.stateEffects={{"phi",false,true}};}
        if (method_=="LevelSetReference") result.stateEffects={{"levelSet.phi0",true,false}};
        if (method_=="LevelSetReinitialization") {
            result.stateUses.push_back({"levelSet.phi0",StateVersion::Frozen});
            result.stateEffects={{"phi",false,true}};
        }
        if (method_=="LevelSetGeometry") result.stateEffects={{"interfaceNormal",true,true},{"interfaceCurvature",true,true}};
        result.backendOperation=operation_;result.providerContract=binding.parameters;
        if (method_=="LevelSetAdvection") {
            result.temporalResidual=true;result.publishesStageToPhysicalTarget=true;
            result.oldTimeWorkspace="explicit.scalar.phi.q0";
            result.reads.push_back("U");
            const auto csf=binding.parameters.at("csf"),ghost=binding.parameters.at("ghostFluid");
            if ((csf!=0 && csf!=1) || (ghost!=0 && ghost!=1) || (csf && ghost))
                throw std::runtime_error("Invalid frozen interface source/jump selection.");
            if (csf) {
                if (binding.parameters.at("width")<=0) throw std::runtime_error("CSF requires positive interface thickness.");
                auto force=E::multiply(E::multiply(E::constantValue(-binding.parameters.at("sigma")),E::symbol("interfaceCurvature")),
                    E::multiply(E::op("regularizedDelta",{E::symbol("phi"),E::constantValue(binding.parameters.at("width"))}),E::symbol("interfaceNormal")));
                Equation momentum{"interfaceMomentumSource",E::symbol("interfaceForce"),force};
                Equation energy{"interfaceEnergySource",E::symbol("interfacePower"),E::op("dot",{E::symbol("U"),force})};
                if (canonicalFormula(system.registry.at(momentum.id))!=canonicalFormula(momentum)
                    || canonicalFormula(system.registry.at(energy.id))!=canonicalFormula(energy))
                    throw std::runtime_error("Unsupported compiled CSF momentum/energy contribution.");
                result.calls.push_back({momentum.id,"rhoU"});result.calls.push_back({energy.id,"rhoE"});
                result.writes.push_back("rhoU");result.writes.push_back("rhoE");
                result.reads.push_back("interfaceNormal");result.reads.push_back("interfaceCurvature");
                result.operatorBindings.push_back("fused CSF residual source and U.dot(force) energy exchange");
            }
            if (ghost) {
                Equation jump{"interfacePressureJump",E::symbol("pLiquidMinusGas"),
                    E::multiply(E::constantValue(-binding.parameters.at("sigma")),E::symbol("interfaceCurvature"))};
                if (canonicalFormula(system.registry.at(jump.id))!=canonicalFormula(jump))
                    throw std::runtime_error("Unsupported frozen ghost-fluid pressure-jump relation.");
                result.requirements.push_back("declared jump query; a pressure operator must explicitly consume it");
            }
            result.operatorBindings={"HJ-WENO"+std::to_string((int)binding.parameters.at("order"))+" Hamilton-Jacobi advection"};
            result.requirements={"registered HamiltonJacobi scalar uses the same RK tableau as the conservative state",
                "phi BC -> halo -> geometry -> HJ-WENO RHS before stage publication"};
        } else {
            if (method_=="LevelSetReference") {
                result.target.viewOwner=StateViewOwner::NumericalProvider;
                result.target.workspace="levelSet.phi0";
                result.workspaceProvides={"levelSet.phi0"};
            }
            if (method_=="LevelSetReinitialization") {
                result.workspaceRequires={"levelSet.phi0"};result.reads.push_back("phi0");
                result.requirements={"frozen phi0 throughout pseudo-time loop; physical clock unchanged"};
            }
            result.operations.push_back({operation_,method_,OperationStage::StepCommit,
                {OperationCapability::LevelSetExecution},{OriginKind::Generated,method_}});
            if (method_=="LevelSetGeometry") {
                Equation normal{"interfaceNormal",E::symbol("interfaceNormal"),E::op("normalizedGradient",{E::symbol("phi")})};
                if (canonicalFormula(system.registry.at("interfaceNormal"))!=canonicalFormula(normal))
                    throw std::runtime_error("Unsupported interface normal mathematical contract.");
                const auto& view=system.state.at("interfaceNormal");
                if (view.components!=3 || view.storageKey!="levelSetNormal")
                    throw std::runtime_error("Interface geometry requires existing normal STATE view.");
                result.calls.push_back({"interfaceNormal","interfaceNormal"});
                result.writes.push_back("interfaceNormal");
            }
            if (method_=="LevelSetGeometry") result.operations.push_back({OpIds::FlowFinalizeBegin,
                "publish corrected flow before interface finalization",OperationStage::StepCommit,
                {OperationCapability::ConservativeExplicit},{OriginKind::Generated,method_}});
        }
        return result;
    }
};
}
void addLevelSetMethods(ProviderRegistry& providers) {
    static const LevelSetMethod advection("LevelSetAdvection","phi",OpIds::ExplicitStageExecute);
    static const LevelSetMethod reference("LevelSetReference","phi0",OpIds::LevelSetReference);
    static const LevelSetMethod reinit("LevelSetReinitialization","phi",OpIds::LevelSetPseudoStage);
    static const LevelSetMethod geometry("LevelSetGeometry","interfaceCurvature",OpIds::LevelSetGeometry);
    for (const auto* method:{&advection,&reference,&reinit,&geometry}) providers.add(*method);
}
}
