#include "SF_multiphaseMethods.h"
#include "core/system/SF_operationIds.h"
#include <stdexcept>
namespace SF::System {
namespace {
class MultiphaseBalance final : public IProvider {
    std::string id_,provider_;
public:
    MultiphaseBalance(const char* id,const char* provider):id_(id),provider_(provider) {}
    std::string_view id() const override {return id_;}
    std::string_view runtimeProvider() const override {return provider_;}
    void lowerLifecycle(const ExecutionScope& scope,bool,
            const std::vector<CompiledEquationCall>&,SolvePlanNode& node) const override {
        if (scope.kind!=ExecutionKind::Commit || id_!="HomogeneousBalance") return;
        for (const auto* op:{OpIds::FlowStepCommit,OpIds::TimeCommit}) {
            SolvePlanNode leaf;leaf.kind=PlanNodeKind::Commit;leaf.id=leaf.operation=op;leaf.provider=provider_;
            node.children.push_back(std::move(leaf));
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding&) const override {
        const auto& formula=system.registry.at(call.equation);
        const auto& target=system.state.at(call.target.symbol);
        if (call.target.kind!=TargetKind::Physical || !dependsOn(formula.lhs,call.target.symbol)
            || (target.storageBinding!=StorageBinding::PackedDistributed && target.storageBinding!=StorageBinding::NamedDistributed))
            throw std::runtime_error("Invalid native multiphase balance/state contract.");
        CompiledEquationCall result;result.source=call;result.equationMethod=id_;
        result.temporalPreparation={OpIds::FlowStepPrepare};
        result.temporalStepSize=OpIds::FlowDtCompute;
        result.temporalSnapshot=OpIds::FlowStepBegin;
        result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.calls={{call.equation,call.target.symbol}};result.writes={call.target.symbol};
        if (id_=="HomogeneousBalance") {
            result.temporalResidual=true;result.publishesStageToPhysicalTarget=true;
            result.oldTimeWorkspace="explicit.homogeneous.q0";
            result.backendOperation=OpIds::ExplicitStageExecute;
            result.operatorBindings={"generic component flux implementation pending"};
            for (const auto* op:{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,OpIds::FlowStepBegin,
                                OpIds::ExplicitStageExecute,OpIds::FlowStepCommit,OpIds::TimeCommit})
                result.operations.push_back({op,"homogeneous component execution",OperationStage::Prepare,
                    {OperationCapability::ConservativeExplicit},{OriginKind::Generated,id_}});
        } else {
            result.temporalResidual=true;result.publishesStageToPhysicalTarget=true;
            result.oldTimeWorkspace="explicit.scalar."+call.target.symbol+".q0";
            result.backendOperation=OpIds::ExplicitStageExecute;
            result.operatorBindings={"scalar transport kernel; native thermodynamic/closure binding pending"};
        }
        return result;
    }
};
}
void addMultiphaseMethods(ProviderRegistry& providers) {
    static const MultiphaseBalance homogeneous("HomogeneousBalance","flow.homogeneous");
    static const MultiphaseBalance mixture("MixtureBalance","flow.conservative");
    providers.add(homogeneous);providers.add(mixture);
}
}
