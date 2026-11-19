#pragma once
#include "SF_scalarBinding.h"
#include "solver/system/SF_cell.h"
namespace SF::SolverAlgorithm {
/// @brief 单次 HOW occurrence 的 cell 数值数组；借用物理 STATE/边界，不拥有 clock 或循环。
class CellEquationExecution {
    const System::CompiledEquationCall& call_;
    const System::CompiledCellEquation& contract_;
    const System::CompiledNumericalSystem& numerical_;
    State::StateBundle& state_;
    System::StateRealization& realized_;
    FDM::IExecutionRuntime* runtime_;
    Field& mesh_;
    struct FieldPort {
        State::DistributedFieldView* view;
        std::weak_ptr<const FDM::ScalarEquationData> data;
        std::vector<unsigned char> boundary;
        std::array<bool,6> sides{};
        unsigned long long prepared=~0ULL;double time=0;
    };
    std::vector<FieldPort> fields_;
    using ViewReader=std::function<State::DistributedFieldView&()>;
    struct SpatialTerm {
        const System::CompiledSpatialBinding* compiled;
        System::PrimitiveMomentumSource source;
        ViewReader coefficient;
        bool prepareCoefficient=false;
    };
    std::vector<SpatialTerm> terms_;
    std::vector<ViewReader> gradients_;
    State::DistributedFieldView workingView_;
    State::DistributedFieldView* physical_;
    std::vector<double> working_,candidate_,response_,rhs_,old_,stage_,withoutGradient_;
    std::array<std::vector<double>,4> increments_;
    State::DistributedFieldView oldView_,stageView_;
    int nextStage_=-1;double evaluationTime_=0;
    void evaluate(State::DistributedFieldView&);
    ViewReader bindInputView(const System::StateUse&);
    const std::vector<double>* base_=nullptr;
    std::vector<double> previousIterate_;
    int frozenStep_=-1,fixedIteration_=0;
    std::vector<int> cells_,candidates_;
    bool retainWorking_=false,ready_=false,assembled_=false,inputsBound_=false;
    unsigned long long revision_=0,preparedRevision_=~0ULL,predictionRevision_=0,potentialRevision_=0;
    EquationInterval interval_;
    int preparedStep_=-1,updated_=-1;
    unsigned long long generation_=0;
    double preparedTime_=0;
    void prepareField(FieldPort&,State::DistributedFieldView&);
public:
    CellEquationExecution(const System::CompiledEquationCall&,const System::CompiledNumericalSystem&,
        State::StateBundle&,const FDM::SolverServices&,System::StateRealization&,bool responseRequired);
    void snapshot();
    void prepareStage(const Run::ExecutionContext&);
    void evaluateStage(const Run::ExecutionContext&);
    void advanceStage(const Run::ExecutionContext&);
    void validatePublish();
    void publish();
    void bindInputs();
    void begin();
    double stepSize();
    void assemble();
    void advance();
    void accept() const;
    void prepareOutput(State::DistributedFieldView&,bool force=false);
    void prepareInput(std::string_view);
    std::shared_ptr<const FDM::ScalarEquationData> inputData(std::string_view) const;
    void copyToPhysical();
    void invalidateWorking() {ready_=false;}
    void freezeBase();
    void beginFixedIteration();
    bool hasFrozenBase() const {return frozenStep_==state_.step;}
    const std::vector<double>& previousIterate() const {return previousIterate_;}
    void requireResponse() const;
    void setResponse(double dt,const std::vector<int>* region=nullptr);
    const System::CompiledEquationCall& call() const {return call_;}
    const std::vector<double>& candidate() const {return candidate_;}
    const std::vector<double>& response() const {return response_;}
    std::vector<double>& response() {return response_;}
    unsigned long long generation() const {return generation_;}
    const std::vector<int>& cells() const {return cells_;}
    const std::vector<int>& candidates() const {return candidates_;}
    const std::vector<unsigned char>& boundary() const {return fields_.front().boundary;}
    const std::array<bool,6>& sides() const {return fields_.front().sides;}
    State::DistributedFieldView& physical() {return *physical_;}
    State::DistributedFieldView& output() {return call_.target.kind==System::TargetKind::Working?workingView_:*physical_;}
    State::DistributedFieldView& working() {return workingView_;}
};
/// @brief 接线一次 cell 调用，将确切实例供所声明的面响应借用。
BoundEquation bindCell(const EquationBindingInput&,bool responseRequired,bool fixedBaseRequired,
    std::shared_ptr<CellEquationExecution>&);
}
