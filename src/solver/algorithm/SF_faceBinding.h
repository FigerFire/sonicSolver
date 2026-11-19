#pragma once
#include "scalar/SF_cellBinding.h"
#include "solver/discretization/pressure/SF_rhieChow.h"
namespace SF::SolverAlgorithm {
/// @brief 单份面数组及发布状态；借用显式 cell 响应，不拥有压力求解器或时间循环。
class FaceFluxExecution {
    const System::CompiledEquationCall& call_;
    const System::CompiledFaceFlux& contract_;
    std::shared_ptr<CellEquationExecution> predictor_;
    State::StateBundle& state_;
    FDM::IExecutionRuntime& runtime_;
    State::DistributedFieldView *pressure_,*physical_;
    State::DistributedFieldView predictedView_,correctedView_;
    std::vector<Pressure::RhieChow::Face> faces_;
    std::vector<double> predicted_,corrected_,response_;
    std::vector<unsigned char> pressureDirichlet_;
    unsigned long long consumed_=0,generation_=0,epochGeneration_=0,predictionGeneration_=0;
    EquationInterval interval_;
    bool predictedReady_=false,correctedReady_=false;
    void requireInterval() const;
public:
    FaceFluxExecution(const System::CompiledEquationCall&,std::shared_ptr<CellEquationExecution>,
        State::StateBundle&,FDM::IExecutionRuntime&,System::StateRealization&);
    void beginEpoch();
    void predict();
    void reconstruct(std::vector<double>&,const State::DistributedFieldView&);
    void invalidateCorrected() {correctedReady_=false;}
    void publishCorrected();
    void requirePrediction() const;
    double published(std::string_view,int axis,int lower) const;
    void synchronize(std::vector<double>&,int,State::ExchangeKind);
    CellEquationExecution& predictor() {return *predictor_;}
    const System::CompiledFaceFlux& contract() const {return contract_;}
    State::DistributedFieldView& pressure() {return *pressure_;}
    std::vector<double>& corrected() {return corrected_;}
    const std::vector<double>& predicted() const {return predicted_;}
    const std::vector<double>& response() const {return response_;}
    const std::vector<unsigned char>& pressureDirichlet() const {return pressureDirichlet_;}
};
}
