#pragma once

/// @file SF_residual.h
/// @brief 已冻结残差的采样、范数、归一化和停止判据；不更新物理量。
#include "SF_planExecutor.h"
#include <map>

namespace SF::Run {
/// @brief 每个Loop调用自己的评价历史；Pending不是零残差或收敛。
class ResidualEvaluation {
public:
    using Revisions=std::map<System::OpId,std::size_t>;
    ResidualEvaluation(const System::SolvePlanNode&,const OpRegistry&,const ExecutionContext&);
    bool enter(const Revisions&);
    bool continueIteration(const Revisions& before,const Revisions& after,int iteration);
    void finish(int iterations,bool converged) const;
private:
    struct Metric {
        std::vector<double> previous;
        double initial=0;
        bool initialized=false;
    };
    const System::SolvePlanNode& loop_;
    const OpRegistry& operations_;
    const ExecutionContext& execution_;
    std::vector<Metric> metrics_;
    bool evaluate(std::size_t,int);
    bool sampleValid(std::size_t,const Revisions&) const;
};
}
