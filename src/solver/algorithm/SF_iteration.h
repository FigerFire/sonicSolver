#pragma once

/// @file SF_iteration.h
/// @brief Physical-time transient base and solution-relaxation arithmetic.

#include "core/state/SF_distributedField.h"
#include <vector>

namespace SF::SolverAlgorithm {

/// @brief 借用实际视图保存数值快照；无时钟、发布或数学刷新。
inline void snapshotValues(const State::DistributedFieldView& view,std::vector<double>& values) {
    values.resize((size_t)view.geometry->TotalSize()*view.components);
    for(int cell=0;cell<view.geometry->TotalSize();++cell)
        for(int c=0;c<view.components;++c)values[(size_t)c*view.geometry->TotalSize()+cell]=view.read(cell,c);
}


inline double fixedTimeCandidate(double physicalPrevious,
                                 double temporalIncrementAtIterate) {
    return physicalPrevious+temporalIncrementAtIterate;
}

inline double relaxSolution(double previousIterate,double candidate,
                            double factor) {
    if (factor==1.0) return candidate;
    return previousIterate+factor*(candidate-previousIterate);
}

} // namespace SF::SolverAlgorithm
