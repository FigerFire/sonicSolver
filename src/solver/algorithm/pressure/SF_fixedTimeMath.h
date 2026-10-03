#pragma once

/// @file SF_fixedTimeMath.h
/// @brief Physical-time transient base and solution-relaxation arithmetic.

namespace SF::Pressure {

inline double fixedTimeCandidate(double physicalPrevious,
                                 double temporalIncrementAtIterate) {
    return physicalPrevious+temporalIncrementAtIterate;
}

inline double relaxSolution(double previousIterate,double candidate,
                            double factor) {
    if (factor==1.0) return candidate;
    return previousIterate+factor*(candidate-previousIterate);
}

} // namespace SF::Pressure
