#pragma once

/// @file SF_equationState.h
/// @brief FluidStateModel/EOS-aware 状态合法性判断。

#include "SF_fluidStateModel.h"
#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"

#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SF::Numerics {

struct EquationStateSummary {
    int checkedCells = 0;
    double minDensity = std::numeric_limits<double>::max();
    double minPressure = std::numeric_limits<double>::max();
    double minTemperature = std::numeric_limits<double>::max();
    double minSoundSpeed = std::numeric_limits<double>::max();
};

inline EquationStateSummary validateEquationState(
        const Field& field,
        const Physics::FluidStateModel::Model& equations,
        const char* context) {
    if (field.NVar() != equations.variableCount())
        throw std::runtime_error("Equation-state validation variable-count mismatch.");
    EquationStateSummary summary;
    std::vector<double> q((size_t)field.NVar(), 0.0);
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) q[(size_t)v] = field(i,j,k,v);
        try {
            const auto state = equations.close(q.data(), field.NVar());
            if (!std::isfinite(state.density) || state.density <= 0.0
                || !std::isfinite(state.pressure) || state.pressure <= 0.0
                || !std::isfinite(state.temperature) || state.temperature <= 0.0
                || !std::isfinite(state.soundSpeed) || state.soundSpeed <= 0.0)
                throw std::runtime_error("EOS returned a non-positive derived state.");
            ++summary.checkedCells;
            summary.minDensity = std::min(summary.minDensity, state.density);
            summary.minPressure = std::min(summary.minPressure, state.pressure);
            summary.minTemperature = std::min(summary.minTemperature, state.temperature);
            summary.minSoundSpeed = std::min(summary.minSoundSpeed, state.soundSpeed);
        } catch (const std::exception& error) {
            std::ostringstream message;
            message << context << ": invalid FluidStateModel state at ("
                    << i << "," << j << "," << k << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    });
    if (summary.checkedCells == 0)
        throw std::runtime_error(std::string(context) + ": no solved fluid cells.");
    return summary;
}

} // namespace SF::Numerics
