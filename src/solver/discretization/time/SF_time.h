#pragma once

/// @file SF_time.h
/// @brief 瞬态项离散调度入口。

#include "SF_config.h"
#include "SF_scalarField.h"
#include "methods/numerics/time/SF_time.h"
#include "core/state/SF_state.h"

#include <functional>
#include <stdexcept>
#include <vector>

namespace SF {

/// @brief 按强类型配置离散 `ddt`，并推进整个显式耦合方程组。
inline void ddtDispatch(
        Field& field,
        Residual& residual,
        double timeStep,
        FDM::TimeScheme scheme,
        const std::function<void(Field&)>& assemble,
        const std::function<void(Field&, const char*)>& postStage = {},
        ScalarField* auxiliary = nullptr,
        std::vector<double>* auxiliaryRHS = nullptr) {
    if (scheme == FDM::TimeScheme::Euler) {
        Euler::ddt(field, residual, timeStep, assemble, postStage,
                   auxiliary, auxiliaryRHS);
        return;
    }
    if (scheme == FDM::TimeScheme::SSPRK3) {
        if (auxiliary || auxiliaryRHS) {
            throw std::runtime_error(
                "SSPRK3 auxiliary scalars require VariableRegistry.");
        }
        SSPRK3::ddt(field, residual, timeStep, assemble, postStage, nullptr);
        return;
    }
    RK4::ddt(field, residual, timeStep, assemble, postStage,
             auxiliary, auxiliaryRHS);
}

/// @brief 使用 VariableRegistry 同步推进任意数量辅助未知量。
inline void ddtDispatch(
        Field& field,
        Residual& residual,
        double timeStep,
        FDM::TimeScheme scheme,
        const std::function<void(Field&)>& assemble,
        const std::function<void(Field&, const char*)>& postStage,
        State::VariableRegistry* variables) {
    if (scheme == FDM::TimeScheme::Euler) {
        Euler::ddt(field, residual, timeStep, assemble, postStage, variables);
    } else if (scheme == FDM::TimeScheme::SSPRK3) {
        SSPRK3::ddt(field, residual, timeStep, assemble, postStage, variables);
    } else {
        RK4::ddt(field, residual, timeStep, assemble, postStage, variables);
    }
}

/// @brief 按时间格式把每次 RHS 求值对应的真实 stage 时间传给调用者。
///
/// Euler 的节点为 `{0}`，SSP-RK3 的 Shu--Osher RHS 节点为
/// `{0,1,1/2}`，经典 RK4 的节点为 `{0,1/2,1/2,1}`。
inline void ddtDispatch(
        Field& field,
        Residual& residual,
        double startTime,
        double timeStep,
        FDM::TimeScheme scheme,
        const std::function<void(Field&, double)>& assemble,
        const std::function<void(Field&, const char*)>& postStage,
        State::VariableRegistry* variables) {
    int evaluation = 0;
    const auto stageFraction = [&](int index) {
        if (scheme == FDM::TimeScheme::Euler) {
            if (index == 0) return 0.0;
        } else if (scheme == FDM::TimeScheme::SSPRK3) {
            constexpr double nodes[] = {0.0, 1.0, 0.5};
            if (index >= 0 && index < 3) return nodes[index];
        } else {
            constexpr double nodes[] = {0.0, 0.5, 0.5, 1.0};
            if (index >= 0 && index < 4) return nodes[index];
        }
        throw std::runtime_error(
            "Time integrator requested an unexpected RHS evaluation.");
    };
    const std::function<void(Field&)> staged = [&](Field& state) {
        const double time = startTime
            + stageFraction(evaluation++) * timeStep;
        assemble(state, time);
    };
    ddtDispatch(field, residual, timeStep, scheme, staged, postStage, variables);
}

} // namespace SF
