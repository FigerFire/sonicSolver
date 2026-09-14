#pragma once

/// @file SF_time.h
/// @brief 所有方程系统共享的 runTime/timeStep/endTime 驱动器。

#include "core/config/SF_runControl.h"

#include <functional>
#include <string>

namespace SF::Time {

/// @brief 方程专用 stepper 返回给统一时间驱动器的推进结果。
struct AdvanceResult {
    double dt = 0.0;
    double time = 0.0;
    std::string detail;
};

/// @brief IO、日志和 MPI 终止归约回调。
struct DriverCallbacks {
    std::function<AdvanceResult(double, int)> advance;
    std::function<void(int)> saveStep;
    std::function<void(double)> saveTime;
    std::function<void(int, const AdvanceResult&)> report;
    std::function<bool(bool)> globallyFinished;
};

/// @brief 统一处理时间上限、输出对齐、endStep 和跨 rank 终止。
class Driver {
public:
    explicit Driver(RunControl control);
    void run(const DriverCallbacks& callbacks) const;

private:
    RunControl control_;
};

} // namespace SF::Time
