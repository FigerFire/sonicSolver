#pragma once

/// @file SF_runControl.h
/// @brief 所有求解工作流共用的运行与输出控制值对象。

namespace SF::Time {

/// @brief OpenFOAM 风格运行与输出控制的中立配置。
struct RunControl {
    double startTime = 0.0;
    double endTime = 0.0;
    int endStep = 0;
    bool writeByStep = false;
    int writeIntervalSteps = 0;
    double writeIntervalTime = 0.0;
};

} // namespace SF::Time
