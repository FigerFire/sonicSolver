/// @file SF_time.cpp
/// @brief runTime/timeStep/maxDeltaT 的统一时间循环驱动。

#include "solver/algorithm/time/SF_time.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::Time {

namespace {

/// @brief 仅吸收时间累加产生的机器精度尾差，不改变用户给定的物理时间步。
double endTimeRoundoffTolerance(double time, double endTime) {
    return 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({std::abs(time), std::abs(endTime), 1.0});
}

} // namespace

Driver::Driver(RunControl control) : control_(control) {
    if (!std::isfinite(control_.startTime)
        || !std::isfinite(control_.endTime)
        || control_.endTime < control_.startTime
        || control_.endStep < 0) {
        throw std::runtime_error("Time::Driver received invalid runTime bounds.");
    }
    if (control_.writeByStep) {
        if (control_.writeIntervalSteps <= 0) {
            throw std::runtime_error(
                "timeStep output requires positive writeIntervalSteps.");
        }
    } else if (!std::isfinite(control_.writeIntervalTime)
               || control_.writeIntervalTime <= 0.0) {
        throw std::runtime_error(
            "runTime output requires positive writeIntervalTime.");
    }
}

void Driver::run(const DriverCallbacks& callbacks) const {
    if (!callbacks.advance) {
        throw std::runtime_error("Time::Driver requires an advance callback.");
    }
    double time = control_.startTime;
    int step = 0;
    int nextOutputStep = control_.writeIntervalSteps;
    double nextOutputTime = control_.startTime
                          + control_.writeIntervalTime;
    const auto globallyFinished = [&](bool local) {
        return callbacks.globallyFinished
            ? callbacks.globallyFinished(local) : local;
    };
    while (true) {
        const double remainingTime = control_.endTime - time;
        const bool reachedEndTime = remainingTime
            <= endTimeRoundoffTolerance(time, control_.endTime);
        const bool localFinished = reachedEndTime
            || (control_.endStep > 0 && step >= control_.endStep);
        if (globallyFinished(localFinished)) break;

        ++step;
        double limit = remainingTime;
        if (!control_.writeByStep) {
            limit = std::min(limit, nextOutputTime - time);
        }
        if (!std::isfinite(limit) || limit <= 0.0) {
            throw std::runtime_error(
                "Time::Driver produced a non-positive time-step limit.");
        }
        const AdvanceResult result = callbacks.advance(limit, step);
        const double tolerance = 1.0e-11
            * std::max({std::abs(time), std::abs(result.time), 1.0});
        if (!std::isfinite(result.dt) || result.dt <= 0.0
            || result.dt > limit + tolerance
            || !std::isfinite(result.time)
            || std::abs(result.time - (time + result.dt)) > tolerance) {
            throw std::runtime_error(
                "Time::Driver stepper returned inconsistent dt/time.");
        }
        time = result.time;
        if (callbacks.report) callbacks.report(step, result);
        if (control_.writeByStep && step == nextOutputStep) {
            if (callbacks.saveStep) callbacks.saveStep(step);
            nextOutputStep += control_.writeIntervalSteps;
        } else if (!control_.writeByStep
                   && time >= nextOutputTime - tolerance) {
            if (callbacks.saveTime) callbacks.saveTime(nextOutputTime);
            nextOutputTime += control_.writeIntervalTime;
        }
    }
}

} // namespace SF::Time
