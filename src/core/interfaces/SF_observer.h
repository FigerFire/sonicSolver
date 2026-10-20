#pragma once

/// @file SF_observer.h
/// @brief 求解器结构化消息及观察者接口。

#include <functional>
#include <string>
#include <utility>

namespace SF::FDM {

/// @brief 求解器结构化消息类别。
enum class SolverMessageKind {
    Setup,
    StateClosure,
    TimeStep,
    Warning,
    Fatal
};

/// @brief 求解器向终端、GUI 或测试发送的值对象消息。
struct SolverMessage {
    SolverMessageKind kind = SolverMessageKind::Setup;
    std::string topic;
    std::string detail;
    int step = 0;
    double time = 0.0;
    double dt = 0.0;
};

/// @brief 求解器消息观察者接口。
class ISolverObserver {
public:
    virtual ~ISolverObserver() = default;
    virtual void onSolverMessage(const SolverMessage& message) = 0;
};

/// @brief 用可调用对象实现的观察者适配器。
class FunctionObserver final : public ISolverObserver {
public:
    explicit FunctionObserver(std::function<void(const SolverMessage&)> sink)
        : sink_(std::move(sink)) {}

    void onSolverMessage(const SolverMessage& message) override {
        if (sink_) sink_(message);
    }

private:
    std::function<void(const SolverMessage&)> sink_;
};

} // namespace SF::FDM
