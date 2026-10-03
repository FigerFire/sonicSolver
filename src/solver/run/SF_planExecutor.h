#pragma once

/// @file SF_planExecutor.h
/// @brief Runtime-neutral traversal of structured compiled control flow.

#include "core/system/SF_solveProgram.h"

#include <functional>
#include <unordered_map>
#include <vector>

namespace SF::Run {

/// @brief Optional immutable values printed when SF_PLAN_TRACE is enabled.
struct PlanTraceContext {
    int step = 0;
    double time = 0.0;
    const double* dt = nullptr;
};

/// @brief Signals belong to one plan execution and never outlive its step.
class LoopSignalState {
public:
    void reset(const System::LoopSignalId& id);
    void publish(const System::LoopSignalId& id, bool value);
    [[nodiscard]] bool value(const System::LoopSignalId& id) const;

private:
    std::unordered_map<System::LoopSignalId,bool> values_;
};

/// @brief Minimal dynamic position within structured plan execution.
struct ExecutionContext {
    int stageIndex = -1;
    int stageCount = 0;
    LoopSignalState* signals = nullptr;
    const System::CompiledTarget* target = nullptr;
};

class OpRegistry {
public:
    using Operation = std::function<void()>;
    using ContextOperation = std::function<void(const ExecutionContext&)>;
    void bind(System::OpId id, Operation operation);
    void bind(System::OpId id, ContextOperation operation);
    void bind(System::OpId id,std::string provider,Operation operation);
    void bind(System::OpId id,std::string provider,ContextOperation operation);
    void invoke(const System::OpId& id,
                const ExecutionContext& context,std::string_view provider = {}) const;
    [[nodiscard]] bool contains(const System::OpId& id) const;
    [[nodiscard]] bool contains(const System::OpId& id,std::string_view provider) const;
    /// @brief 只保留编译期已分配给当前 runtime provider 的 OpId。
    void retain(const std::vector<System::OpId>& assigned);
    /// @brief 复用同一 registry 时丢弃上一轮的 binding。
    void clear() { entries_.clear(); }

private:
    struct Entry { System::OpId id; std::string provider; ContextOperation operation; };
    std::vector<Entry> entries_;
};

class PlanExecutor {
public:
    /// @brief Fail before execution if any Plan leaf lacks a provider.
    static void validateBindings(const System::CompiledSolvePlan& plan,
                                 const OpRegistry& operations);

    static void execute(const System::CompiledSolvePlan& plan,
                        const OpRegistry& operations,
                        const PlanTraceContext* trace = nullptr);

};

} // namespace SF::Run
