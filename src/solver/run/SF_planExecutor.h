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
    using ResidualOperation = std::function<std::vector<double>(const ExecutionContext&)>;
    using ResidualAvailability = std::function<bool(const ExecutionContext&)>;
    void bind(System::OpId id, Operation operation);
    void bind(System::OpId id, ContextOperation operation);
    void bind(System::OpId id,std::string provider,Operation operation);
    void bind(System::OpId id,std::string provider,ContextOperation operation);
    void invoke(const System::OpId& id,
                const ExecutionContext& context,std::string_view provider = {}) const;
    [[nodiscard]] bool contains(const System::OpId& id) const;
    [[nodiscard]] bool contains(const System::OpId& id,std::string_view provider) const;
    /// @brief 借用已绑定方法的真实状态/残差；范数、归一化和继续迭代判断由残差模块负责。
    void bindResidual(std::string id,std::string provider,ResidualOperation operation,ResidualAvailability available={});
    [[nodiscard]] bool residualAvailable(std::string_view id,std::string_view provider,const ExecutionContext&) const;
    [[nodiscard]] bool containsResidual(std::string_view id,std::string_view provider) const;
    [[nodiscard]] std::vector<double> sampleResidual(std::string_view id,
        std::string_view provider,const ExecutionContext&) const;
    /// @brief 只保留编译期已分配给当前 runtime provider 的 OpId。
    void retain(const std::vector<System::OpId>& assigned);
    /// @brief 复用同一 registry 时丢弃上一轮的 binding。
    void clear();

private:
    friend class PlanExecutor;
    /// @brief Single-fluid 持有对同一不可变计划的借用；封存后禁止重新绑定。
    const System::CompiledSolvePlan* sealedPlan_ = nullptr;
    void requireMutable() const;
    struct Entry { System::OpId id; std::string provider; ContextOperation operation; };
    std::vector<Entry> entries_;
    struct ResidualEntry { std::string id,provider; ResidualOperation operation; ResidualAvailability available; };
    std::vector<ResidualEntry> residuals_;
};

class PlanExecutor {
public:
    /// @brief Fail before execution if any Plan leaf lacks a provider.
    static void validateBindings(const System::CompiledSolvePlan& plan,
                                 const OpRegistry& operations);

    /// @brief 启动时验证并封存，调用者必须保持该计划不变且存活至操作表释放。
    static void sealBindings(const System::CompiledSolvePlan& plan,OpRegistry& operations);

    static void execute(const System::CompiledSolvePlan& plan,
                        const OpRegistry& operations,
                        const PlanTraceContext* trace = nullptr);

};

} // namespace SF::Run
