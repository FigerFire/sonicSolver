/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_interfaces.h
/// @brief interfaces 目录总入口及求解器服务接口。
///
/// Solver Algorithm owns numerical step orchestration, but it does not own MPI,
/// halo exchange, IBM implementation details, output, or logging. These
/// interfaces let those concerns be injected as independent modules. Observer、
/// transport 和 equation coupling 已按职责拆到各自头文件；本文件保留同步、IBM、
/// stepper 和服务表，并作为兼容总入口。

#include "core/field/SF_field.h"
#include "SF_equationCoupling.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "SF_observer.h"
#include "SF_transportModel.h"
#include "SF_immersedSystem.h"
#include "SF_immersedConstraint.h"
#include "SF_solveStrategy.h"
#include "core/state/SF_state.h"

#include <functional>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace SF {
namespace Boundary { class Applicator; }
namespace FDM {

/// @brief 已装配边界管线的统一 stage 入口。
///
/// 管线实现拥有 fitted/IBM geometry、物理边界律和 algebraic/ILW 重构的组合；
/// flow algorithm 只请求“准备可供离散读取的边界状态”。
class IBoundaryPipeline {
public:
    virtual ~IBoundaryPipeline() = default;
    virtual void prepare(const std::vector<Field*>& fields,
                         double time, double dt) = 0;
};

/// @brief Interface for immersed-boundary ghost-cell enforcement.
///
/// IBM setup remains owned by the IBM module. The solver only calls this hook
/// when a step or RK sub-step needs fresh ghost-cell values.
class IImmersedBoundary {
public:
    virtual ~IImmersedBoundary() = default;

    /// @brief Apply immersed-boundary constraints to the field.
    /// @param field Field modified in-place.
    virtual void apply(Field& field, double time, double dt) = 0;

    /// @brief 多 patch 入口；默认逐 patch 应用。
    virtual void apply(const std::vector<Field*>& fields,
                       double time, double dt) {
        for (Field* field : fields) {
            if (field) apply(*field, time, dt);
        }
    }
};

/// @brief IBM 对求解器公开的组合端口。
///
/// 几何 ghost 重建、约束代数和方法描述属于同一个 IBM 模块，但进入求解流程的
/// 阶段不同。用一个值对象传递三类端口，避免 algorithm/application 保存三组可能
/// 相互不一致的裸指针；未启用的端口保持为空。
struct ImmersedCouplingPorts {
    IImmersedBoundary* boundary = nullptr;
    IImmersedConstraint* constraint = nullptr;
    IImmersedSystem* system = nullptr;

    bool hasBoundaryReconstruction() const { return boundary != nullptr; }
    bool hasVariationalConstraint() const {
        return constraint != nullptr || system != nullptr;
    }
};

/// @brief Solver 可见的唯一并行阶段接口。
///
/// MPI communicator、patch graph、rank 和消息缓冲均由实现类拥有。Solver 只在
/// 明确的 stage 边界请求状态同步、canonical 接口通量装配和全局归约。
class IParallelCoordinator {
public:
    virtual ~IParallelCoordinator() = default;

    /// @brief 绑定本次推进的唯一状态包；实现不接管其所有权。
    virtual void attachState(State::StateBundle& state) = 0;

    /// @brief 按注册契约同步当前阶段需要的全部状态。
    virtual void synchronizeRegistered(State::HaloSyncStage stage) = 0;

    /// @brief 保证指定注册字段至少具有 requiredDepth 层新鲜 halo。
    virtual void ensureHalo(const std::string& name,
                            int requiredDepth) = 0;

    /// @brief 声明注册字段刚被本地算法写入，使既有 halo 失效。
    virtual void markModified(const std::string& name) = 0;

    /// @brief 同步未长期注册的算法 workspace；用于压力矩阵编号等短生命周期状态。
    virtual void synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) = 0;

    /// @brief 广播唯一 F*，按 face incidence 对 owner/neighbour 装配相反残差。
    virtual void assembleCanonicalInterfaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) = 0;

    /// @brief 将各 patch 的积分贡献汇总到 GlobalDof owner，再广播唯一结果。
    virtual void accumulateGlobalDof(const std::string& field,
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) = 0;

    enum class Reduction { Minimum, Maximum, Sum };

    /// @brief 执行与字段类型无关的全局标量归约。
    virtual double reduce(double localValue, Reduction operation) = 0;

    /// @brief 对显式临时向量作逐项全局 SUM；不得从 Field 访问器隐式调用。
    virtual void reduceSum(std::vector<double>& values) = 0;

    /// @brief 将 canonical owner 的稀疏实体值复制到所有 replicas。
    /// @param ownerMask 与 entityIds 等长；仅 owner 项为 1，其余项必须为 0。
    virtual void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) = 0;

    /// @brief 以各执行单元的本地数量分配连续、无重叠的全局编号区间。
    virtual DistributedIndexRange distributeIndices(
        std::int64_t localCount) = 0;

    /// @brief 从 root 向全部执行单元广播一段通用数值数据。
    virtual void broadcast(std::vector<double>& values, int root) = 0;

    /// @brief 所有 rank 对布尔完成状态取逻辑与。
    virtual bool allRanksAgree(bool localValue) = 0;

    virtual bool active() const = 0;
    virtual int rank() const = 0;
    virtual int size() const = 0;
    virtual void barrier() = 0;
};

/// @brief 时间驱动器传给 stepper 的非拥有 state 请求。
///
/// 物理 state、patch、clock 和 equation binding 全部在 `bundle` 中；此对象只
/// 携带 bundle 与本步硬上限。
struct SolverState {
    /// @brief 唯一参与推进的物理 state bundle。
    State::StateBundle* bundle = nullptr;
    /// @brief 外层时间/输出驱动器给出的本步硬上限。
    double maximumTimeStep = std::numeric_limits<double>::max();
};

/// @brief Solver workflow 可调用的外围服务表。
///
/// 这些指针均可为空；为空时 stepper 使用自身构造时保存的默认配置或 no-op 行为。
/// 目标是让 MPI、IBM、湍流、边界应用和观察者通过小接口注入，而不是让 solver
/// 继续读取越来越多的全局 parser 状态。
struct SolverServices {
    /// @brief 统一边界状态管线；提供时优先于分散的边界/halo/IBM 服务。
    IBoundaryPipeline* boundaryPipeline = nullptr;
    /// @brief 可选边界应用器；为空时使用 stepper 内部配置生成的应用器。
    Boundary::Applicator* boundary = nullptr;
    /// @brief IBM 边界重建、变分约束和方法描述的组合端口。
    ImmersedCouplingPorts immersed;
    /// @brief 显式 contract/freshness 驱动的确定性执行 Runtime。
    IExecutionRuntime* executionRuntime = nullptr;
    /// @brief 湍流/有效传输模型服务。
    ITransportModel* transportModel = nullptr;
    /// @brief 可选附加方程系统。
    IEquationSystemCoupling* equationSystem = nullptr;
    /// @brief 结构化日志、GUI 或测试观察者。
    ISolverObserver* observer = nullptr;
};

/// @brief 一个 NS stepper 完成一次推进后的结构化结果。
struct StepResult {
    /// @brief 本次推进是否成功完成。
    bool accepted = false;
    /// @brief 实际采用的 dt。
    double dt = 0.0;
    /// @brief 推进后的物理时间。
    double time = 0.0;
    /// @brief 推进后的步编号。
    int step = 0;
    /// @brief 是否建议外层输出；当前过渡实现由外层输出调度决定。
    bool outputDue = false;
    /// @brief 是否已经到达终止条件；当前过渡实现由外层 stop 条件决定。
    bool finished = false;
    /// @brief 失败或诊断信息；成功时可为空。
    std::string message;
};

/// @brief Resolved workflow 中一个有序、可执行的 solve block。
///
/// 该值对象只描述 equation/constraint identity，不复制物理模型或离散参数。
/// Workflow builder 与 stepper 共用同一表示，避免 explain plan 与 runtime plan
/// 再次形成两份 authority。
enum class SolveStageKind {
    Predictor,
    Correction,
    Constraint,
    Commit
};

struct SolveStage {
    std::string id;
    SolveStageKind kind = SolveStageKind::Commit;
    SolveStrategyKind strategy = SolveStrategyKind::AlgebraicUpdate;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
};

/// @brief Navier-Stokes 方程推进唯一接口。
///
/// 密度基、压力基、单块、多块、未来 GPU backend 都应以这个接口作为外部入口。
class INavierStokesStepper {
public:
    virtual ~INavierStokesStepper() = default;

    /// @brief 在时间循环开始前绑定全部稳定服务。
    virtual void bindServices(SolverServices services) = 0;

    /// @brief 在时间循环开始前绑定已解析的有序 solve blocks。
    ///
    /// Specialized stepper 可以把一个 block lowering 为内部嵌套算法，但必须
    /// fail fast 拒绝缺失或无法执行的 block；advance 期间不得替换此计划。
    virtual void bindSolveStages(const std::vector<SolveStage>& stages) = 0;

    /// @brief 在进入时间循环前绑定并验证 resolved runtime state。
    virtual void prepare(SolverState& state) = 0;

    /// @brief 推进一个时间步或一个 workflow iteration。
    /// @param state 可变求解状态包。
    /// @return 结构化推进结果。
    virtual StepResult advance(SolverState& state) = 0;
};

/// @brief No-op immersed-boundary hook for non-IBM runs.
class NoopImmersedBoundary final : public IImmersedBoundary {
public:
    /// @brief Leave the field unchanged.
    /// @param field Unused field reference.
    void apply(Field&, double, double) override {}
};

/// @brief 无物理边界、IBM 或 halo 的空边界管线。
class NoopBoundaryPipeline final : public IBoundaryPipeline {
public:
    void prepare(const std::vector<Field*>&, double, double) override {}
};

/// @brief 串行或单 patch 工作流的空并行协调器。
class LocalParallelCoordinator final : public IParallelCoordinator {
public:
    void attachState(State::StateBundle& state) override { state_ = &state; }
    void synchronizeRegistered(State::HaloSyncStage) override {}
    void ensureHalo(const std::string&, int) override {}
    void markModified(const std::string&) override {}
    void synchronizeTransient(
        const std::vector<State::DistributedFieldView>&) override {}
    void assembleCanonicalInterfaceFluxes(const std::vector<Field*>&,
                                          const std::vector<FluxField*>&,
                                          const std::vector<Residual*>&) override {}
    void accumulateGlobalDof(const std::string&, const std::vector<Field*>&,
                             const std::vector<Residual*>&) override {}
    double reduce(double localValue, Reduction) override { return localValue; }
    void reduceSum(std::vector<double>&) override {}
    void copyCanonicalEntities(
        const std::vector<std::int64_t>&,
        std::vector<double>&,
        const std::vector<unsigned char>&) override {}
    DistributedIndexRange distributeIndices(
            std::int64_t localCount) override {
        return {0, localCount - 1, localCount};
    }
    void broadcast(std::vector<double>&, int) override {}
    bool allRanksAgree(bool localValue) override { return localValue; }
    bool active() const override { return false; }
    int rank() const override { return 0; }
    int size() const override { return 1; }
    void barrier() override {}
private:
    State::StateBundle* state_ = nullptr;
};

/// @brief No-op transport model for laminar/test execution.
class NoopTransportModel final : public ITransportModel {
public:
    /// @brief 不做任何湍流边界处理。
    /// @param field 未使用。
    void applyBoundary(const Field&) override {}

    /// @brief 不更新任何湍流状态。
    /// @param field 未使用。
    /// @param dt 未使用。
    void correct(const Field&, double) override {}

    /// @brief 返回基础层流粘度。
    /// @param field 未使用。
    /// @param i 未使用。
    /// @param j 未使用。
    /// @param k 未使用。
    /// @param laminarMu 基础层流粘度。
    /// @return `laminarMu`。
    double dynamicViscosity(const Field&, int, int, int, double laminarMu) const override {
        return laminarMu;
    }
};

} // namespace FDM
} // namespace SF
