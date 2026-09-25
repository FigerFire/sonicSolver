#pragma once

/// @file SF_method.h
/// @brief 统一 IBM 方法调度 facade 与约束生命周期适配。

#include "SF_ibmConfig.h"
#include "SF_immersedConstraint.h"
#include "SF_immersedSystem.h"
#include "transfer/body/SF_bodyOperator.h"
#include "solid/SF_bodyModel.h"
#include "solid/rigid/SF_rigidBody.h"
#include "transfer/surface/SF_surfaceOperator.h"
#include "geoProcessing/SF_STLGeometry.h"

#include <memory>
#include <utility>
#include <vector>

namespace SF::IBM::Forcing {

/// @brief 扩展流体域上的离散变分约束 facade。
///
/// 当前实现求解 prescribed rigid body 的体约束。对每个受约束自由度，离散驻值
/// 条件为 `M(u-u*)-dt*C^T*lambda=0, C*u=U_s`；乘子由该系统自动恢复，
/// 而不是由经验 forcing 公式预先指定。表面乘子与 coupled-solid KKT 使用同一
/// `IImmersedConstraint` 生命周期，但尚未实现时会明确拒绝。
class ImmersedForcingSystem final : public FDM::IImmersedConstraint {
public:
    ImmersedForcingSystem();

    /// @brief 注入 body constraint operator；空指针直接拒绝。
    void setBodyConstraintOperator(
        std::unique_ptr<IBodyConstraintOperator> operatorInstance);
    /// @brief 注入 surface constraint operator；空指针直接拒绝。
    void setSurfaceConstraintOperator(
        std::unique_ptr<ISurfaceConstraintOperator> operatorInstance);
    /// @brief 注入固体运动/动力学模型；空指针直接拒绝。
    void setBodyModel(std::unique_ptr<IBodyModel> bodyModel);

    /// @brief 绑定参考构型 STL 与只读 forcing 配置。
    void configure(
        const GeoProcessing::STLGeometry& geometry,
        IBMRuntimeConfig config);

    FDM::ImmersedConstraintResult projectPredictedState(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) override;
    void setExecutionRuntime(FDM::IExecutionRuntime* runtime) override;

    const FDM::ImmersedSurfaceSystem& prepareMonolithicSystem(
        Field& field, double targetTime, double dt) override;
    FDM::ImmersedConstraintResult acceptMonolithicSolution(
        Field& field, double targetTime, double dt,
        const FDM::ImmersedKKTState& state) override;

    /// @brief 查询最近一次校正的体约束 mask。
    double constraintMask(const Field& field, int i, int j, int k) const;
    /// @brief 查询最近一次校正的乘子/流体力密度，单位 N/m3。
    Vector3 multiplier(const Field& field, int i, int j, int k) const;
    const FDM::ImmersedConstraintResult& lastResult() const {
        return lastResult_;
    }
    /// @brief 返回配置阶段冻结的 IBM 数学描述。
    const FDM::ImmersedAlgorithmDescriptor& algorithmDescriptor() const {
        return descriptor_;
    }

private:
    FDM::ImmersedConstraintResult applyPeskinOriginal(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    FDM::ImmersedConstraintResult applyDFMExplicitSelfPropelled(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    FDM::ImmersedConstraintResult applyDFMFractionalStepSelfPropelled(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    FDM::ImmersedConstraintResult applyDFMFractionalStepPrescribed(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    FDM::ImmersedConstraintResult applyVelocityForcingFTS(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    FDM::ImmersedConstraintResult applyVelocityForcingBP(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt);
    void validateImplicitAlgorithmContract() const;
    bool usesMonolithicKKT() const;
    /// @brief 用 Runtime 的显式归并完成分区一致的 marker 归一化。
    void buildSurfaceSystem(Field& field, double targetTime, double dt);
    /// @brief 对 partial J 结果执行 Runtime 的逐 marker 向量 SUM。
    void reduceSurfaceVectors(std::vector<Vector3>& values) const;
    /// @brief 将 constraint owner 写出的表面 `lambda` COPY 给 J^T consumer。
    ///
    /// 每个 marker 分量经过 Runtime 的 canonical-copy 边界。此处不执行 SUM
    /// 或 AVG；consumer 只读取其唯一 GlobalConstraintDof owner 的值。
    void copySurfaceMultipliers(std::vector<Vector3>& values) const;
    /// @brief 一个 marker 的全局唯一诊断/载荷 owner，不暴露 rank 给 IBM 算法。
    bool ownsSurfacePoint(const FDM::ImmersedSurfacePoint& point) const;
    /// @brief 将 owner-local IBM 载荷、泛函和残差归并成全局诊断。
    void finalizeDistributedResult(FDM::ImmersedConstraintResult& result) const;

    const GeoProcessing::STLGeometry* geometry_ = nullptr;
    const Field* lastField_ = nullptr;
    IBMRuntimeConfig config_;
    std::vector<unsigned char> mask_;
    std::vector<Vector3> multiplier_;
    FDM::ImmersedConstraintResult lastResult_;
    /// @brief Peskin/显式 DFM 下一时间步使用的滞后 Eulerian 乘子。
    std::vector<Vector3> laggedMultiplier_;
    std::unique_ptr<ISurfaceConstraintOperator> surfaceOperator_;
    std::unique_ptr<IBodyConstraintOperator> bodyOperator_;
    std::unique_ptr<IBodyModel> bodyModel_;
    FDM::ImmersedSurfaceSystem surfaceSystem_;
    /// @brief 最近一次建图时的全局 raw `J` 权重和，仅用于分区一致性诊断。
    std::vector<double> surfaceNormalizations_;
    FDM::ImmersedAlgorithmDescriptor descriptor_;
    FDM::IExecutionRuntime* runtime_ = nullptr;

};

} // namespace SF::IBM::Forcing
