#pragma once

/// @file SF_ibmAdapters.h
/// @brief 将 concrete IBM implementation 暴露为 solver immersed interfaces。
///
/// Data flow:
///   IBM::IB / CompositeIB
///       -> application adapter
///       -> IImmersedBoundary / IImmersedConstraint
///       -> existing numerical caller
///
/// Adapter 不选择 solver runner、IBM 数学方法或 timestep lifecycle。

#include "SF_interfaces.h"

#include <vector>

namespace SF {
class MultiBlockMesh;
struct MeshBlockField;
namespace IBM {
class IB;
class CompositeIB;
}
namespace Application::Adapters {

class IBBoundaryAdapter final : public FDM::IImmersedBoundary {
public:
    explicit IBBoundaryAdapter(IBM::IB* ibm);
    void apply(Field& field, double time, double dt) override;
private:
    IBM::IB* ibm_ = nullptr;
};

/// @brief 把 IBM manager 的变分乘子阶段注入统一 flow algorithm。
class IBConstraintAdapter final : public FDM::IImmersedConstraint {
public:
    explicit IBConstraintAdapter(IBM::IB* ibm);
    const FDM::IImmersedSystem* systemProvider() const override;
    void setExecutionRuntime(FDM::IExecutionRuntime* runtime) override;
    FDM::ImmersedConstraintResult projectPredictedState(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) override;
    const FDM::ImmersedSurfaceSystem& prepareMonolithicSystem(
        Field& field, double targetTime, double dt) override;
    FDM::ImmersedConstraintResult acceptMonolithicSolution(
        Field& field, double targetTime, double dt,
        const FDM::ImmersedKKTState& state) override;
private:
    IBM::IB* ibm_ = nullptr;
};

class MultiPatchIBAdapter final
    : public FDM::IImmersedBoundary {
public:
    MultiPatchIBAdapter(
        IBM::CompositeIB* ibm,
        MultiBlockMesh* mesh,
        const std::vector<int>* localPatchIds);
    void apply(Field&, double time, double dt) override;
    void apply(const std::vector<Field*>& fields,
               double time, double dt) override;
private:
    IBM::CompositeIB* ibm_ = nullptr;
    MultiBlockMesh* mesh_ = nullptr;
    const std::vector<int>* localPatchIds_ = nullptr;
};

} // namespace Application::Adapters
} // namespace SF
