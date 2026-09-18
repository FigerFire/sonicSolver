/// @file SF_ibmAdapters.cpp
/// @brief 转发 solver immersed ports 到现有 IB/CompositeIB 数值实现。
///
/// Data flow:
///   solver boundary/constraint request
///       -> typed IBM adapter
///       -> existing ghost、projection 或 KKT helper
///
/// 本文件只做 interface adaptation，不创建 MPI communicator 或执行循环。

#include "app/application/adapters/SF_ibmAdapters.h"

#include "SF_IBM.h"
#include "SF_MultiBlockMesh.h"
#include "SF_compositeIBM.h"

#include <stdexcept>

namespace SF::Application::Adapters {

IBBoundaryAdapter::IBBoundaryAdapter(
        IBM::IB* ibm)
    : ibm_(ibm) {}

void IBBoundaryAdapter::apply(
        Field& field, double time, double dt) {
    if (ibm_) ibm_->applyGhostCells(field, time, dt);
}

IBConstraintAdapter::IBConstraintAdapter(
        IBM::IB* ibm)
    : ibm_(ibm) {}

const FDM::IImmersedSystem*
IBConstraintAdapter::systemProvider() const {
    return ibm_;
}

void IBConstraintAdapter::setExecutionRuntime(
        FDM::IExecutionRuntime* runtime) {
    if (ibm_) ibm_->setExecutionRuntime(runtime);
}

FDM::ImmersedConstraintResult
IBConstraintAdapter::projectPredictedState(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    return ibm_ ? ibm_->applyConstraint(fields, targetTime, dt)
                : FDM::ImmersedConstraintResult{};
}

const FDM::ImmersedSurfaceSystem&
IBConstraintAdapter::prepareMonolithicSystem(
        Field& field, double targetTime, double dt) {
    if (!ibm_) {
        throw std::runtime_error("IBM constraint adapter is not configured.");
    }
    return ibm_->prepareMonolithicSystem(field,targetTime,dt);
}

FDM::ImmersedConstraintResult
IBConstraintAdapter::acceptMonolithicSolution(
        Field& field, double targetTime, double dt,
        const FDM::ImmersedKKTState& state) {
    if (!ibm_) {
        throw std::runtime_error("IBM constraint adapter is not configured.");
    }
    return ibm_->acceptMonolithicSolution(field,targetTime,dt,state);
}

MultiPatchIBAdapter::MultiPatchIBAdapter(
        IBM::CompositeIB* ibm,
        MultiBlockMesh* mesh,
        const std::vector<int>* localPatchIds)
    : ibm_(ibm), mesh_(mesh), localPatchIds_(localPatchIds) {}

void MultiPatchIBAdapter::apply(Field&, double, double) {
    throw std::runtime_error(
        "MultiPatchIBAdapter requires the multi-patch entry.");
}

void MultiPatchIBAdapter::apply(
        const std::vector<Field*>&, double time, double dt) {
    if (ibm_ && mesh_ && localPatchIds_) {
        ibm_->applyLocalPatches(*mesh_, *localPatchIds_, time, dt);
    }
}

} // namespace SF::Application::Adapters
