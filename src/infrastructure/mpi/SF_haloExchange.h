/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.08.21-----------*/

#pragma once

/// @file SF_haloExchange.h
/// @brief 基于 Mesh communication plan 的 Field halo 与 processor-face 通量通信。

#include "SF_communicationPlan.h"
#include "core/field/SF_field.h"
#include "core/flux/SF_flux.h"
#include "core/residual/SF_residual.h"
#include "backend/SF_communicationBackend.h"

#include <string>
#include <vector>

namespace SF {
struct MeshBlockField;
namespace Parallel {

/// @brief 一个与结构块 Field 同尺寸的标量数组视图。
struct ScalarBlockValues {
    int blockId = -1;
    const Field* field = nullptr;
    std::vector<double>* values = nullptr;
    std::string name;
};

/// @brief MPI 无关的局部通信能力；实际传输由 backend 注入。
class HaloExchange {
public:
    void configure(
        const MeshCommunication::HaloExchangePlan* plan,
        int localBlockId,
        const Backend::CommunicationBackend* backend);

    bool active() const {
        return backend_ && backend_->active()
            && plan_ != nullptr && localBlockId_ >= 0;
    }
    int localBlockId() const { return localBlockId_; }

    void exchange(Field& field) const;
    void exchange(std::vector<MeshBlockField>& blocks) const;
    void exchangeScalarValues(
        const Field& field, std::vector<double>& values) const;
    void exchangeScalarValues(
        const std::vector<ScalarBlockValues>& views) const;
    void exchangeScalarIdentifiers(
        const std::vector<ScalarBlockValues>& views) const;
    void synchronizeCanonicalFaceFlux(
        const Field& field, std::vector<double>& values) const;
    void assembleCanonicalInterfaceFluxes(
        std::vector<MeshBlockField>& blocks,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) const;
    /// @brief 汇总各 patch 的 GlobalDof 积分残差并广播唯一结果。
    void assembleGlobalDofResiduals(
        std::vector<MeshBlockField>& blocks,
        const std::vector<Residual*>& residuals) const;
    void assembleGlobalDofResiduals(Field& localField, Residual& residual) const;

private:
    const MeshCommunication::HaloExchangePlan* plan_ = nullptr;
    const Backend::CommunicationBackend* backend_ = nullptr;
    int localBlockId_ = -1;
    int rank_ = 0;
    int size_ = 1;
    // 以下三个值是通信计划算法的本地执行上下文，不是 MPI 对象。
    int mpiRank_ = 0;
    int mpiSize_ = 1;
    bool mpiEnabled_ = false;

    int ownerRank(int blockId) const;
    const MeshCommunication::HaloBlockPlan* localPlan() const;
    const MeshCommunication::HaloBlockPlan* blockPlan(int blockId) const;
    std::vector<int> communicationRanks() const;
    std::vector<int> rankInteriorCounts(int components) const;
    void exchangeNeighbourPayload(
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& displacements,
        int tag,
        std::vector<double>& received) const;
    void exchangeScalarBlockValues(
        const std::vector<ScalarBlockValues>& views,
        bool synchronizeSharedPoints) const;
};

} // namespace Parallel
} // namespace SF
