/// @file SF_parallelContext.cpp
/// @brief halo、GlobalDof 汇总与并行协调基础设施实现。

#include "SF_parallelContext.h"

#include "SF_MultiBlockMesh.h"

namespace SF::Parallel {

bool mpiCompiled() {
#if defined(SF_USE_MPI) && SF_USE_MPI
    return true;
#else
    // 没有 MPI 时 SF_mpi 仍会构建 serial stub；stub 不是并行能力。
    return false;
#endif
}

namespace {

std::unique_ptr<Backend::CommunicationBackend> makeBackend(
        int& argc, char**& argv, bool requested) {
    if (requested) {
        return std::make_unique<Backend::MPIBackend>(argc, argv, true);
    }
    return std::make_unique<Backend::SerialBackend>();
}

} // namespace

ParallelContext::ParallelContext(
        int& argc, char**& argv, bool requested)
    : backend_(makeBackend(argc, argv, requested)),
      processorFlux_(&halo_),
      reduction_(backend_.get()),
      decomposition_(backend_.get()),
      coordinator_(&halo_, nullptr, backend_.get(),
                   &reduction_, &decomposition_, &processorFlux_, false) {}

void ParallelContext::configureSingle(
        const MeshCommunication::HaloExchangePlan* plan,
        int localBlockId,
        std::vector<MeshBlockField>* blocks,
        bool canonicalInterfaceFlux) {
    topology_.configure(plan, localBlockId, rank(), size());
    halo_.configure(topology_.plan(), topology_.localBlockId(), backend_.get());
    // 单 Field 求解仍可访问全局 block 拓扑。Runtime 仅在 processor face
    // 处通过这份 storage 找到本 rank 所有的 canonical flux owner；算法
    // 继续只推进一个 Field，不需要了解 MPI 分块。
    coordinator_.configureBlocks(blocks, canonicalInterfaceFlux);
}

void ParallelContext::configureBlocks(
        const MeshCommunication::HaloExchangePlan* plan,
        std::vector<MeshBlockField>* blocks,
        bool canonicalInterfaceFlux) {
    topology_.configure(plan, rank(), rank(), size());
    halo_.configure(topology_.plan(), topology_.localBlockId(), backend_.get());
    coordinator_.configureBlocks(blocks, canonicalInterfaceFlux);
}

void ParallelContext::exchangeInitial(Field& field) {
    if (halo_.active()) halo_.exchange(field);
}

void ParallelContext::exchangeInitial(
        std::vector<MeshBlockField>& blocks) {
    if (halo_.active()) halo_.exchange(blocks);
}

} // namespace SF::Parallel
