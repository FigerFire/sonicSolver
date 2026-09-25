#pragma once

/// @file SF_parallelContext.h
/// @brief Application 装配并行能力的唯一入口。

#include "SF_domainDecomposition.h"
#include "SF_haloExchange.h"
#include "SF_parallelCoordinator.h"
#include "SF_parallelTopology.h"
#include "SF_processorFlux.h"
#include "SF_reduction.h"
#include "backend/SF_mpiBackend.h"
#include "backend/SF_serialBackend.h"

#include <memory>
#include <vector>

namespace SF {
struct MeshBlockField;
namespace Parallel {

/// @brief 本 binary 是否链接了真实 MPI 通信后端。
///
/// 这是构建事实，不是 case 属性：`[parallel].enabled` 只是 case 的请求，
/// capability 报告必须回答“这个可执行文件能不能真的通信”。
bool mpiCompiled();

class ParallelContext {
public:
    ParallelContext(int& argc, char**& argv, bool requested);

    bool available() const { return backend_->available(); }
    bool active() const { return backend_->active(); }
    bool isRoot() const { return rank() == 0; }
    int rank() const { return backend_->rank(); }
    int size() const { return backend_->size(); }
    bool allRanksAgree(bool localValue) {
        return reduction_.allRanksAgree(localValue);
    }

    FDM::IParallelCoordinator& coordinator() { return coordinator_; }

    void configureSingle(
        const MeshCommunication::HaloExchangePlan* plan,
        int localBlockId,
        std::vector<MeshBlockField>* blocks = nullptr,
        bool canonicalInterfaceFlux = false);
    void configureBlocks(
        const MeshCommunication::HaloExchangePlan* plan,
        std::vector<MeshBlockField>* blocks,
        bool canonicalInterfaceFlux);
    void exchangeInitial(Field& field);
    void exchangeInitial(std::vector<MeshBlockField>& blocks);

private:
    std::unique_ptr<Backend::CommunicationBackend> backend_;
    ParallelTopology topology_;
    HaloExchange halo_;
    ProcessorFlux processorFlux_;
    ReductionService reduction_;
    DomainDecomposition decomposition_;
    ParallelCoordinator coordinator_;
};

} // namespace Parallel
} // namespace SF
