/// @file SF_parallelCoordinator.cpp
/// @brief halo、GlobalDof 汇总与并行协调基础设施实现。

#include "SF_parallelCoordinator.h"

#include "SF_MultiBlockMesh.h"

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::Parallel {
namespace {

bool highOrderTraceEnabled() {
    static const bool enabled = [] {
        const char* env = std::getenv("SF_HIGH_ORDER_TRACE");
        return env && std::string(env) == "1";
    }();
    return enabled;
}

bool highOrderTraceActiveFor(std::int64_t step) {
    if (!highOrderTraceEnabled()) return false;
    static const std::int64_t selectedStep = [] {
        const char* env = std::getenv("SF_HIGH_ORDER_TRACE_STEP");
        if (!env || *env == '\0') return std::int64_t{0};
        char* end = nullptr;
        errno = 0;
        const long long parsed = std::strtoll(env, &end, 10);
        if (errno != 0 || end == env || *end != '\0' || parsed < 0) {
            throw std::runtime_error(
                "SF_HIGH_ORDER_TRACE_STEP must be a non-negative integer.");
        }
        return static_cast<std::int64_t>(parsed);
    }();
    return step == selectedStep;
}

void traceWorkspaceMapping(
        const std::vector<MeshBlockField>& blocks,
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals,
        int rank,
        std::int64_t step,
        const char* barrier) {
    if (!highOrderTraceActiveFor(step)) return;
    for (std::size_t block = 0; block < blocks.size(); ++block) {
        std::size_t matchCount = 0;
        std::size_t patchIndex = 0;
        for (std::size_t patch = 0; patch < fields.size(); ++patch) {
            if (fields[patch] == &blocks[block].field) {
                ++matchCount;
                patchIndex = patch;
            }
        }
        std::cout << "[SF TRACE] step=" << step
                  << " stage=runtime-barrier " << barrier
                  << " rank=" << rank
                  << " block=" << block
                  << " owner-rank=" << blocks[block].ownerRank
                  << " matching-workspaces=" << matchCount;
        if (matchCount == 1) {
            std::cout << " patch=" << patchIndex
                      << " field=" << static_cast<const void*>(fields[patchIndex])
                      << " flux=" << static_cast<const void*>(fluxes[patchIndex])
                      << " residual=" << static_cast<const void*>(residuals[patchIndex]);
        }
        std::cout << '\n';
        if (matchCount > 1 || (blocks[block].ownerRank == rank && matchCount != 1)) {
            throw std::runtime_error(
                "high-order trace found a non-bijective local block/workspace mapping.");
        }
    }
}

} // namespace

ParallelCoordinator::ParallelCoordinator(
        HaloExchange* halo,
        std::vector<MeshBlockField>* blocks,
        const Backend::CommunicationBackend* backend,
        const ReductionService* reduction,
        const DomainDecomposition* decomposition,
        const ProcessorFlux* processorFlux,
        bool canonicalInterfaceFlux)
    : halo_(halo), blocks_(blocks), backend_(backend),
      reduction_(reduction), decomposition_(decomposition),
      processorFlux_(processorFlux),
      canonicalInterfaceFlux_(canonicalInterfaceFlux) {}

void ParallelCoordinator::attachState(State::StateBundle& state) {
    state_ = &state;
}

void ParallelCoordinator::synchronizeRegistered(
        State::HaloSyncStage stage) {
    if (!state_) {
        throw std::runtime_error(
            "ParallelCoordinator has no attached StateBundle.");
    }
    for (const auto& name : state_->distributed.names(stage)) {
        auto fields = state_->distributed.select(name, stage);
        exchangeViews(fields);
        int depth = 0;
        for (const auto* field : fields) {
            if (field && field->haloDepth > depth) depth = field->haloDepth;
        }
        state_->distributed.markSynchronized(name, depth);
    }
}

void ParallelCoordinator::ensureHalo(
        const std::string& name, int requiredDepth) {
    if (!state_) {
        throw std::runtime_error(
            "ParallelCoordinator::ensureHalo has no attached StateBundle.");
    }
    // Halo freshness 是 owner-local 状态，但一次 MPI halo exchange 是所有
    // 拓扑邻居必须共同参与的 collective phase。若仅有一个 rank 的 owned
    // state 被写入，它会请求交换；其余 rank 即使本地 cache 尚新鲜也不能
    // 跳过，否则发送方会在 MPI_Waitall 中等待一个不存在的配对通信。
    const bool localExchange = state_->distributed.needsExchange(
        name, requiredDepth);
    const bool globalExchange = reduction_
        ? reduction_->maximum(localExchange ? 1.0 : 0.0) > 0.5
        : localExchange;
    if (!globalExchange) return;
    auto fields = state_->distributed.select(
        name, State::HaloSyncStage::None, requiredDepth);
    exchangeViews(fields);
    state_->distributed.markSynchronized(name, requiredDepth);
}

void ParallelCoordinator::markModified(const std::string& name) {
    if (!state_) {
        throw std::runtime_error(
            "ParallelCoordinator::markModified has no attached StateBundle.");
    }
    state_->distributed.markModified(name);
}

void ParallelCoordinator::synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) {
    std::vector<State::DistributedFieldView*> views;
    views.reserve(fields.size());
    for (const auto& field : fields) {
        // 同步过程只改写视图引用的数据，不修改视图描述本身。
        views.push_back(const_cast<State::DistributedFieldView*>(&field));
    }
    exchangeViews(views);
}

void ParallelCoordinator::exchangeViews(
        const std::vector<State::DistributedFieldView*>& fields) {
    if (fields.empty() || !halo_ || !halo_->active()) return;

    // Pressure corrector 等算法只拥有本 rank 的 patch workspace，而 halo plan
    // 仍以全局 block 列表构造固定消息段。由协调器补齐非本地 scratch，保证所有
    // rank 的 payload 形状一致；物理/方程模块无需认识全局拓扑。
    if (blocks_ && fields.size() != blocks_->size()) {
        std::vector<std::vector<double>> scratch(blocks_->size());
        std::vector<State::DistributedFieldView> expandedViews;
        expandedViews.reserve(blocks_->size());
        std::vector<State::DistributedFieldView*> expanded;
        expanded.reserve(blocks_->size());
        const auto& prototype = *fields.front();
        for (size_t block = 0; block < blocks_->size(); ++block) {
            Field& geometry = (*blocks_)[block].field;
            State::DistributedFieldView* matched = nullptr;
            for (auto* field : fields) {
                if (field->geometry == &geometry) {
                    matched = field;
                    break;
                }
            }
            if (matched) {
                expandedViews.push_back(*matched);
                expandedViews.back().blockId = (int)block;
            } else {
                scratch[block].assign(
                    (size_t)geometry.TotalSize()
                        * (size_t)prototype.components,
                    0.0);
                expandedViews.push_back(State::workspaceView(
                    prototype.name, (int)block, geometry,
                    scratch[block], prototype.components,
                    prototype.haloDepth, prototype.stages,
                    prototype.exchange, prototype.location));
            }
        }
        for (auto& field : expandedViews) expanded.push_back(&field);
        exchangeViews(expanded);
        return;
    }
    const auto kind = fields.front()->exchange;
    const int components = fields.front()->components;
    for (const auto* field : fields) {
        field->validate();
        if (field->exchange != kind || field->components != components) {
            throw std::runtime_error(
                "A distributed field group has inconsistent exchange metadata.");
        }
    }
    if (kind == State::ExchangeKind::None) return;

    if (kind == State::ExchangeKind::CanonicalFaceFlux) {
        if (fields.size() != 1 || components != 3) {
            throw std::runtime_error(
                "CanonicalFaceFlux transient exchange requires one three-component view.");
        }
        auto* field = fields.front();
        std::vector<double> packed(
            (size_t)field->geometry->TotalSize() * 3u);
        for (int component = 0; component < 3; ++component) {
            for (int cell = 0; cell < field->geometry->TotalSize(); ++cell) {
                packed[(size_t)component * (size_t)field->geometry->TotalSize()
                       + (size_t)cell] = field->read(cell, component);
            }
        }
        halo_->synchronizeCanonicalFaceFlux(*field->geometry, packed);
        for (int component = 0; component < 3; ++component) {
            for (int cell = 0; cell < field->geometry->TotalSize(); ++cell) {
                field->write(cell, component,
                    packed[(size_t)component * (size_t)field->geometry->TotalSize()
                           + (size_t)cell]);
            }
        }
        return;
    }

    // 主守恒状态保留原有多分量 Field halo 路径，数值和共享点顺序不变。
    if (fields.front()->name == "conservative") {
        if (blocks_) halo_->exchange(*blocks_);
        else if (fields.size() == 1) halo_->exchange(*fields.front()->geometry);
        else {
            throw std::runtime_error(
                "Conservative state group requires a block container.");
        }
        for (auto* field : fields) field->geometry->invalidateThermodynamicCache();
        return;
    }

    // HaloExchange 当前以标量 block view 为底层传输单元。统一接口在这里
    // 对多分量字段逐分量打包，物理模块不再维护专用 MPI 回调。
    for (int component = 0; component < components; ++component) {
        std::vector<std::vector<double>> buffers(fields.size());
        std::vector<ScalarBlockValues> scalarViews;
        scalarViews.reserve(fields.size());
        for (size_t n = 0; n < fields.size(); ++n) {
            auto* field = fields[n];
            auto& buffer = buffers[n];
            buffer.resize((size_t)field->geometry->TotalSize());
            for (int cell = 0; cell < field->geometry->TotalSize(); ++cell) {
                buffer[(size_t)cell] = field->read(cell, component);
            }
            scalarViews.push_back({
                field->blockId, field->geometry, &buffer, field->name});
        }
        if (fields.size() == 1
            && kind == State::ExchangeKind::State) {
            halo_->exchangeScalarValues(
                *fields.front()->geometry, buffers.front());
        } else if (kind == State::ExchangeKind::Identifier) {
            if (fields.size() == 1) {
                scalarViews.front().blockId = halo_->localBlockId();
            }
            halo_->exchangeScalarIdentifiers(scalarViews);
        } else {
            halo_->exchangeScalarValues(scalarViews);
        }
        for (size_t n = 0; n < fields.size(); ++n) {
            auto* field = fields[n];
            for (int cell = 0; cell < field->geometry->TotalSize(); ++cell) {
                field->write(cell, component, buffers[n][(size_t)cell]);
            }
        }
    }
}

void ParallelCoordinator::assembleCanonicalInterfaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) {
    if (canonicalInterfaceFlux_ && processorFlux_ && blocks_) {
        if (fields.size() != fluxes.size() || fields.size() != residuals.size()) {
            throw std::runtime_error("canonical workspace patch count mismatch.");
        }
        traceWorkspaceMapping(*blocks_, fields, fluxes, residuals, rank(),
                              state_ ? state_->step : -1,
                              "canonical workspace mapping");
        std::vector<FluxField*> byBlock(blocks_->size(), nullptr);
        std::vector<Residual*> residualByBlock(blocks_->size(), nullptr);
        for (size_t block = 0; block < blocks_->size(); ++block) {
            for (size_t patch = 0; patch < fields.size(); ++patch) {
                if (fields[patch] == &(*blocks_)[block].field) {
                    byBlock[block] = fluxes[patch];
                    residualByBlock[block] = residuals[patch];
                    break;
                }
            }
        }
        processorFlux_->assemble(*blocks_, byBlock, residualByBlock);
    }
}

void ParallelCoordinator::accumulateGlobalDof(
        const std::string& field, const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) {
    if (field != "conservativeResidual") {
        throw std::runtime_error(
            "ParallelCoordinator has no GlobalDof accumulator for field '"
            + field + "'.");
    }
    if (!halo_) {
        throw std::runtime_error(
            "GlobalDof residual assembly requires a configured halo service.");
    }
    if (fields.size() != residuals.size()) {
        throw std::runtime_error("GlobalDof workspace patch count mismatch.");
    }
    if (blocks_) {
        std::vector<FluxField*> noFlux(fields.size(), nullptr);
        traceWorkspaceMapping(*blocks_, fields, noFlux, residuals, rank(),
                              state_ ? state_->step : -1,
                              "GlobalDof workspace mapping");
        std::vector<Residual*> byBlock(blocks_->size(), nullptr);
        for (size_t block = 0; block < blocks_->size(); ++block)
            for (size_t patch = 0; patch < fields.size(); ++patch)
                if (fields[patch] == &(*blocks_)[block].field) { byBlock[block] = residuals[patch]; break; }
        halo_->assembleGlobalDofResiduals(*blocks_, byBlock);
    } else if (fields.size() == 1 && residuals.front()) {
        halo_->assembleGlobalDofResiduals(*fields.front(), *residuals.front());
    } else {
        throw std::runtime_error(
            "GlobalDof residual assembly has no local Field storage.");
    }
}

double ParallelCoordinator::reduce(
        double localValue, Reduction operation) {
    if (!reduction_) return localValue;
    switch (operation) {
        case Reduction::Minimum: return reduction_->minimum(localValue);
        case Reduction::Maximum: return reduction_->maximum(localValue);
        case Reduction::Sum: return reduction_->sum(localValue);
    }
    throw std::runtime_error("Unknown parallel reduction operation.");
}

void ParallelCoordinator::reduceSum(std::vector<double>& values) {
    if (reduction_) reduction_->sum(values);
}

void ParallelCoordinator::copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) {
    if (backend_) backend_->copyCanonical(entityIds, values, ownerMask);
}

FDM::DistributedIndexRange ParallelCoordinator::distributeIndices(
        std::int64_t localCount) {
    if (!decomposition_) return {0, localCount - 1, localCount};
    return decomposition_->distribute(localCount);
}

void ParallelCoordinator::broadcast(
        std::vector<double>& values, int root) {
    if (reduction_) reduction_->broadcast(values, root);
}

bool ParallelCoordinator::allRanksAgree(bool localValue) {
    return reduction_ ? reduction_->allRanksAgree(localValue) : localValue;
}

bool ParallelCoordinator::active() const {
    return backend_ && backend_->active();
}

int ParallelCoordinator::rank() const {
    return backend_ ? backend_->rank() : 0;
}

int ParallelCoordinator::size() const {
    return backend_ ? backend_->size() : 1;
}

void ParallelCoordinator::barrier() {
    if (reduction_) reduction_->barrier();
}

void ParallelCoordinator::configureBlocks(
        std::vector<MeshBlockField>* blocks,
        bool canonicalInterfaceFlux) {
    blocks_ = blocks;
    canonicalInterfaceFlux_ = canonicalInterfaceFlux;
}

} // namespace SF::Parallel
