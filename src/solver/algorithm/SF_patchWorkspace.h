/// @file SF_patchWorkspace.h
/// @brief 一个 density-based patch 在一次算法生命周期内使用的临时数值工作区。

#pragma once

#include "core/field/SF_field.h"
#include "core/flux/SF_flux.h"
#include "core/residual/SF_residual.h"
#include "solver/algorithm/SF_stageFlux.h"

#include <cstddef>

namespace SF::SolverAlgorithm {

/// @brief 不拥有物理 state 的 patch-local residual/face-flux scratch storage。
///
/// `SingleFluidStepper` 以 `StateBundle::patches` 的同一索引拥有此对象。
/// 它不进入 Field、StateBundle 或 ExecutionRuntime；后两者仅在同步阶段借用
/// 其中的具体 storage。
struct PatchWorkspace {
    PatchWorkspace()=default;
    PatchWorkspace(const PatchWorkspace&)=delete;
    PatchWorkspace& operator=(const PatchWorkspace&)=delete;
    PatchWorkspace(PatchWorkspace&& other) noexcept
        :convectiveFlux(std::move(other.convectiveFlux)),residual(std::move(other.residual)) {other.massFlux->invalidate();}
    PatchWorkspace& operator=(PatchWorkspace&& other) noexcept {
        if(this!=&other){massFlux->invalidate();other.massFlux->invalidate();convectiveFlux=std::move(other.convectiveFlux);residual=std::move(other.residual);}return *this;
    }
    FluxField convectiveFlux;
    Residual residual;
    std::shared_ptr<StageFluxPublication> massFlux=std::make_shared<StageFluxPublication>();

    void ensureFor(const Field& field) {
        const auto faceCount = static_cast<std::size_t>(3)
            * static_cast<std::size_t>(field.TotalSize());
        if (convectiveFlux.faceCount() != faceCount
            || convectiveFlux.variableCount() != field.NVar()) {
            massFlux->invalidate();
            convectiveFlux.setup(faceCount, field.NVar());
        }
        if (!residual.isCompatibleWith(
                field.MX(), field.MY(), field.MZ(), field.NVar())) {
            residual.setup(field.MX(), field.MY(), field.MZ(), field.NVar());
        }
    }

    void clear() {
        massFlux->invalidate();
        convectiveFlux.clear();
        residual.clear();
    }
};

} // namespace SF::SolverAlgorithm
