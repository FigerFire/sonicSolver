#pragma once

/// @file SF_pipeline.h
/// @brief fitted/IBM geometry 与 algebraic/ILW 重构的统一边界管线。

#include "SF_applicator.h"
#include "core/interfaces/SF_boundaryPipeline.h"
#include "core/interfaces/SF_executionRuntime.h"

namespace SF::Boundary {

/// @brief 将物理边界、halo 同步和 IBM ghost 闭合按稳定顺序组合。
///
/// `Applicator` 内部按配置选择 algebraic 或 ILW 物理边界重构；IBM 服务拥有
/// body geometry 和 ghost-cell 重构。求解器只调用 `prepare()`。
class Pipeline final : public FDM::IBoundaryPipeline {
public:
    struct Services {
        Applicator* physical = nullptr;
        FDM::IExecutionRuntime* runtime = nullptr;
        FDM::IImmersedBoundary* immersed = nullptr;
        bool synchronizeAfterImmersed = false;
        int conservativeHaloDepth = 1;
    };

    explicit Pipeline(Services services) : services_(services) {}

    /// @brief 应用 physical BC -> halo -> IBM，并按配置刷新 IBM 后 halo。
    void prepare(const std::vector<Field*>& fields,
                 double time, double dt) override;

private:
    Services services_;
    void synchronize();
};

} // namespace SF::Boundary
