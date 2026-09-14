/// @file SF_pipeline.cpp
/// @brief 按边界集合调度 algebraic 或 ILW 重建闭合。

#include "SF_pipeline.h"

#include <stdexcept>

namespace SF::Boundary {

void Pipeline::synchronize() {
    if (!services_.runtime) {
        throw std::runtime_error(
            "BoundaryPipeline requires ExecutionRuntime for halo freshness.");
    }
    services_.runtime->prepare({
        "boundary pipeline halo",
        {Execution::readHalo(
            "conservative", services_.conservativeHaloDepth)}});
}

void Pipeline::prepare(const std::vector<Field*>& fields,
                       double time, double dt) {
    if (services_.physical) {
        for (Field* field : fields) {
            if (field) services_.physical->apply(*field);
        }
    }
    synchronize();
    if (services_.immersed) {
        services_.immersed->apply(fields, time, dt);
        if (services_.runtime) {
            services_.runtime->finalize({
                "immersed boundary state",
                {Execution::writeOwned("conservative")}});
        }
    }
    if (services_.immersed && services_.synchronizeAfterImmersed) {
        synchronize();
    }
    for (Field* field : fields) {
        if (field) field->invalidateThermodynamicCache();
    }
}

} // namespace SF::Boundary
