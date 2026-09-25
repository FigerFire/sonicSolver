#pragma once

/// @file SF_boundaryPipeline.h
/// @brief Boundary and immersed-boundary ports used by solver execution.

#include "core/field/SF_field.h"
#include "SF_immersedConstraint.h"
#include "SF_immersedSystem.h"

#include <vector>

namespace SF::FDM {

class IBoundaryPipeline {
public:
    virtual ~IBoundaryPipeline() = default;
    virtual void prepare(
        const std::vector<Field*>& fields, double time, double dt) = 0;
};

class IImmersedBoundary {
public:
    virtual ~IImmersedBoundary() = default;
    virtual void apply(Field& field, double time, double dt) = 0;
    virtual void apply(
            const std::vector<Field*>& fields, double time, double dt) {
        for (Field* field : fields) if (field) apply(*field,time,dt);
    }
};

struct ImmersedCouplingPorts {
    IImmersedBoundary* boundary = nullptr;
    IImmersedConstraint* constraint = nullptr;
    IImmersedSystem* system = nullptr;

    bool hasBoundaryReconstruction() const { return boundary != nullptr; }
    bool hasVariationalConstraint() const {
        return constraint != nullptr || system != nullptr;
    }
};

} // namespace SF::FDM
