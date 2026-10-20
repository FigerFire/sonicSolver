#pragma once

/// @file SF_physicsTemplate.h
/// @brief Physical-system preset provenance shared by composition and explain.

namespace SF::System {

/// Physical state family records preset provenance; it never selects a
/// runtime lifecycle.
enum class PhysicsTemplateKind {
    SingleFluid,
    HomogeneousMixture,
    OneFluidInterface,
    EulerianEulerian
};

const char* toString(PhysicsTemplateKind kind);

} // namespace SF::System
