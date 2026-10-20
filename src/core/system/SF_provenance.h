#pragma once

/// @file SF_provenance.h
/// @brief Composition origin is diagnostic metadata, never runtime dispatch.

#include <string>

namespace SF::System {

enum class OriginKind {
    BuiltinDefault,
    BuiltinPreset,
    Model,
    User,
    Generated
};

struct Provenance {
    OriginKind kind = OriginKind::BuiltinDefault;
    std::string source;
};

} // namespace SF::System
