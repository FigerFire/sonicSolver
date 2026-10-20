#pragma once

/// @file SF_transformationTypes.h
/// @brief Formulation request/status 的中立值对象，不执行 transformation。

#include "core/system/SF_equationIR.h"

namespace SF::System {

enum class TransformationState {
    NotRegistered,
    RegisteredAndActive,
    RegisteredButNotApplicable,
    RegisteredButInvalid,
    Applied
};

struct TransformationDescriptor {
    std::string id;
    std::string name;
    int priority = 0;
    bool explicitlyRequested = false;
    Provenance origin;
};

struct TransformationMatch {
    TransformationState state = TransformationState::RegisteredButNotApplicable;
    std::string reason;
};

struct TransformationRecord {
    TransformationDescriptor descriptor;
    TransformationState state = TransformationState::NotRegistered;
    std::string reason;
    std::vector<std::string> generatedEquations;
    std::vector<std::string> generatedOperators;
    std::vector<std::string> generatedOperations;
};

const char* toString(TransformationState state);

} // namespace SF::System
