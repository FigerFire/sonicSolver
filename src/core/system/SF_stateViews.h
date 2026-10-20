#pragma once

/// @file SF_stateViews.h
/// @brief Compiler output: lazily requested versions of base STATE symbols.
#include "core/system/SF_stateRegistry.h"

namespace SF::System {
/// @brief Stage/OldTime are temporal lowering output, never user registrations.
enum class StateViewKind { Physical, Working, Correction, OldTime, Stage, Workspace };
/// @brief Only one storage authority supplies a view; aliases never own another Q.
enum class StateViewOwner { PhysicalState, CompilerWorkspace, NumericalProvider };

/// @brief No ordering data lives in the state-view allocation contract.
struct CompiledStateView {
    std::string symbol;
    StateViewKind kind=StateViewKind::Physical;
    int stage=-1;
    std::string storage;
    int componentOffset=0;
    int components=1;
    VariableLocation location=VariableLocation::EulerianCell;
    OwnershipKind ownership=OwnershipKind::EulerianGlobalDof;
    StateViewOwner owner=StateViewOwner::PhysicalState;
};

inline const char* toString(StateViewKind kind) {
    switch (kind) {
    case StateViewKind::Physical:return "Physical";
    case StateViewKind::Working:return "Working";
    case StateViewKind::Correction:return "Correction";
    case StateViewKind::OldTime:return "OldTime";
    case StateViewKind::Stage:return "Stage";
    case StateViewKind::Workspace:return "Workspace";
    }
    return "invalid";
}
} // namespace SF::System
