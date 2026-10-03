#pragma once

/// @file SF_builtinState.h
/// @brief Solver-owned catalog of ordinary physical variables. Models still declare
/// their requirements; users declare only new variables, initialization and boundaries.
#include "core/system/SF_stateRegistry.h"
#include "core/system/SF_solveProgram.h"

namespace SF::System {

/// @brief Typed entries prevent field-name heuristics in composition and cell loops.
enum class BuiltinState {
    Density, Momentum, TotalEnergy, Velocity, Pressure, Temperature,
    Enthalpy, FaceFlux, VolumeFraction, TurbulentEnergy, SpecificDissipation,
    DynamicViscosity, KinematicViscosity, ThermalConductivity
};

/// @brief A catalog declaration has no execution/solver-family identity.
StateSymbol builtinState(BuiltinState variable);

/// Known metadata is independent of the active case and has no runtime storage.
class BuiltinStateCatalog {
public:
    BuiltinStateCatalog();
    bool contains(std::string_view id) const;
    const StateSymbol& at(std::string_view id) const;
    /// Activate only a requested base symbol; explicit module metadata wins.
    void require(StateRegistry& active,std::string_view id) const;
    /// Writable solution metadata; dependency metadata remains a closure/view.
    StateSymbol solution(std::string_view id) const;
private:
    std::vector<StateSymbol> symbols_;
};

/// HOW requests its base target, never qualified/versioned symbols or algorithms.
void requireTargetStates(StateRegistry& active,const ExecutionScope& scope);

} // namespace SF::System
