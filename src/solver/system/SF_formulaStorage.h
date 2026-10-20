#pragma once

/// @file SF_formulaStorage.h
/// @brief Non-owning, single-patch cell binding for compiled Equation methods.

#include "SF_stateRealizer.h"
#include "solver/discretization/SF_formulaOperator.h"

namespace SF::System {

/// Binds a declared Output to already-realized Field/ScalarField storage.
/// Multi-patch GlobalDof numbering remains unsupported by this generic path;
/// callers must supply the numerical boundary value closure explicitly.
FormulaValues bindFormulaCellValues(
    const StateRealization& state, const CompiledTarget& output,
    std::function<double(const std::string&,int,int,int,int)> boundaryValue);

} // namespace SF::System
