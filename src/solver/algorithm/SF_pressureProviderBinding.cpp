#include "core/system/SF_operationIds.h"
/// @file SF_pressureProviderBinding.cpp
/// @brief One-time construction from a compiled operation binding.

#include "solver/algorithm/SF_pressureProviderBinding.h"

#include <stdexcept>
#include <utility>

namespace SF::SolverAlgorithm {

PressureProviderBinding bindPressureProvider(
        const System::ResolvedOperationBinding& binding,
        const System::CompiledNumericalSystem& numerics,
        const LegacyPressureInputs& legacy) {
    if (binding.operation!=System::OpIds::PressurePrepare
        || binding.status!=System::BindingStatus::Resolved)
        throw std::runtime_error("Pressure prepare has no resolved provider.");
    PressureProviderBinding result;
    if (binding.provider=="flow.pressure-operators") {
        if (!numerics.pressureOperator)
            throw std::runtime_error(
                "Resolved pressure operator has no compiled configuration.");
        result.constant=std::make_unique<Pressure::PressureOperators>(
            *numerics.pressureOperator,numerics);
    } else if (binding.provider=="flow.conservative") {
        result.conservative=std::make_unique<PressureBased::Corrector>(
            legacy.boundaries,legacy.pressure,legacy.idealGasGamma);
    } else {
        throw std::runtime_error(
            "No pressure operation binder for resolved provider '"
            +binding.provider+"'.");
    }
    return result;
}

} // namespace SF::SolverAlgorithm
