#include "core/system/SF_operationIds.h"
/// @file SF_pressureProviderBinding.cpp
/// @brief One-time construction from the selected method and actual resource contract.

#include "solver/algorithm/SF_pressureProviderBinding.h"

#include <algorithm>
#include <stdexcept>

namespace SF::SolverAlgorithm {

PressureProviderBinding bindPressureProvider(
        const System::Program& program,
        const LegacyPressureInputs& legacy) {
    const auto binding=std::find_if(program.runtime.operationBindings.begin(),program.runtime.operationBindings.end(),
        [](const auto& item){return item.operation==System::OpIds::PressurePrepare;});
    if (binding==program.runtime.operationBindings.end() || binding->status!=System::BindingStatus::Resolved)
        throw std::runtime_error("Pressure prepare has no resolved operation binding.");
    PressureProviderBinding result;
    if (const auto* inputs=System::requiredPressureInputs(program)) {
        result.constant=std::make_unique<Pressure::PressureOperators>(*inputs,program.numerics);
    } else if (binding->provider=="flow.conservative") {
        result.conservative=std::make_unique<PressureBased::Corrector>(
            legacy.boundaries,legacy.pressure,legacy.idealGasGamma);
    } else {
        throw std::runtime_error("Pressure operation binding has no supported resource contract: "+binding->operation);
    }
    return result;
}

} // namespace SF::SolverAlgorithm
