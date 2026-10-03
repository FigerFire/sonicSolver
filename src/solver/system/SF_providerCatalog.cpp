#include "core/system/SF_operationIds.h"
#include "SF_providerCatalog.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {
namespace {
bool supports(const ProviderDescriptor& provider,
              const ExecutableOperation& operation) {
    return !operation.requirements.empty()
        && std::all_of(operation.requirements.begin(),operation.requirements.end(),
            [&](OperationCapability need) {
                return std::find(provider.capabilities.begin(),
                    provider.capabilities.end(),need)!=provider.capabilities.end();
            });
}

} // namespace

void ProviderCatalog::add(ProviderDescriptor descriptor) {
    if (descriptor.id.empty())
        throw std::invalid_argument("Provider registration needs an implementation id.");
    if (std::any_of(descriptors_.begin(),descriptors_.end(),
            [&](const ProviderDescriptor& existing) { return existing.id==descriptor.id; }))
        throw std::invalid_argument("Duplicate provider registration: " + descriptor.id);
    descriptors_.push_back(std::move(descriptor));
}

ResolvedOperationBinding ProviderCatalog::resolve(
        const ExecutableOperation& operation,
        std::string_view selectedProvider) const {
    ResolvedOperationBinding result;
    result.operation=operation.operation;
    const auto found=std::find_if(descriptors_.begin(),descriptors_.end(),
        [&](const auto& descriptor) { return descriptor.id==selectedProvider; });
    if (selectedProvider.empty())
        result.reason="selected method supplies no runtime provider for '"+operation.operation+"'";
    else if (found==descriptors_.end())
        result.reason="selected runtime provider '"+std::string(selectedProvider)+"' is not registered";
    else if (!supports(*found,operation))
        result.reason="selected runtime provider '"+found->id+"' lacks capability for '"+operation.operation+"'";
    else {
        result.provider=found->id;
        result.status=BindingStatus::Resolved;
    }
    return result;
}

ProviderCatalog ProviderCatalog::builtIn() {
    ProviderCatalog catalog;
    catalog.add({"flow.turbulence-closure",{OperationCapability::SingleFluidTurbulenceClosure}});
    catalog.add({"flow.turbulence",{OperationCapability::SingleFluidTurbulenceTransport}});
    catalog.add({"flow.conservative",{
        OperationCapability::ConservativeExplicit,
        OperationCapability::PressureSchedule,
        OperationCapability::MomentumPredictor,
        OperationCapability::PressureCorrection,
        OperationCapability::PressureLinearSolve,
        OperationCapability::PressureBoundary,
        OperationCapability::VelocityCorrection,
        OperationCapability::FluxCorrection}});
    const std::vector<OperationCapability> constantPressure{
        OperationCapability::PressureSchedule,
        OperationCapability::FixedTimeIteration,
        OperationCapability::MomentumPredictor,
        OperationCapability::PressureCorrection,
        OperationCapability::PressureLinearSolve,
        OperationCapability::PressureBoundary,
        OperationCapability::VelocityCorrection,
        OperationCapability::FluxCorrection};
    catalog.add({"flow.pressure-operators",constantPressure});
    catalog.add({"flow.rhie-chow",constantPressure});
    catalog.add({"flow.eulerian-pressure",{
        OperationCapability::EulerianPhaseExecution,
        OperationCapability::PressureLinearSolve}});
    catalog.add({"ibm.constraint",{OperationCapability::ImmersedConstraint}});
    return catalog;
}
} // namespace SF::System
