#include "SF_transportBinding.h"
#include "core/system/SF_operationIds.h"
#include <stdexcept>
namespace SF::SolverAlgorithm {
namespace {
void correctTransport(const std::vector<Field*>& fields, State::StateBundle& state,
                      const FDM::SolverServices& services) {
    if (!services.transportModel) return;
    const auto reads = services.transportModel->distributedReadFields();
    const auto writes = services.transportModel->distributedWriteFields();
    const int haloDepth = services.transportModel->distributedHaloDepth();
    auto writeContract = [&](const char* name) {
        std::vector<Execution::FieldAccess> accesses;
        for (const auto& fieldName : writes) {
            accesses.push_back(Execution::writeOwned(fieldName));
        }
        if (services.executionRuntime && !accesses.empty()) {
            services.executionRuntime->finalize({name, accesses});
        }
    };
    auto readContract = [&](const char* name) {
        if (reads.empty()) return;
        if (haloDepth <= 0) {
            throw std::runtime_error(
                "Transport model declared distributed reads without "
                "a positive stencil halo depth.");
        }
        std::vector<Execution::FieldAccess> accesses;
        for (const auto& fieldName : reads) {
            accesses.push_back(Execution::readHalo(fieldName, haloDepth));
        }
        if (services.executionRuntime) {
            services.executionRuntime->prepare({name, accesses});
        }
    };
    for (Field* field : fields) {
        if (!field) continue;
        services.transportModel->applyBoundary(*field);
        writeContract("transport boundary update");
        readContract("transport correction stencil");
        services.transportModel->correct(*field, state.dt);
        services.transportModel->applyBoundary(*field);
        writeContract("transport corrected state");
        readContract("transport diffusion stencil");
    }
}
}
void bindTransportOperations(Run::OpRegistry& operations,
                            const System::RuntimeRequirements& requirements,
                            State::StateBundle& state,
                            const FDM::SolverServices& services) {
    for (const auto& binding:requirements.operationBindings) {
        if (binding.status!=System::BindingStatus::Resolved) continue;
        if (binding.provider!="flow.turbulence" && binding.provider!="flow.turbulence-closure") continue;
        if (!services.transportModel)
            throw std::runtime_error("Missing frozen turbulence closure service for provider '"+binding.provider+"', operation '"+binding.operation+"'.");
        if (binding.operation==System::OpIds::TurbulenceAdvance) {
            if (state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
                throw std::runtime_error("flow.turbulence requires its bound transport service and one serial patch.");
        } else if (binding.operation!=System::OpIds::TurbulenceClosureRefresh)
            throw std::runtime_error("Unimplemented transport operation '"+binding.operation+"' for provider '"+binding.provider+"'.");
        // Bundle and borrowed services outlive the registry; dt is read at invocation.
        operations.bind(binding.operation,binding.provider,[&state,&services] {
            correctTransport(state.patches,state,services);
        });
    }
}
}
