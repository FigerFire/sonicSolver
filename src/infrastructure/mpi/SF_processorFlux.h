#pragma once

/// @file SF_processorFlux.h
/// @brief processor interface 上唯一守恒数值通量的装配能力。

#include "SF_haloExchange.h"

namespace SF::Parallel {

class ProcessorFlux {
public:
    explicit ProcessorFlux(HaloExchange* halo = nullptr) : halo_(halo) {}
    void attach(HaloExchange* halo) { halo_ = halo; }

    void assemble(std::vector<MeshBlockField>& blocks,
                  const std::vector<FluxField*>& fluxes,
                  const std::vector<Residual*>& residuals) const {
        if (halo_) halo_->assembleCanonicalInterfaceFluxes(blocks, fluxes, residuals);
    }

private:
    HaloExchange* halo_ = nullptr;
};

} // namespace SF::Parallel
