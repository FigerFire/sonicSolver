#pragma once

/// @file SF_sourceTerm.h
/// @brief Execute already-compiled conservative source terms in declaration order.

#include "core/interfaces/SF_termKernel.h"
#include "core/residual/SF_residual.h"

#include <stdexcept>
#include <vector>

namespace SF::SourceTerm {

inline void Sp(
        Field& field,Residual& residual,
        const std::vector<System::ConservativeSourceKernel>& sources) {
    // One RHS stage has one source write-set; model kernels only contribute.
    residual.clearSource();
    for (const auto& source:sources) {
        if (!source)
            throw std::runtime_error("Compiled conservative source kernel is empty.");
        source(field,residual);
    }
}

} // namespace SF::SourceTerm
