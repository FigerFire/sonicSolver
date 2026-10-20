/// @file SF_saturationProperties.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"

#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

void addSaturationPropertyRates(const ModelContext& ctx,
                                std::vector<double>& mdot) {
    (void)ctx;
    (void)mdot;
    throw std::runtime_error(
        "SaturationProperties is a property closure, not a standalone "
        "mass-transfer source.");
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
