/// @file SF_band.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_band.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

void NarrowBand::update(const Field& field,
                        LevelSetField& levelSet,
                        double width) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "NarrowBand::update: LevelSetField dimensions do not match Field.");
    }
    if (!std::isfinite(width) || width <= 0.0) {
        throw std::runtime_error("NarrowBand::update: width must be finite and > 0.");
    }

    std::fill(levelSet.narrowBandMask().begin(),
              levelSet.narrowBandMask().end(),
              (unsigned char)0);
    const int ng = field.NG();
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (field.CellFlag(i, j, k) != FLUID_CELL) continue;
                const int id = levelSet.getIdx(i, j, k);
                levelSet.narrowBandMask()[(size_t)id] =
                    std::abs(levelSet.phi(i, j, k)) <= width ? 1 : 0;
            }
        }
    }
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
