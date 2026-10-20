/// @file SF_mass.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_mass.h"

#include "SF_heaviside.h"

#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

double MassCorrection::liquidIndicatorSum(const Field& field,
                                          const LevelSetField& levelSet,
                                          double epsilon) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "MassCorrection::liquidIndicatorSum: LevelSetField dimensions do not match Field.");
    }
    if (!std::isfinite(epsilon) || epsilon < 0.0) {
        throw std::runtime_error(
            "MassCorrection::liquidIndicatorSum: epsilon must be finite and >= 0.");
    }

    double sum = 0.0;
    const int ng = field.NG();
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (field.CellFlag(i, j, k) != FLUID_CELL) continue;
                const double phi = levelSet.phi(i, j, k);
                sum += epsilon > 0.0 ? Heaviside::regularized(phi, epsilon)
                                     : Heaviside::sharp(phi);
            }
        }
    }
    return sum;
}

void MassCorrection::applyGlobalShift(Field&,
                                      LevelSetField&,
                                      double,
                                      double) {
    throw std::runtime_error(
        "MassCorrection::applyGlobalShift is an explicit placeholder. "
        "Global Level Set mass correction has not been implemented yet.");
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
