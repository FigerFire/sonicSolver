/// @file SF_surfaceTension.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_surfaceTension.h"

#include "SF_heaviside.h"

#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

Vector3 SurfaceTension::csfForce(const LevelSetField& levelSet,
                                 int i, int j, int k,
                                 double sigma,
                                 double epsilon) {
    if (!std::isfinite(sigma) || sigma < 0.0) {
        throw std::runtime_error(
            "SurfaceTension::csfForce: sigma must be finite and >= 0.");
    }
    const double delta = Heaviside::delta(levelSet.phi(i, j, k), epsilon);
    // phi>0 为液相，normal 从气相指向液相；该符号与
    // p_liquid-p_gas=-sigma*kappa 的 Young-Laplace 约定一致。
    const double factor = -sigma * levelSet.curvature(i, j, k) * delta;
    const Vector3& n = levelSet.normal(i, j, k);
    return Vector3(factor * n.x, factor * n.y, factor * n.z);
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
