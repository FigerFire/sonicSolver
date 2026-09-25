/// @file SF_advect.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.08.10-----------*/

#include "SF_advect.h"

#include "SF_hjWeno.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::Physics::Multiphase {
namespace {

std::string cellText(int i, int j, int k) {
    std::ostringstream os;
    os << "(" << i << "," << j << "," << k << ")";
    return os.str();
}

Vector3 velocityAt(const Field& field, int i, int j, int k) {
    const double rho = field(i, j, k, RHO);
    if (!std::isfinite(rho) || rho <= 0.0) {
        throw std::runtime_error(
            "LevelSet RHS: non-physical density at "
            + cellText(i, j, k) + ", rho=" + std::to_string(rho) + ".");
    }
    const Vector3 velocity(field(i, j, k, RU) / rho,
                           field(i, j, k, RV) / rho,
                           field(i, j, k, RW) / rho);
    if (!std::isfinite(velocity.x) || !std::isfinite(velocity.y)
        || !std::isfinite(velocity.z)) {
        throw std::runtime_error(
            "LevelSet RHS: non-finite stage velocity at "
            + cellText(i, j, k) + ".");
    }
    return velocity;
}

void validate(const Field& field,
              const LevelSetField& levelSet,
              const AdvectOptions& options) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSet RHS: state dimensions do not match Field.");
    }
    if (!std::isfinite(options.wenoEpsilon) || options.wenoEpsilon <= 0.0
        || !std::isfinite(options.wenoPower) || options.wenoPower <= 0.0) {
        throw std::runtime_error(
            "LevelSet RHS requires explicit wenoEpsilon > 0 and "
            "wenoPower > 0.");
    }
    if (!options.skipSolidCellsDeclared) {
        throw std::runtime_error(
            "LevelSet RHS requires an explicit non-fluid-cell policy.");
    }
    HJWeno::requireStencil(field, options.order, "LevelSet RHS");
}

} // namespace

void Advect::assembleRHS(const Field& field,
                         const LevelSetField& levelSet,
                         const AdvectOptions& options,
                         std::vector<double>& rhs) {
    validate(field, levelSet, options);
    rhs.assign((size_t)levelSet.TotalSize(), 0.0);
    const int ng = field.NG();
    const HJWenoWeightOptions weights{
        options.wenoEpsilon, options.wenoPower};

    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (options.skipSolidCells
                    && field.CellFlag(i, j, k) != FLUID_CELL) {
                    continue;
                }

                const Vector3 velocity = velocityAt(field, i, j, k);
                double hamiltonian = 0.0;
                for (int axis = 0; axis < 3; ++axis) {
                    if (!HJWeno::activeAxis(field, axis)) continue;
                    const OneSidedDerivative derivative = HJWeno::derivative(
                        levelSet, levelSet.values(), i, j, k,
                        axis, options.order, weights);
                    const double speed = HJWeno::contravariantSpeed(
                        field, i, j, k, axis, velocity);
                    hamiltonian += speed * HJWeno::upwind(derivative, speed);
                }

                const double value = -hamiltonian;
                if (!std::isfinite(value)) {
                    throw std::runtime_error(
                        "LevelSet RHS: non-finite value at "
                        + cellText(i, j, k) + ".");
                }
                rhs[(size_t)levelSet.getIdx(i, j, k)] = value;
            }
        }
    }
}

} // namespace SF::Physics::Multiphase
