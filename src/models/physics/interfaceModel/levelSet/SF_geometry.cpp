/// @file SF_geometry.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_geometry.h"

#include "SF_derivative.h"
#include "SF_hjWeno.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

std::string cellText(int i, int j, int k) {
    std::ostringstream os;
    os << "(" << i << "," << j << "," << k << ")";
    return os.str();
}

int offsetI(int axis, int offset) { return axis == 0 ? offset : 0; }
int offsetJ(int axis, int offset) { return axis == 1 ? offset : 0; }
int offsetK(int axis, int offset) { return axis == 2 ? offset : 0; }

bool signDiffers(double a, double b) {
    if (a == b) return false;
    return a == 0.0 || b == 0.0 || a * b < 0.0;
}

bool isInterfaceCandidate(const Field& field,
                          const LevelSetField& levelSet,
                          int i, int j, int k) {
    const double c = levelSet.phi(i, j, k);
    if (HJWeno::activeAxis(field, 0)
        && (signDiffers(c, levelSet.phi(i - 1, j, k))
            || signDiffers(c, levelSet.phi(i + 1, j, k)))) {
        return true;
    }
    if (HJWeno::activeAxis(field, 1)
        && (signDiffers(c, levelSet.phi(i, j - 1, k))
            || signDiffers(c, levelSet.phi(i, j + 1, k)))) {
        return true;
    }
    if (HJWeno::activeAxis(field, 2)
        && (signDiffers(c, levelSet.phi(i, j, k - 1))
            || signDiffers(c, levelSet.phi(i, j, k + 1)))) {
        return true;
    }
    return false;
}

void requireCentralStencil(const Field& field) {
    if (field.NG() < 1
        && (HJWeno::activeAxis(field, 0)
            || HJWeno::activeAxis(field, 1)
            || HJWeno::activeAxis(field, 2))) {
        throw std::runtime_error(
            "LevelSet geometry requires at least one ghost layer for central derivatives.");
    }
}

std::vector<double> centralWeights() {
    std::vector<double> weights;
    if (!Math::Derivative::finiteDifferenceWeights({-1.0, 0.0, 1.0},
                                                   1,
                                                   weights)) {
        throw std::runtime_error(
            "LevelSet geometry: failed to build central derivative weights.");
    }
    return weights;
}

double phiDerivative(const Field& field,
                     const LevelSetField& levelSet,
                     const std::vector<double>& weights,
                     int i, int j, int k,
                     int axis) {
    if (!HJWeno::activeAxis(field, axis)) return 0.0;
    const int offsets[3] = {-1, 0, 1};
    double value = 0.0;
    for (int n = 0; n < 3; ++n) {
        value += weights[(size_t)n]
               * levelSet.phi(i + offsetI(axis, offsets[n]),
                              j + offsetJ(axis, offsets[n]),
                              k + offsetK(axis, offsets[n]));
    }
    return value;
}

double normalComponent(const LevelSetField& levelSet,
                       int i, int j, int k,
                       int component) {
    const Vector3& n = levelSet.normal(i, j, k);
    if (component == 0) return n.x;
    if (component == 1) return n.y;
    return n.z;
}

double normalDerivative(const Field& field,
                        const LevelSetField& levelSet,
                        const std::vector<double>& weights,
                        int i, int j, int k,
                        int axis,
                        int component) {
    if (!HJWeno::activeAxis(field, axis)) return 0.0;
    const int offsets[3] = {-1, 0, 1};
    double value = 0.0;
    for (int n = 0; n < 3; ++n) {
        value += weights[(size_t)n]
               * normalComponent(
                     levelSet,
                     i + offsetI(axis, offsets[n]),
                     j + offsetJ(axis, offsets[n]),
                     k + offsetK(axis, offsets[n]),
                     component);
    }
    return value;
}

double divergenceOfNormal(const Field& field,
                          const LevelSetField& levelSet,
                          const std::vector<double>& weights,
                          int i, int j, int k) {
    const double dnxDxi = normalDerivative(field, levelSet, weights, i, j, k, 0, 0);
    const double dnxDeta = normalDerivative(field, levelSet, weights, i, j, k, 1, 0);
    const double dnxDzeta = normalDerivative(field, levelSet, weights, i, j, k, 2, 0);

    const double dnyDxi = normalDerivative(field, levelSet, weights, i, j, k, 0, 1);
    const double dnyDeta = normalDerivative(field, levelSet, weights, i, j, k, 1, 1);
    const double dnyDzeta = normalDerivative(field, levelSet, weights, i, j, k, 2, 1);

    const double dnzDxi = normalDerivative(field, levelSet, weights, i, j, k, 0, 2);
    const double dnzDeta = normalDerivative(field, levelSet, weights, i, j, k, 1, 2);
    const double dnzDzeta = normalDerivative(field, levelSet, weights, i, j, k, 2, 2);

    const double dnxDx =
        field.XiX(i, j, k) * dnxDxi
        + field.EtX(i, j, k) * dnxDeta
        + field.ZeX(i, j, k) * dnxDzeta;
    const double dnyDy =
        field.XiY(i, j, k) * dnyDxi
        + field.EtY(i, j, k) * dnyDeta
        + field.ZeY(i, j, k) * dnyDzeta;
    const double dnzDz =
        field.XiZ(i, j, k) * dnzDxi
        + field.EtZ(i, j, k) * dnzDeta
        + field.ZeZ(i, j, k) * dnzDzeta;
    return dnxDx + dnyDy + dnzDz;
}

} // namespace

void Geometry::compute(const Field& field,
                       LevelSetField& levelSet,
                       const GeometryOptions& options) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSet geometry: LevelSetField dimensions do not match Field.");
    }
    if (!std::isfinite(options.gradientTolerance)
        || options.gradientTolerance <= 0.0) {
        throw std::runtime_error(
            "LevelSet geometry: gradientTolerance must be finite and > 0.");
    }

    requireCentralStencil(field);
    const std::vector<double> weights = centralWeights();

    std::fill(levelSet.normals().begin(), levelSet.normals().end(), Vector3());
    std::fill(levelSet.curvatures().begin(), levelSet.curvatures().end(), 0.0);
    std::fill(levelSet.interfaceMask().begin(), levelSet.interfaceMask().end(), 0);

    const int ng = field.NG();
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (field.CellFlag(i, j, k) != FLUID_CELL) continue;

                const double dXi = phiDerivative(field, levelSet, weights, i, j, k, 0);
                const double dEta = phiDerivative(field, levelSet, weights, i, j, k, 1);
                const double dZeta = phiDerivative(field, levelSet, weights, i, j, k, 2);
                const Vector3 grad =
                    HJWeno::physicalGradient(field, i, j, k, dXi, dEta, dZeta);
                const double mag =
                    std::sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
                if (!std::isfinite(mag)) {
                    throw std::runtime_error(
                        "LevelSet geometry: non-finite gradient at "
                        + cellText(i, j, k) + ".");
                }
                if (mag <= options.gradientTolerance) continue;

                const int id = levelSet.getIdx(i, j, k);
                levelSet.normals()[(size_t)id] =
                    Vector3(grad.x / mag, grad.y / mag, grad.z / mag);
            }
        }
    }

    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (field.CellFlag(i, j, k) != FLUID_CELL) continue;
                if (!isInterfaceCandidate(field, levelSet, i, j, k)) continue;

                const Vector3& n = levelSet.normal(i, j, k);
                const double nMag =
                    std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
                if (!std::isfinite(nMag) || nMag <= options.gradientTolerance) {
                    throw std::runtime_error(
                        "LevelSet geometry: degenerate interface gradient at "
                        + cellText(i, j, k) + ", phi="
                        + std::to_string(levelSet.phi(i, j, k)) + ".");
                }

                const double kappa =
                    divergenceOfNormal(field, levelSet, weights, i, j, k);
                if (!std::isfinite(kappa)) {
                    throw std::runtime_error(
                        "LevelSet geometry: non-finite curvature at "
                        + cellText(i, j, k) + ".");
                }
                const int id = levelSet.getIdx(i, j, k);
                levelSet.curvatures()[(size_t)id] = kappa;
                levelSet.interfaceMask()[(size_t)id] = 1;
            }
        }
    }
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
