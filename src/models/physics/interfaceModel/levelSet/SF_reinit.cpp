/// @file SF_reinit.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_reinit.h"

#include "SF_hjWeno.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

double distanceBetween(const Field& field,
                       int i0, int j0, int k0,
                       int i1, int j1, int k1) {
    const double dx = field.X(i1, j1, k1) - field.X(i0, j0, k0);
    const double dy = field.Y(i1, j1, k1) - field.Y(i0, j0, k0);
    const double dz = field.Z(i1, j1, k1) - field.Z(i0, j0, k0);
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!std::isfinite(d) || d <= 0.0) {
        throw std::runtime_error(
            "LevelSet reinit: invalid mesh spacing between "
            + cellText(i0, j0, k0) + " and " + cellText(i1, j1, k1)
            + ".");
    }
    return d;
}

double minimumSpacing(const Field& field) {
    double h = std::numeric_limits<double>::max();
    const int ng = field.NG();
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (HJWeno::activeAxis(field, 0)) {
                    h = std::min(h, distanceBetween(field, i - 1, j, k, i, j, k));
                    h = std::min(h, distanceBetween(field, i, j, k, i + 1, j, k));
                }
                if (HJWeno::activeAxis(field, 1)) {
                    h = std::min(h, distanceBetween(field, i, j - 1, k, i, j, k));
                    h = std::min(h, distanceBetween(field, i, j, k, i, j + 1, k));
                }
                if (HJWeno::activeAxis(field, 2)) {
                    h = std::min(h, distanceBetween(field, i, j, k - 1, i, j, k));
                    h = std::min(h, distanceBetween(field, i, j, k, i, j, k + 1));
                }
            }
        }
    }
    if (h == std::numeric_limits<double>::max()) {
        throw std::runtime_error(
            "LevelSet reinit: at least one active spatial direction is required.");
    }
    return h;
}

double sqr(double value) { return value * value; }

double godunovSelectPositive(const OneSidedDerivative& d) {
    const double backward = std::max(d.minus, 0.0);
    const double forward = std::min(d.plus, 0.0);
    return sqr(backward) >= sqr(forward) ? backward : forward;
}

double godunovSelectNegative(const OneSidedDerivative& d) {
    const double backward = std::min(d.minus, 0.0);
    const double forward = std::max(d.plus, 0.0);
    return sqr(backward) >= sqr(forward) ? backward : forward;
}

double selectedDerivative(const OneSidedDerivative& d, double signValue) {
    return signValue >= 0.0 ? godunovSelectPositive(d)
                            : godunovSelectNegative(d);
}

double gradientNorm(const Field& field,
                    const LevelSetField& levelSet,
                    const std::vector<double>& phi,
                    int i, int j, int k,
                    int order,
                    const HJWenoWeightOptions& weights,
                    double signValue) {
    double dXi = 0.0;
    double dEta = 0.0;
    double dZeta = 0.0;

    if (HJWeno::activeAxis(field, 0)) {
        dXi = selectedDerivative(
            HJWeno::derivative(levelSet, phi, i, j, k, 0, order, weights),
            signValue);
    }
    if (HJWeno::activeAxis(field, 1)) {
        dEta = selectedDerivative(
            HJWeno::derivative(levelSet, phi, i, j, k, 1, order, weights),
            signValue);
    }
    if (HJWeno::activeAxis(field, 2)) {
        dZeta = selectedDerivative(
            HJWeno::derivative(levelSet, phi, i, j, k, 2, order, weights),
            signValue);
    }

    const Vector3 grad =
        HJWeno::physicalGradient(field, i, j, k, dXi, dEta, dZeta);
    const double norm =
        std::sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
    if (!std::isfinite(norm)) {
        throw std::runtime_error(
            "LevelSet reinit: non-finite gradient norm at "
            + cellText(i, j, k) + ".");
    }
    return norm;
}

} // namespace

void Reinit::advance(const Field& field,
                     LevelSetField& levelSet,
                     const ReinitOptions& options) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSet reinit: LevelSetField dimensions do not match Field.");
    }
    if (options.pseudoSteps <= 0) {
        throw std::runtime_error(
            "LevelSet reinit: pseudoSteps must be > 0.");
    }
    if (!std::isfinite(options.pseudoTimeStep)
        || options.pseudoTimeStep <= 0.0) {
        throw std::runtime_error(
            "LevelSet reinit: pseudoTimeStep must be finite and > 0.");
    }
    if (!std::isfinite(options.wenoEpsilon) || options.wenoEpsilon <= 0.0
        || !std::isfinite(options.wenoPower) || options.wenoPower <= 0.0) {
        throw std::runtime_error(
            "LevelSet reinit requires explicit wenoEpsilon > 0 and wenoPower > 0.");
    }
    if (!std::isfinite(options.signSmoothingFactor)
        || options.signSmoothingFactor <= 0.0) {
        throw std::runtime_error(
            "LevelSet reinit requires explicit signSmoothingFactor > 0.");
    }
    HJWeno::requireStencil(field, options.order, "LevelSet reinit");

    const double h = minimumSpacing(field);
    const double signWidth = options.signSmoothingFactor * h;
    const HJWenoWeightOptions weights{
        options.wenoEpsilon, options.wenoPower};
    const std::vector<double> phi0 = levelSet.values();
    const int ng = field.NG();

    for (int step = 0; step < options.pseudoSteps; ++step) {
        if (options.prepareStage) options.prepareStage();
        const std::vector<double> oldPhi = levelSet.values();
        std::vector<double> newPhi = oldPhi;

        for (int k = ng; k < ng + field.NZ(); ++k) {
            for (int j = ng; j < ng + field.NY(); ++j) {
                for (int i = ng; i < ng + field.NX(); ++i) {
                    if (field.CellFlag(i, j, k) != FLUID_CELL) continue;

                    const int id = levelSet.getIdx(i, j, k);
                    const double base = phi0[(size_t)id];
                    const double signValue =
                        base / std::sqrt(base * base + signWidth * signWidth);
                    const double norm =
                        gradientNorm(field, levelSet, oldPhi,
                                     i, j, k, options.order, weights, signValue);
                    const double value =
                        oldPhi[(size_t)id]
                        - options.pseudoTimeStep * signValue * (norm - 1.0);
                    if (!std::isfinite(value)) {
                        throw std::runtime_error(
                            "LevelSet reinit: non-finite phi update at "
                            + cellText(i, j, k) + ".");
                    }
                    newPhi[(size_t)id] = value;
                }
            }
        }

        levelSet.values().swap(newPhi);
        levelSet.setMaterialPropertiesReady(false);
    }
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
