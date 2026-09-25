/// @file SF_state.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_state.h"
#include "SF_geometry.h"
#include "SF_hjWeno.h"
#include "SF_properties.h"
#include "SF_reinit.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

std::string bcTypeName(BCType type) {
    switch (type) {
    case FIXED_VALUE: return "FIXED_VALUE";
    case ZERO_GRADIENT: return "ZERO_GRADIENT";
    case SYMMETRY: return "SYMMETRY";
    case EMPTY: return "EMPTY";
    }
    return "UNKNOWN";
}

void requireFinitePhi(const std::string& context, double value) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(context + ": phi value must be finite.");
    }
}

int inwardNeighborIdx(const Field& field,
                      int i,
                      int j,
                      int k,
                      const std::string& setName) {
    const int i0 = field.NG();
    const int i1 = field.NG() + field.NX() - 1;
    const int j0 = field.NG();
    const int j1 = field.NG() + field.NY() - 1;
    const int k0 = field.NG();
    const int k1 = field.NG() + field.NZ() - 1;

    auto index = [&](int ii, int jj, int kk) {
        if (ii < 0 || ii >= field.MX()
            || jj < 0 || jj >= field.MY()
            || kk < 0 || kk >= field.MZ()) {
            return -1;
        }
        return field.getIdx(ii, jj, kk);
    };

    if (i <= i0 && field.NX() > 1) return index(i + 1, j, k);
    if (i >= i1 && field.NX() > 1) return index(i - 1, j, k);
    if (j <= j0 && field.NY() > 1) return index(i, j + 1, k);
    if (j >= j1 && field.NY() > 1) return index(i, j - 1, k);
    if (k <= k0 && field.NZ() > 1) return index(i, j, k + 1);
    if (k >= k1 && field.NZ() > 1) return index(i, j, k - 1);

    throw std::runtime_error(
        "LevelSet boundary: ZERO_GRADIENT set '" + setName
        + "' contains point (" + std::to_string(i) + ", "
        + std::to_string(j) + ", " + std::to_string(k)
        + ") that is not on a resolvable structured boundary.");
}

} // namespace

void LevelSetField::setupLike(const Field& field, double fillValue) {
    mx_ = field.MX();
    my_ = field.MY();
    mz_ = field.MZ();
    ng_ = field.NG();
    totalSize_ = field.TotalSize();

    phi_.setupLike(field, "phi", fillValue);
    normals_.assign((size_t)totalSize_, Vector3());
    curvature_.assign((size_t)totalSize_, 0.0);
    interfaceMask_.assign((size_t)totalSize_, 0);
    narrowBandMask_.assign((size_t)totalSize_, 0);
    density_.assign((size_t)totalSize_, 0.0);
    viscosity_.assign((size_t)totalSize_, 0.0);
    hasMaterialProperties_ = false;
}

bool LevelSetField::isCompatibleWith(const Field& field) const {
    return mx_ == field.MX()
        && my_ == field.MY()
        && mz_ == field.MZ()
        && ng_ == field.NG()
        && totalSize_ == field.TotalSize()
        && (int)phi_.values().size() == totalSize_
        && (int)normals_.size() == totalSize_
        && (int)curvature_.size() == totalSize_
        && (int)interfaceMask_.size() == totalSize_
        && (int)narrowBandMask_.size() == totalSize_
        && (int)density_.size() == totalSize_
        && (int)viscosity_.size() == totalSize_;
}

void LevelSetField::initializeFromSets(const Field& field,
                                       const MultiPhaseConfig& config) {
    validateMultiPhaseConfig(config, "multiPhase");

    const double baseSign = phaseSign(config, config.defaultPhase);
    setupLike(field, baseSign * std::abs(config.defaultSignedDistance));

    const auto& sets = field.getAllSets();
    for (const SetPhaseAssignment& assignment : config.setPhases) {
        const auto it = sets.find(assignment.setName);
        if (it == sets.end()) {
            throw std::runtime_error(
                "LevelSet initialization: set '" + assignment.setName
                + "' does not exist in Field.");
        }

        const double sign = phaseSign(config, assignment.phaseName);
        const double value = sign * std::abs(config.defaultSignedDistance);
        for (int idx : it->second) {
            if (idx < 0 || idx >= totalSize_) {
                throw std::runtime_error(
                    "LevelSet initialization: set '" + assignment.setName
                    + "' contains out-of-range index " + std::to_string(idx)
                    + ".");
            }
            phi_.values()[(size_t)idx] = value;
        }
    }

    // OpenFOAM internalField 先写全域，internalSets 中的解析 signed-distance
    // 随后覆盖对应 set；不能反过来把解析界面再次抹平。
    for (const LevelSetScalarCondition& condition
         : config.phiInitialConditions) {
        if (condition.type != FIXED_VALUE) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] set '"
                + condition.setName + "' uses " + bcTypeName(condition.type)
                + ". Initial phi supports only FIXED_VALUE; put "
                "ZERO_GRADIENT/EMPTY in 0/phi boundaryField.");
        }
        requireFinitePhi("LevelSet initialization set '" + condition.setName + "'",
                         condition.value);
        const auto it = sets.find(condition.setName);
        if (it == sets.end()) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] set '" + condition.setName
                + "' does not exist in Field.");
        }
        for (int idx : it->second) {
            if (idx < 0 || idx >= totalSize_) {
                throw std::runtime_error(
                    "LevelSet initialization: IC [phi] set '"
                    + condition.setName + "' contains out-of-range index "
                    + std::to_string(idx) + ".");
            }
            phi_.values()[(size_t)idx] = condition.value;
        }
    }

    for (const LevelSetPlaneInitializer& plane
         : config.phiPlaneInitializers) {
        const auto it = sets.find(plane.setName);
        if (it == sets.end()) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] signed-distance plane set '"
                + plane.setName + "' does not exist in Field.");
        }
        const double normalNorm =
            std::sqrt(plane.normal.x * plane.normal.x
                      + plane.normal.y * plane.normal.y
                      + plane.normal.z * plane.normal.z);
        if (!std::isfinite(normalNorm) || normalNorm <= 0.0) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] signed-distance plane set '"
                + plane.setName + "' has invalid normal.");
        }
        if (!std::isfinite(plane.scale) || plane.scale <= 0.0) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] signed-distance plane set '"
                + plane.setName + "' has invalid scale.");
        }
        const Vector3 n(plane.normal.x / normalNorm,
                        plane.normal.y / normalNorm,
                        plane.normal.z / normalNorm);
        for (int idx : it->second) {
            if (idx < 0 || idx >= totalSize_) {
                throw std::runtime_error(
                    "LevelSet initialization: IC [phi] signed-distance plane set '"
                    + plane.setName + "' contains out-of-range index "
                    + std::to_string(idx) + ".");
            }
            int i = 0;
            int j = 0;
            int k = 0;
            field.getIJK(idx, i, j, k);
            const double dx = field.X(i, j, k) - plane.point.x;
            const double dy = field.Y(i, j, k) - plane.point.y;
            const double dz = field.Z(i, j, k) - plane.point.z;
            const double value = plane.scale * (dx * n.x + dy * n.y + dz * n.z);
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "LevelSet initialization: IC [phi] signed-distance plane set '"
                    + plane.setName + "' produced non-finite phi.");
            }
            phi_.values()[(size_t)idx] = value;
        }
    }

    for (const LevelSetSphereInitializer& sphere
         : config.phiSphereInitializers) {
        const auto it = sets.find(sphere.setName);
        if (it == sets.end()) {
            throw std::runtime_error(
                "LevelSet initialization: IC [phi] signed-distance sphere set '"
                + sphere.setName + "' does not exist in Field.");
        }
        const double insideSign = phaseSign(config, sphere.insidePhase);
        for (int idx : it->second) {
            if (idx < 0 || idx >= totalSize_) {
                throw std::runtime_error(
                    "LevelSet initialization: IC [phi] signed-distance sphere set '"
                    + sphere.setName + "' contains out-of-range index "
                    + std::to_string(idx) + ".");
            }
            int i = 0, j = 0, k = 0;
            field.getIJK(idx, i, j, k);
            const double dx = field.X(i,j,k) - sphere.center.x;
            const double dy = field.Y(i,j,k) - sphere.center.y;
            const double dz = field.Z(i,j,k) - sphere.center.z;
            const double radiusSquared =
                (HJWeno::activeAxis(field, 0) ? dx*dx : 0.0)
                + (HJWeno::activeAxis(field, 1) ? dy*dy : 0.0)
                + (HJWeno::activeAxis(field, 2) ? dz*dz : 0.0);
            const double distance = std::sqrt(radiusSquared);
            const double value = insideSign * sphere.scale
                * (sphere.radius - distance);
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "LevelSet initialization: signed-distance sphere '"
                    + sphere.setName + "' produced non-finite phi.");
            }
            phi_.values()[(size_t)idx] = value;
        }
    }

    std::fill(normals_.begin(), normals_.end(), Vector3());
    std::fill(curvature_.begin(), curvature_.end(), 0.0);
    std::fill(interfaceMask_.begin(), interfaceMask_.end(), 0);
    std::fill(narrowBandMask_.begin(), narrowBandMask_.end(), 0);
    std::fill(density_.begin(), density_.end(), 0.0);
    std::fill(viscosity_.begin(), viscosity_.end(), 0.0);
    hasMaterialProperties_ = false;
}

void LevelSetField::applyBoundaryConditions(
        const Field& field,
        const MultiPhaseConfig& config) {
    if (!isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSet boundary: field size does not match LevelSetField.");
    }

    const auto& sets = field.getAllSets();
    auto applyCondition = [&](const LevelSetScalarCondition& condition) {
        const auto it = sets.find(condition.setName);
        if (it == sets.end()) {
            throw std::runtime_error(
                "LevelSet boundary: BC [phi] set '" + condition.setName
                + "' does not exist in Field.");
        }

        if (condition.type == FIXED_VALUE) {
            requireFinitePhi(
                "LevelSet boundary set '" + condition.setName + "'",
                condition.value);
        }

        for (int idx : it->second) {
            if (idx < 0 || idx >= totalSize_) {
                throw std::runtime_error(
                    "LevelSet boundary: BC [phi] set '" + condition.setName
                    + "' contains out-of-range index " + std::to_string(idx)
                    + ".");
            }

            switch (condition.type) {
            case FIXED_VALUE:
                phi_.values()[(size_t)idx] = condition.value;
                break;
            case ZERO_GRADIENT: {
                int i = 0;
                int j = 0;
                int k = 0;
                field.getIJK(idx, i, j, k);
                const int donor =
                    inwardNeighborIdx(field, i, j, k, condition.setName);
                if (donor < 0 || donor >= totalSize_) {
                    throw std::runtime_error(
                        "LevelSet boundary: ZERO_GRADIENT set '"
                        + condition.setName + "' produced invalid donor index "
                        + std::to_string(donor) + ".");
                }
                phi_.values()[(size_t)idx] =
                    phi_.values()[(size_t)donor];
                break;
            }
            case EMPTY:
                break;
            case SYMMETRY:
                throw std::runtime_error(
                    "LevelSet boundary: SYMMETRY for [phi] is not implemented; "
                    "use ZERO_GRADIENT or EMPTY explicitly.");
            }
        }
    };

    for (const LevelSetScalarCondition& condition
         : config.phiBoundaryConditions) {
        if (condition.type == FIXED_VALUE) continue;
        applyCondition(condition);
    }
    for (const LevelSetScalarCondition& condition
         : config.phiBoundaryConditions) {
        if (condition.type != FIXED_VALUE) continue;
        applyCondition(condition);
    }

    hasMaterialProperties_ = false;
}

void LevelSetField::computeGeometry(const Field& field,
                                    double gradientTolerance) {
    GeometryOptions options;
    options.gradientTolerance = gradientTolerance;
    Geometry::compute(field, *this, options);
}

void LevelSetField::reinitialize(const Field& field,
                                 int pseudoSteps,
                                 double pseudoTimeStep,
                                 int order,
                                 double wenoEpsilon,
                                 double wenoPower,
                                 double signSmoothingFactor,
                                 const std::function<void()>& prepareStage) {
    ReinitOptions options;
    options.pseudoSteps = pseudoSteps;
    options.pseudoTimeStep = pseudoTimeStep;
    options.order = order;
    options.wenoEpsilon = wenoEpsilon;
    options.wenoPower = wenoPower;
    options.signSmoothingFactor = signSmoothingFactor;
    options.prepareStage = prepareStage;
    Reinit::advance(field, *this, options);
}

void LevelSetField::updateMaterialProperties(const Field& field,
                                             const MultiPhaseConfig& config) {
    Properties::update(field, *this, config);
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
