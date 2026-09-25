/// @file SF_turbulence.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_turbulence.h"

#include "SF_DNS.h"
#include "SF_LES.h"
#include "SF_RAS.h"
#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "core/mesh/SF_dimension.h"
#include "SF_kEpsilon.h"
#include "SF_kOmegaSST.h"
#include "SF_viscous.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>

namespace SF {
namespace Turbulence {

namespace {

std::unique_ptr<IModel> makeModel(const FDM::TurbulenceConfig& config) {
    if (!config.enabled) return nullptr;

    switch (config.family) {
        case FDM::TurbulenceFamily::RAS:
            if (config.model == FDM::TurbulenceModelKind::kOmegaSST) {
                return std::make_unique<RAS::KOmegaSSTModel>();
            }
            return std::make_unique<RAS::KEpsilonModel>();

        case FDM::TurbulenceFamily::LES:
            return std::make_unique<LES::SmagorinskyModel>();

        case FDM::TurbulenceFamily::DNS:
            return std::make_unique<DNS::DirectNumericalSimulationModel>();

        case FDM::TurbulenceFamily::None:
            break;
    }

    return nullptr;
}

int clampToInterior(int index, int lower, int upper) {
    return std::max(lower, std::min(index, upper));
}

void setFixedScalar(ScalarFields& state, ScalarSlot slot,
                    int i, int j, int k,
                    double bcValue,
                    int ng, int nx, int ny, int nz) {
    auto applyFixed = [&](int gI, int gJ, int gK, int rI, int rJ, int rK) {
        double realValue = state.slot(slot, rI, rJ, rK);
        state.slot(slot, gI, gJ, gK) = 2.0 * bcValue - realValue;
    };

    bool isGhost = (i < ng || i >= (nx + ng)
                 || j < ng || j >= (ny + ng)
                 || k < ng || k >= (nz + ng));

    if (isGhost) {
        int ri = clampToInterior(i, ng, nx + ng - 1);
        int rj = clampToInterior(j, ng, ny + ng - 1);
        int rk = clampToInterior(k, ng, nz + ng - 1);
        applyFixed(i, j, k, ri, rj, rk);
        return;
    }

    state.slot(slot, i, j, k) = bcValue;

    if (i == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(b, j, k, i, j, k);
    }
    if (i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i + b, j, k, i, j, k);
    }
    if (j == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(i, b, k, i, j, k);
    }
    if (j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i, j + b, k, i, j, k);
    }
    if (k == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(i, j, b, i, j, k);
    }
    if (k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i, j, k + b, i, j, k);
    }
}

void setZeroGradientScalar(ScalarFields& state, ScalarSlot slot,
                           int i, int j, int k,
                           int ng, int nx, int ny, int nz) {
    bool isGhost = (i < ng || i >= (nx + ng)
                 || j < ng || j >= (ny + ng)
                 || k < ng || k >= (nz + ng));

    if (isGhost) {
        int ri = clampToInterior(i, ng, nx + ng - 1);
        int rj = clampToInterior(j, ng, ny + ng - 1);
        int rk = clampToInterior(k, ng, nz + ng - 1);
        state.slot(slot, i, j, k) = state.slot(slot, ri, rj, rk);
        return;
    }

    int srcI = i;
    int srcJ = j;
    int srcK = k;
    if (i == ng && nx > 1) srcI = i + 1;
    else if (i == nx + ng - 1 && nx > 1) srcI = i - 1;
    if (j == ng && ny > 1) srcJ = j + 1;
    else if (j == ny + ng - 1 && ny > 1) srcJ = j - 1;
    if (k == ng && nz > 1) srcK = k + 1;
    else if (k == nz + ng - 1 && nz > 1) srcK = k - 1;

    state.slot(slot, i, j, k) = state.slot(slot, srcI, srcJ, srcK);

    if (i == ng) {
        for (int b = 0; b < ng; ++b) state.slot(slot, b, j, k) = state.slot(slot, i, j, k);
    }
    if (i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) state.slot(slot, i + b, j, k) = state.slot(slot, i, j, k);
    }
    if (j == ng) {
        for (int b = 0; b < ng; ++b) state.slot(slot, i, b, k) = state.slot(slot, i, j, k);
    }
    if (j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) state.slot(slot, i, j + b, k) = state.slot(slot, i, j, k);
    }
    if (k == ng) {
        for (int b = 0; b < ng; ++b) state.slot(slot, i, j, b) = state.slot(slot, i, j, k);
    }
    if (k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) state.slot(slot, i, j, k + b) = state.slot(slot, i, j, k);
    }
}

void setEmptyScalar(ScalarFields& state, ScalarSlot slot,
                    int i, int j, int k,
                    int ng, int nx, int ny, int nz, int axis) {
    int di = 0, dj = 0, dk = 0;
    if (axis == 0) di = 1;
    else if (axis == 1) dj = 1;
    else dk = 1;

    auto fill = [&](int gi, int gj, int gk) {
        state.slot(slot, gi, gj, gk) = state.slot(slot, i, j, k);
    };

    if (axis == 0 && i == ng) {
        for (int b = 1; b <= ng; ++b) fill(i - b * di, j, k);
    }
    if (axis == 0 && i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i + b * di, j, k);
    }
    if (axis == 1 && j == ng) {
        for (int b = 1; b <= ng; ++b) fill(i, j - b * dj, k);
    }
    if (axis == 1 && j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i, j + b * dj, k);
    }
    if (axis == 2 && k == ng) {
        for (int b = 1; b <= ng; ++b) fill(i, j, k - b * dk);
    }
    if (axis == 2 && k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i, j, k + b * dk);
    }
}

} // namespace

void ScalarFields::resizeLike(const Field& field) {
    mx_ = field.MX();
    my_ = field.MY();
    mz_ = field.MZ();
    ng_ = field.NG();
    totalSize_ = mx_ * my_ * mz_;

    k_.assign((size_t)totalSize_, 0.0);
    epsilon_.assign((size_t)totalSize_, 0.0);
    omega_.assign((size_t)totalSize_, 0.0);
    eddyMu_.assign((size_t)totalSize_, 0.0);
}

int ScalarFields::idx(int i, int j, int k) const {
    return (k * my_ + j) * mx_ + i;
}

double& ScalarFields::K(int i, int j, int k) { return k_[(size_t)idx(i, j, k)]; }
double ScalarFields::K(int i, int j, int k) const { return k_[(size_t)idx(i, j, k)]; }
double& ScalarFields::Epsilon(int i, int j, int k) { return epsilon_[(size_t)idx(i, j, k)]; }
double ScalarFields::Epsilon(int i, int j, int k) const { return epsilon_[(size_t)idx(i, j, k)]; }
double& ScalarFields::Omega(int i, int j, int k) { return omega_[(size_t)idx(i, j, k)]; }
double ScalarFields::Omega(int i, int j, int k) const { return omega_[(size_t)idx(i, j, k)]; }
double& ScalarFields::EddyMu(int i, int j, int k) { return eddyMu_[(size_t)idx(i, j, k)]; }
double ScalarFields::EddyMu(int i, int j, int k) const { return eddyMu_[(size_t)idx(i, j, k)]; }

double ScalarFields::clamped(ScalarSlot slot, int i, int j, int k) const {
    int ii = clampToInterior(i, ng_, mx_ - ng_ - 1);
    int jj = clampToInterior(j, ng_, my_ - ng_ - 1);
    int kk = clampToInterior(k, ng_, mz_ - ng_ - 1);
    return this->slot(slot, ii, jj, kk);
}

double& ScalarFields::slot(ScalarSlot slot, int i, int j, int k) {
    const size_t flatIndex = (size_t)idx(i, j, k);
    switch (slot) {
        case ScalarSlot::K: return k_[flatIndex];
        case ScalarSlot::Epsilon: return epsilon_[flatIndex];
        case ScalarSlot::Omega: return omega_[flatIndex];
        case ScalarSlot::EddyMu: return eddyMu_[flatIndex];
    }
    return k_[flatIndex];
}

double ScalarFields::slot(ScalarSlot slot, int i, int j, int k) const {
    const size_t flatIndex = (size_t)idx(i, j, k);
    switch (slot) {
        case ScalarSlot::K: return k_[flatIndex];
        case ScalarSlot::Epsilon: return epsilon_[flatIndex];
        case ScalarSlot::Omega: return omega_[flatIndex];
        case ScalarSlot::EddyMu: return eddyMu_[flatIndex];
    }
    return k_[flatIndex];
}

std::vector<double>& ScalarFields::values(ScalarSlot slot) {
    switch (slot) {
        case ScalarSlot::K: return k_;
        case ScalarSlot::Epsilon: return epsilon_;
        case ScalarSlot::Omega: return omega_;
        case ScalarSlot::EddyMu: return eddyMu_;
    }
    throw std::runtime_error("Unknown turbulence scalar slot.");
}

const std::vector<double>& ScalarFields::values(ScalarSlot slot) const {
    switch (slot) {
        case ScalarSlot::K: return k_;
        case ScalarSlot::Epsilon: return epsilon_;
        case ScalarSlot::Omega: return omega_;
        case ScalarSlot::EddyMu: return eddyMu_;
    }
    throw std::runtime_error("Unknown turbulence scalar slot.");
}

Manager::Manager(FDM::TurbulenceConfig config)
    : config_(std::move(config)) {}

bool Manager::initialize(const Field& field) {
    model_ = makeModel(config_);
    if (!model_) return false;

    state_.resizeLike(field);
    model_->initialize(field, state_, config_);
    return true;
}

bool Manager::active() const {
    return model_ != nullptr;
}

std::string Manager::description() const {
    if (!model_) return "None";
    return std::string(model_->familyName()) + "/" + model_->modelName();
}

void Manager::applyBoundary(const Field& field) {
    if (model_) model_->applyBoundary(field, state_, config_);
}

void Manager::correct(const Field& field, double dt) {
    if (model_) model_->correct(field, state_, config_, dt);
}

double Manager::dynamicViscosity(const Field& field,
                                 int i, int j, int k,
                                 double laminarMu) const {
    if (!model_) return laminarMu;
    int ii = clampToInterior(i, field.NG(), field.NG() + field.NX() - 1);
    int jj = clampToInterior(j, field.NG(), field.NG() + field.NY() - 1);
    int kk = clampToInterior(k, field.NG(), field.NG() + field.NZ() - 1);
    return laminarMu + std::max(model_->eddyDynamicViscosity(field, state_, config_, ii, jj, kk), 0.0);
}

std::vector<std::string> Manager::distributedReadFields() const {
    if (!model_ || config_.family != FDM::TurbulenceFamily::RAS) return {};
    if (config_.model == FDM::TurbulenceModelKind::kEpsilon) {
        return {"k", "epsilon", "mu_t"};
    }
    if (config_.model == FDM::TurbulenceModelKind::kOmegaSST) {
        return {"k", "omega", "mu_t"};
    }
    return {};
}

std::vector<std::string> Manager::distributedWriteFields() const {
    if (!model_) return {};
    auto result = distributedReadFields();
    if (result.empty()) result.push_back("mu_t");
    return result;
}

int Manager::distributedHaloDepth() const {
    return distributedReadFields().empty() ? 0 : 1;
}

void applyScalarInitialConditions(const Field& flow,
                                  const std::vector<BCSetting<double>>& settings,
                                  ScalarFields& state,
                                  ScalarSlot slot) {
    for (const auto& ic : settings) {
        const auto& indices = flow.getSet(ic.name);
        for (int flatIndex : indices) {
            int i, j, k;
            flow.getIJK(flatIndex, i, j, k);
            state.slot(slot, i, j, k) = ic.value;
        }
    }
}

void applyScalarBoundaryConditions(const Field& flow,
                                   const std::vector<BCSetting<double>>& settings,
                                   ScalarFields& state,
                                   ScalarSlot slot) {
    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (const auto& bc : settings) {
        auto applyAt = [&](int i, int j, int k, int axis) {

            switch (bc.type) {
                case FIXED_VALUE:
                    setFixedScalar(state, slot, i, j, k, bc.value, ng, nx, ny, nz);
                    break;
                case ZERO_GRADIENT:
                case SYMMETRY:
                    setZeroGradientScalar(state, slot, i, j, k, ng, nx, ny, nz);
                    break;
                case EMPTY:
                    setEmptyScalar(state, slot, i, j, k, ng, nx, ny, nz, axis);
                    break;
            }
        };

        const auto& allSets = flow.getAllSets();
        if (allSets.find(bc.name) == allSets.end()) {
            std::cerr << "[SF FATAL] Turbulence boundary zone '" << bc.name
                      << "' is not a mesh set." << std::endl;
            std::exit(1);
        }

        const auto& indices = flow.getSet(bc.name);
        const int axis = (bc.type == EMPTY)
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(flow, bc.name)
            : -1;
        for (int flatIndex : indices) {
            int i, j, k;
            flow.getIJK(flatIndex, i, j, k);
            applyAt(i, j, k, axis);
        }
    }
}

double clampPositive(double value, double floorValue) {
    return std::max(value, floorValue);
}

double strainRateMagnitude(const Field& flow, int i, int j, int k) {
    using Viscous::CentralOrder;

    // Use physical gradients (chain-rule with metric terms) instead of
    // raw computational-space derivatives, so the strain-rate tensor is
    // correct on stretched / curvilinear grids.
    const SymmTensor3 strain = Math::symm(Viscous::velocityGradientAt(
        flow, i, j, k, CentralOrder::SECOND));
    return std::sqrt(std::max(
        2.0 * Math::doubleDot(strain, strain), 0.0));
}

double characteristicFilterWidth(const Field& flow,
                                 int i, int j, int k,
                                 double filterScale) {
    double jac = std::abs(flow.Jac(i, j, k));
    double base = std::cbrt(1.0 / std::max(jac, 1.0e-12));
    return std::max(filterScale, 1.0e-12) * base;
}

double scalarTurbulentDiffusion(const Field& flow,
                                const ScalarFields& state,
                                ScalarSlot slot,
                                double invSigma,
                                int i, int j, int k) {
    double result = 0.0;
    double Jcell = 1.0 / std::max(std::abs(flow.Jac(i, j, k)), 1.0e-30);

    auto gammaAt = [&](int ii, int jj, int kk) -> double {
        return std::max(state.EddyMu(ii, jj, kk), 0.0) * invSigma;
    };

    auto phiAt = [&](int ii, int jj, int kk) -> double {
        return state.slot(slot, ii, jj, kk);
    };

    // ── XI direction ──
    if (Math::isDirectionActive(Math::XI)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i - 1, j, k, Math::XI, cfL);
        Math::faceMetrics(flow, i,     j, k, Math::XI, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i - 1, j, k) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i + 1, j, k));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i - 1, j, k));
        double fluxR = gamR * JgR * (phiAt(i + 1, j, k) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    // ── ETA direction ──
    if (Math::isDirectionActive(Math::ETA)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i, j - 1, k, Math::ETA, cfL);
        Math::faceMetrics(flow, i, j,     k, Math::ETA, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i, j - 1, k) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i, j + 1, k));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i, j - 1, k));
        double fluxR = gamR * JgR * (phiAt(i, j + 1, k) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    // ── ZETA direction ──
    if (Math::isDirectionActive(Math::ZETA)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i, j, k - 1, Math::ZETA, cfL);
        Math::faceMetrics(flow, i, j, k,     Math::ZETA, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i, j, k - 1) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i, j, k + 1));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i, j, k - 1));
        double fluxR = gamR * JgR * (phiAt(i, j, k + 1) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    return result;
}

} // namespace Turbulence
} // namespace SF
