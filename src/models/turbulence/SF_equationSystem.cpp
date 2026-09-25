/// @file SF_equationSystem.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

#include "SF_equationSystem.h"

#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "SF_equationModels.h"
#include "SF_equationModelOps.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace SF::Turbulence {
namespace {

void setupEquation(const Field& geometry,
                   TransportEquationState& equation,
                   const std::string& name) {
    equation.variable.setupLike(geometry, name);
    equation.previousConserved.setupLike(
        geometry, "previous." + name);
    equation.diffusivity.setupLike(
        geometry, "diffusivity." + name);
    equation.explicitSource.setupLike(
        geometry, "source." + name);
    equation.implicitSink.setupLike(
        geometry, "sink." + name);
}

void applyInitial(const Field& geometry,
                  const std::vector<BCSetting<double>>& settings,
                  ScalarField& field) {
    if (settings.empty()) {
        throw std::runtime_error(
            "Transported turbulence equation requires an explicit 0/ field.");
    }
    for (const auto& setting : settings) {
        if (setting.name == "all") {
            std::fill(
                field.values().begin(), field.values().end(),
                setting.value);
            continue;
        }
        const auto& indices = geometry.getSet(setting.name);
        for (int cell : indices) {
            field.values()[(size_t)cell] = setting.value;
        }
    }
}

void applyBoundary(const Field& geometry,
                   const std::vector<BCSetting<double>>& settings,
                   ScalarField& field) {
    if (settings.empty()) {
        throw std::runtime_error(
            "Transported turbulence equation requires boundaryField entries.");
    }
    const int ng = geometry.NG();
    const std::array<int,3> lower{ng, ng, ng};
    const std::array<int,3> upper{
        ng + geometry.NX() - 1,
        ng + geometry.NY() - 1,
        ng + geometry.NZ() - 1};
    for (const auto& setting : settings) {
        const int axis = setting.type == EMPTY
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(
                geometry, setting.name)
            : StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(
                geometry, setting.name);
        const auto& indices = geometry.getSet(setting.name);
        for (int cell : indices) {
            int i = 0, j = 0, k = 0;
            geometry.getIJK(cell, i, j, k);
            if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(
                    geometry, i, j, k)) {
                continue;
            }
            std::array<int,3> point{i,j,k};
            const int sign = point[(size_t)axis]
                == lower[(size_t)axis] ? -1 : 1;
            std::array<int,3> adjacent = point;
            adjacent[(size_t)axis] -= sign;
            const int adjacentCell = geometry.getIdx(
                adjacent[0], adjacent[1], adjacent[2]);
            if (setting.type == FIXED_VALUE) {
                field.values()[(size_t)cell] = setting.value;
            } else {
                field.values()[(size_t)cell] =
                    field.values()[(size_t)adjacentCell];
            }
            for (int layer = 1; layer <= ng; ++layer) {
                std::array<int,3> ghost = point;
                ghost[(size_t)axis] += sign * layer;
                const int ghostCell = geometry.getIdx(
                    ghost[0], ghost[1], ghost[2]);
                if (setting.type == FIXED_VALUE) {
                    std::array<int,3> mirror = point;
                    mirror[(size_t)axis] -= sign * layer;
                    const int mirrorCell = geometry.getIdx(
                        mirror[0], mirror[1], mirror[2]);
                    field.values()[(size_t)ghostCell] =
                        2.0 * setting.value
                        - field.values()[(size_t)mirrorCell];
                } else {
                    field.values()[(size_t)ghostCell] =
                        field.values()[(size_t)cell];
                }
            }
        }
    }
}

size_t phaseIndex(
        const Physics::PhaseSystems::PhaseSystem& system,
        const std::string& name) {
    for (size_t phase = 0; phase < system.phases().size(); ++phase) {
        if (Physics::Multiphase::normalizePhaseName(
                system.phases()[phase].name)
            == Physics::Multiphase::normalizePhaseName(name)) {
            return phase;
        }
    }
    throw std::runtime_error(
        "turbulenceProperties references unknown phase '" + name + "'.");
}

} // namespace

EquationSystem::EquationSystem(FDM::TurbulenceConfig config)
    : config_(std::move(config)) {}

void EquationSystem::initialize(
        const Physics::PhaseSystems::PhaseSystem& system) {
    active_ = config_.enabled
        && config_.family != FDM::TurbulenceFamily::None;
    if (!active_ || config_.family == FDM::TurbulenceFamily::DNS) return;
    if (config_.phaseNames.empty()) {
        throw std::runtime_error(
            "Eulerian turbulence requires explicit phases (...) in "
            "constant/turbulenceProperties.");
    }
    const bool kEpsilon =
        config_.family == FDM::TurbulenceFamily::RAS
        && config_.model == FDM::TurbulenceModelKind::kEpsilon;
    const bool kOmegaSST =
        config_.family == FDM::TurbulenceFamily::RAS
        && config_.model == FDM::TurbulenceModelKind::kOmegaSST;
    const bool smagorinsky =
        config_.family == FDM::TurbulenceFamily::LES
        && config_.model == FDM::TurbulenceModelKind::Smagorinsky;
    if (!kEpsilon && !kOmegaSST && !smagorinsky) {
        throw std::runtime_error(
            "Eulerian turbulence supports RAS/kEpsilon, RAS/kOmegaSST, "
            "LES/Smagorinsky, and DNS; the requested combination is not "
            "implemented and will not be substituted.");
    }
    const Field& geometry = system.geometry();
    for (const std::string& name : config_.phaseNames) {
        PhaseEquationState state;
        state.phaseIndex = phaseIndex(system, name);
        const std::string suffix = system.phases()[state.phaseIndex].name;
        state.eddyViscosity.setupLike(
            geometry, "mut." + suffix);
        if (kEpsilon || kOmegaSST) {
            setupEquation(
                geometry, state.kineticEnergy, "k." + suffix);
            applyInitial(
                geometry, config_.scalars.kInitial,
                state.kineticEnergy.variable);
        }
        if (kEpsilon) {
            setupEquation(
                geometry, state.dissipation, "epsilon." + suffix);
            applyInitial(
                geometry, config_.scalars.epsilonInitial,
                state.dissipation.variable);
        }
        if (kOmegaSST) {
            setupEquation(
                geometry, state.specificDissipation, "omega." + suffix);
            applyInitial(
                geometry, config_.scalars.omegaInitial,
                state.specificDissipation.variable);
            ModelOps::buildWallDistance(system, config_, state);
        }
        states_.push_back(std::move(state));
    }
    applyBoundary(system);
    prepare(system);
    commit(system);
    validate(system, "initialization");
}

bool EquationSystem::hasTransportEquations() const {
    return active_
        && config_.family == FDM::TurbulenceFamily::RAS
        && (config_.model == FDM::TurbulenceModelKind::kEpsilon
            || config_.model == FDM::TurbulenceModelKind::kOmegaSST);
}

std::string EquationSystem::description() const {
    if (!active_) return "Disabled";
    const std::string model=std::string(FDM::toString(config_.family))
        +"/"+FDM::toString(config_.model);
    if (config_.family == FDM::TurbulenceFamily::LES) {
        return model+" algebraic closure";
    }
    if (config_.family == FDM::TurbulenceFamily::DNS) {
        return model+" (no modeled turbulence equations)";
    }
    return model+" equations";
}

void EquationSystem::applyBoundary(
        const Physics::PhaseSystems::PhaseSystem& system) {
    if (!hasTransportEquations()) return;
    for (auto& state : states_) {
        ::SF::Turbulence::applyBoundary(
            system.geometry(), config_.scalars.kBoundary,
            state.kineticEnergy.variable);
        ScalarField* second = nullptr;
        const std::vector<BCSetting<double>>* secondBoundary = nullptr;
        if (config_.model == FDM::TurbulenceModelKind::kEpsilon) {
            second = &state.dissipation.variable;
            secondBoundary = &config_.scalars.epsilonBoundary;
        } else {
            second = &state.specificDissipation.variable;
            secondBoundary = &config_.scalars.omegaBoundary;
        }
        ::SF::Turbulence::applyBoundary(
            system.geometry(), *secondBoundary, *second);
    }
}

void EquationSystem::prepare(
        const Physics::PhaseSystems::PhaseSystem& system) {
    if (!active_ || config_.family == FDM::TurbulenceFamily::DNS) return;
    for (auto& state : states_) {
        if (config_.family == FDM::TurbulenceFamily::LES) {
            prepareSmagorinsky(system,config_,state);
        } else if (
                config_.model == FDM::TurbulenceModelKind::kEpsilon) {
            prepareKEpsilon(system,config_,state);
        } else {
            prepareKOmegaSST(system,config_,state);
        }
    }
}

void EquationSystem::commit(
        const Physics::PhaseSystems::PhaseSystem& system) {
    for (auto& state : states_) {
        const auto& mass = system.phases()[state.phaseIndex]
                               .primary.phaseMass.values();
        for (auto* equation : transportEquations(state)) {
            for (size_t cell = 0; cell < mass.size(); ++cell) {
                equation->previousConserved.values()[cell] =
                    mass[cell]*equation->variable.values()[cell];
            }
        }
    }
}

void EquationSystem::validate(
        const Physics::PhaseSystems::PhaseSystem& system,
        const std::string& stage) const {
    if (!hasTransportEquations()) return;
    const auto& geometry = system.geometry();
    const auto& c = config_.coefficients;
    const bool epsilonModel =
        config_.model == FDM::TurbulenceModelKind::kEpsilon;
    const double secondFloor =
        epsilonModel ? c.epsilonFloor : c.omegaFloor;
    const char* secondName = epsilonModel ? "epsilon" : "omega";
    for (const auto& state : states_) {
        const auto equations=transportEquations(state);
        const auto& kinetic=*equations[0];
        const auto& second=*equations[1];
        for (int k = geometry.NG();
             k < geometry.NG()+geometry.NZ(); ++k) {
            for (int j = geometry.NG();
                 j < geometry.NG()+geometry.NY(); ++j) {
                for (int i = geometry.NG();
                     i < geometry.NG()+geometry.NX(); ++i) {
                    const double kValue =
                        kinetic.variable(i,j,k);
                    const double secondValue =
                        second.variable(i,j,k);
                    if (!std::isfinite(kValue) || kValue < c.kFloor
                        || !std::isfinite(secondValue)
                        || secondValue < secondFloor) {
                        throw std::runtime_error(
                            std::string(FDM::toString(config_.model))
                            +" state invalid during " + stage
                            + " for phase '"
                            + system.phases()[state.phaseIndex].name
                            + "': k="+std::to_string(kValue)+", "
                            +secondName+"="+std::to_string(secondValue)+".");
                    }
                }
            }
        }
    }
}

PhaseEquationState* EquationSystem::find(size_t phase) {
    for (auto& state : states_) {
        if (state.phaseIndex == phase) return &state;
    }
    return nullptr;
}

const PhaseEquationState* EquationSystem::find(size_t phase) const {
    for (const auto& state : states_) {
        if (state.phaseIndex == phase) return &state;
    }
    return nullptr;
}

std::vector<TransportEquationState*> EquationSystem::transportEquations(
        PhaseEquationState& state) {
    if (!hasTransportEquations()) return {};
    if (config_.model == FDM::TurbulenceModelKind::kEpsilon) {
        return {&state.kineticEnergy,&state.dissipation};
    }
    return {&state.kineticEnergy,&state.specificDissipation};
}

std::vector<const TransportEquationState*>
EquationSystem::transportEquations(
        const PhaseEquationState& state) const {
    if (!hasTransportEquations()) return {};
    if (config_.model == FDM::TurbulenceModelKind::kEpsilon) {
        return {&state.kineticEnergy,&state.dissipation};
    }
    return {&state.kineticEnergy,&state.specificDissipation};
}

double EquationSystem::momentumEddyViscosity(
        size_t phase, int cell) const {
    if (!config_.coupleMomentum) return 0.0;
    const auto* state = find(phase);
    return state ? state->eddyViscosity.values()[(size_t)cell] : 0.0;
}

double EquationSystem::energyEddyDiffusivity(
        size_t phase, int cell) const {
    if (!config_.coupleEnergy) return 0.0;
    if (!std::isfinite(config_.turbulentPrandtl)
        || config_.turbulentPrandtl <= 0.0) {
        throw std::runtime_error(
            "turbulentPrandtl must be finite and positive.");
    }
    const auto* state = find(phase);
    return state
        ? state->eddyViscosity.values()[(size_t)cell]
            / config_.turbulentPrandtl
        : 0.0;
}

double EquationSystem::sourceTimeStep(
        const Physics::PhaseSystems::PhaseSystem& system,
        double sourceCfl) const {
    if (!hasTransportEquations()) {
        return std::numeric_limits<double>::max();
    }
    double result = std::numeric_limits<double>::max();
    for (const auto& state : states_) {
        const auto& mass = system.phases()[state.phaseIndex]
                               .primary.phaseMass.values();
        const Field& geometry = system.geometry();
        for (int k = geometry.NG();
             k < geometry.NG()+geometry.NZ(); ++k) {
            for (int j = geometry.NG();
                 j < geometry.NG()+geometry.NY(); ++j) {
                for (int i = geometry.NG();
                     i < geometry.NG()+geometry.NX(); ++i) {
                    const int cell = geometry.getIdx(i,j,k);
                    for (const auto* equation :
                         transportEquations(state)) {
                        const double primary =
                            mass[(size_t)cell]
                            *equation->variable.values()[(size_t)cell];
                        const double source =
                            equation->explicitSource.values()[(size_t)cell];
                        if (source > 0.0) {
                            result = std::min(
                                result,sourceCfl*primary/source);
                        }
                    }
                }
            }
        }
    }
    return result;
}

} // namespace SF::Turbulence
