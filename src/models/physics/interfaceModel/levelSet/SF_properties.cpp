/// @file SF_properties.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_properties.h"

#include "SF_heaviside.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

std::string phaseLabel(PhaseRole role) {
    return role == PhaseRole::Liquid ? "liquid" : "gas";
}

const PhaseProperties* findByRole(const MultiPhaseConfig& config,
                                  PhaseRole role) {
    for (const PhaseProperties& phase : config.phases) {
        if (phase.role == role) return &phase;
    }
    return nullptr;
}

void requirePhaseProperties(const PhaseProperties* phase,
                            PhaseRole role,
                            const char* property) {
    if (phase == nullptr) {
        throw std::runtime_error(
            std::string("LevelSet properties: missing ")
            + phaseLabel(role) + " phase declaration.");
    }
    if (std::string(property) == "density") {
        if (!std::isfinite(phase->density) || phase->density <= 0.0) {
            throw std::runtime_error(
                "LevelSet properties: phase '" + phase->name
                + "' density must be finite and > 0.");
        }
        return;
    }
    if (!std::isfinite(phase->viscosity) || phase->viscosity < 0.0) {
        throw std::runtime_error(
            "LevelSet properties: phase '" + phase->name
            + "' viscosity must be finite and >= 0.");
    }
}

double indicator(double phi, double epsilon) {
    return epsilon > 0.0 ? Heaviside::regularized(phi, epsilon)
                         : Heaviside::sharp(phi);
}

} // namespace

void Properties::update(const Field& field,
                        LevelSetField& levelSet,
                        const MultiPhaseConfig& config) {
    if (!levelSet.isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSet properties: LevelSetField dimensions do not match Field.");
    }

    const PhaseProperties* liquid = findByRole(config, PhaseRole::Liquid);
    const PhaseProperties* gas = findByRole(config, PhaseRole::Gas);
    requirePhaseProperties(liquid, PhaseRole::Liquid, "density");
    requirePhaseProperties(gas, PhaseRole::Gas, "density");
    requirePhaseProperties(liquid, PhaseRole::Liquid, "viscosity");
    requirePhaseProperties(gas, PhaseRole::Gas, "viscosity");

    const double epsilon = config.levelSet.interfaceThickness;
    if (epsilon < 0.0 || !std::isfinite(epsilon)) {
        throw std::runtime_error(
            "LevelSet properties: interfaceThickness must be finite and >= 0.");
    }

    for (int k = 0; k < field.MZ(); ++k) {
        for (int j = 0; j < field.MY(); ++j) {
            for (int i = 0; i < field.MX(); ++i) {
                const int id = levelSet.getIdx(i, j, k);
                const double h = indicator(levelSet.phi(i, j, k), epsilon);
                levelSet.densities()[(size_t)id] =
                    h * liquid->density + (1.0 - h) * gas->density;
                levelSet.viscosities()[(size_t)id] =
                    h * liquid->viscosity + (1.0 - h) * gas->viscosity;
            }
        }
    }

    levelSet.setMaterialPropertiesReady(true);
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
