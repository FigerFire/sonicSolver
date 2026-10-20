/// @file SF_phaseGravity.cpp
/// @brief 重力源项物理模型实现。

#include "SF_phaseGravity.h"

#include "methods/numerics/structured/SF_structured.h"
#include "SF_sourceAssembly.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::Physics::PhaseSystems {
namespace {

class PhaseGravitySource final : public PhaseEquationSource {
public:
    explicit PhaseGravitySource(std::vector<ZoneVectorSetting> settings)
        : settings_(std::move(settings)) {
        if (settings_.empty()) {
            throw std::runtime_error(
                "Gravity is enabled but no acceleration setting was loaded.");
        }
        for (const auto& setting : settings_) {
            if (!std::isfinite(setting.value.x)
                || !std::isfinite(setting.value.y)
                || !std::isfinite(setting.value.z)) {
                throw std::runtime_error(
                    "Gravity acceleration contains a non-finite component.");
            }
        }
    }

    std::string name() const override { return "Gravity"; }

    void add(
            const PhaseSystem& system,
            double,
            PhaseEquationSources& sources) const override {
        const Field& field = system.geometry();
        Math::forFluidInterior(field, [&](int i, int j, int k) {
            const int cell = field.getIdx(i,j,k);
            for (const auto& setting : settings_) {
                if (!Source::MRF::appliesToZone(
                        field, setting.zone, i, j, k)) {
                    continue;
                }
                for (size_t phase = 0;
                     phase < system.phases().size(); ++phase) {
                    const auto& state = system.phases()[phase];
                    const double mass =
                        state.primary.phaseMass.values()[(size_t)cell];
                    const std::array<double, 3> acceleration{
                        setting.value.x, setting.value.y, setting.value.z};
                    double power = 0.0;
                    for (int component = 0; component < 3; ++component) {
                        sources.momentum[phase][(size_t)component]
                            .values()[(size_t)cell] +=
                            mass * acceleration[(size_t)component];
                        power += state.primitive.velocity[(size_t)component]
                                     .values()[(size_t)cell]
                               * acceleration[(size_t)component];
                    }
                    sources.energy[phase].values()[(size_t)cell] += mass * power;
                }
            }
        });
    }

private:
    std::vector<ZoneVectorSetting> settings_;
};

} // namespace

std::unique_ptr<PhaseEquationSource> makePhaseGravitySource(
        std::vector<ZoneVectorSetting> settings) {
    return std::make_unique<PhaseGravitySource>(std::move(settings));
}

} // namespace SF::Physics::PhaseSystems
