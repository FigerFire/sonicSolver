/// @file SF_phaseMRF.cpp
/// @brief MRF 平动/旋转参考系源项模型实现。

#include "SF_phaseMRF.h"

#include "methods/numerics/structured/SF_structured.h"
#include "SF_rotation.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::Physics::PhaseSystems {
namespace {

class PhaseMRFSource final : public PhaseEquationSource {
public:
    explicit PhaseMRFSource(std::vector<RotatingSetting> settings)
        : settings_(std::move(settings)) {
        if (settings_.empty()) {
            throw std::runtime_error(
                "MRF is enabled but no rotating setting was loaded.");
        }
        for (const auto& setting : settings_) {
            const double axisMagnitude = std::sqrt(
                dot(setting.axis, setting.axis));
            if (!std::isfinite(setting.omega)
                || !std::isfinite(axisMagnitude)
                || axisMagnitude <= 0.0) {
                throw std::runtime_error(
                    "MRF setting requires finite omega and a non-zero axis.");
            }
        }
    }

    std::string name() const override { return "MRF"; }

    void add(
            const PhaseSystem& system,
            double,
            PhaseEquationSources& sources) const override {
        const Field& field = system.geometry();
        Math::forFluidInterior(field, [&](int i, int j, int k) {
            const int cell = field.getIdx(i,j,k);
            const Vector3 position(
                field.X(i,j,k), field.Y(i,j,k), field.Z(i,j,k));
            for (const auto& setting : settings_) {
                if (!Source::MRF::appliesToZone(
                        field, setting.zone, i, j, k)) {
                    continue;
                }
                const Vector3 omega =
                    Source::MRF::angularVelocity(setting);
                const Vector3 radius = position - setting.center;
                const Vector3 frameVelocity =
                    setting.hasVelocity ? setting.velocity : Vector3();
                for (size_t phase = 0;
                     phase < system.phases().size(); ++phase) {
                    const auto& state = system.phases()[phase];
                    const Vector3 velocity(
                        state.primitive.velocity[0].values()[(size_t)cell],
                        state.primitive.velocity[1].values()[(size_t)cell],
                        state.primitive.velocity[2].values()[(size_t)cell]);
                    const Vector3 relative = velocity - frameVelocity;
                    const Vector3 acceleration =
                        Source::Rotating::coriolisAcceleration(relative, omega)
                        + Source::Rotating::centrifugalAcceleration(
                            radius, omega);
                    const double mass =
                        state.primary.phaseMass.values()[(size_t)cell];
                    sources.momentum[phase][0].values()[(size_t)cell]
                        += mass * acceleration.x;
                    sources.momentum[phase][1].values()[(size_t)cell]
                        += mass * acceleration.y;
                    sources.momentum[phase][2].values()[(size_t)cell]
                        += mass * acceleration.z;
                    sources.energy[phase].values()[(size_t)cell]
                        += mass * dot(velocity, acceleration);
                }
            }
        });
    }

private:
    std::vector<RotatingSetting> settings_;
};

} // namespace

std::unique_ptr<PhaseEquationSource> makePhaseMRFSource(
        std::vector<RotatingSetting> settings) {
    return std::make_unique<PhaseMRFSource>(std::move(settings));
}

} // namespace SF::Physics::PhaseSystems
