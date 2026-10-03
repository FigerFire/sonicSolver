#include "SF_frameProvider.h"

#include "SF_rotation.h"
#include "methods/numerics/structured/SF_iteration.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Physics::MRF {
namespace {

bool matches(const System::TermMatchContext& context,const char* unknown) {
    return context.expression.kind==System::FormulaExpr::Kind::Operator
        && context.expression.name=="source"
        && context.expression.arguments.size()==1
        && context.expression.arguments.front().name=="MRF"
        && context.output==unknown;
}

struct CompiledRotation {
    std::string zone;
    Vector3 center;
    Vector3 omega;
    Vector3 frameVelocity;
};

} // namespace

std::vector<System::SourceTermProviderDescriptor> termProviders(
        const std::vector<RotatingSetting>& settings) {
    return {
        {"source.mrf.conservative",
            [](const auto& context) { return matches(context,"rhoU"); },
            FDM::builtInSourceRecipe(FDM::SourceKind::MRF),{},
            "models.physics.mrf",
            [settings] {
                if (settings.empty()) throw std::runtime_error(
                    "MRF source is enabled but no rotating setting was loaded.");
                std::vector<CompiledRotation> compiled;
                compiled.reserve(settings.size());
                for (const auto& setting:settings)
                    compiled.push_back({setting.zone,setting.center,
                        Source::MRF::angularVelocity(setting),
                        setting.hasVelocity ? setting.velocity : Vector3()});
                return System::ConservativeSourceKernel{
                    [compiled=std::move(compiled)](
                            Field& field,Residual& residual) {
                        Math::forFluidInterior(field,[&](int i,int j,int k) {
                            for (const auto& rotation:compiled) {
                                const auto& zone=rotation.zone;
                                if (!(zone.empty() || zone=="all"
                                      || zone=="All" || zone=="ALL"
                                      || Source::MRF::appliesToZone(
                                          field,zone,i,j,k))) continue;
                                Source::Rotating::addRotatingFrame(
                                    field,residual,i,j,k,rotation.center,
                                    rotation.omega,rotation.frameVelocity);
                            }
                        });
                    }};
            }},
        {"source.mrf.primitive",
            [](const auto& context) { return matches(context,"U"); },
            FDM::builtInSourceRecipe(FDM::SourceKind::MRF),
            [settings] {
                if (settings.empty()) throw std::runtime_error(
                    "MRF source is enabled but no rotating setting was loaded.");
                std::vector<CompiledRotation> compiled;
                compiled.reserve(settings.size());
                for (const auto& setting:settings)
                    compiled.push_back({setting.zone,setting.center,
                        Source::MRF::angularVelocity(setting),
                        setting.hasVelocity ? setting.velocity : Vector3()});
                return System::PrimitiveMomentumSource{
                    [compiled=std::move(compiled)](
                            const Field& geometry,int i,int j,int k,
                            const Vector3& velocity) {
                        Vector3 acceleration;
                        const Vector3 position=
                            Source::MRF::cellPosition(geometry,i,j,k);
                        for (const auto& rotation:compiled)
                            if (Source::MRF::appliesToZone(
                                    geometry,rotation.zone,i,j,k))
                                acceleration=acceleration+
                                    Source::Rotating::rotatingFrameAcceleration(
                                        position,velocity,rotation.center,
                                        rotation.omega,rotation.frameVelocity);
                        return acceleration;
                    }};
            },"models.physics.mrf"}
    };
}

} // namespace SF::Physics::MRF
