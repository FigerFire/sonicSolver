#include "SF_accelerationProvider.h"

#include "models/physics/mrf/SF_rotation.h"
#include "SF_gravity.h"
#include "methods/numerics/structured/SF_iteration.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Physics::Gravity {
namespace {

bool matches(const System::TermMatchContext& context,const char* unknown) {
    return context.expression.kind==System::FormulaExpr::Kind::Operator
        && context.expression.name=="source"
        && context.expression.arguments.size()==1
        && context.expression.arguments.front().name=="gravity"
        && context.output==unknown;
}

} // namespace

std::vector<System::SourceTermProviderDescriptor> termProviders(
        const std::vector<ZoneVectorSetting>& settings) {
    return {
        {"source.gravity.conservative",
            [](const auto& context) { return matches(context,"rhoU"); },
            FDM::builtInSourceRecipe(FDM::SourceKind::Gravity),{},
            "models.physics.gravity",
            [settings] {
                if (settings.empty()) throw std::runtime_error(
                    "Gravity source is enabled but no acceleration setting was loaded.");
                return System::ConservativeSourceKernel{
                    [settings](Field& field,Residual& residual) {
                        Math::forFluidInterior(field,[&](int i,int j,int k) {
                            for (const auto& setting:settings) {
                                const auto& zone=setting.zone;
                                if (zone.empty() || zone=="all" || zone=="All"
                                    || zone=="ALL"
                                    || Source::MRF::appliesToZone(
                                        field,zone,i,j,k))
                                    Source::Gravity::addTranslation(
                                        field,residual,i,j,k,setting.value);
                            }
                        });
                    }};
            }},
        {"source.gravity.primitive",
            [](const auto& context) { return matches(context,"U"); },
            FDM::builtInSourceRecipe(FDM::SourceKind::Gravity),
            [settings] {
                if (settings.empty()) throw std::runtime_error(
                    "Gravity source is enabled but no acceleration setting was loaded.");
                return System::PrimitiveMomentumSource{
                    [settings](const Field& geometry,int i,int j,int k,
                               const Vector3&) {
                        Vector3 acceleration;
                        for (const auto& setting:settings)
                            if (Source::MRF::appliesToZone(
                                    geometry,setting.zone,i,j,k))
                                acceleration=acceleration+setting.value;
                        return acceleration;
                    }};
            },"models.physics.gravity"}
    };
}

} // namespace SF::Physics::Gravity
