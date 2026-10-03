#include "SF_wallFluxProvider.h"
#include "SF_wallHeatSource.h"

#include <algorithm>

namespace SF::Physics::WallHeat {

System::SourceTermProviderDescriptor termProvider(
        const std::vector<WallHeatSetting>& settings) {
    return {"source.wallHeat.energy",
        [](const System::TermMatchContext& context) {
            return context.expression.kind==System::FormulaExpr::Kind::Operator
                && context.expression.name=="source"
                && context.expression.arguments.size()==1
                && context.expression.arguments.front().name=="wallHeat"
                && context.output=="rhoE";
        },FDM::builtInSourceRecipe(FDM::SourceKind::WallHeat),{},
        "models.physics.wallHeat",
        [settings] {
            return System::ConservativeSourceKernel{
                [settings](Field& field,Residual& residual) {
                    Source::WallHeat::addSource(field,residual,settings);
                }};
        }};
}

} // namespace SF::Physics::WallHeat
