#include "SF_levelSetBinding.h"
#include "core/system/SF_operationIds.h"
#include <map>
#include <memory>
#include <stdexcept>
namespace SF::SolverAlgorithm {
void bindLevelSetOperations(Run::OpRegistry& operations,
                           const System::CompiledSolvePlan& plan,
                           State::StateBundle& state,
                           const FDM::SolverServices& services,
                           System::StateRealization& realized) {
    // This handle aliases the port's original frozen workspace; it owns no phi array.
    auto reference=std::make_shared<State::DistributedFieldView>();
    for (const auto& call:plan.compiledProgram.steps) {
        if (call.equationMethod=="LevelSetAdvection") {
            if (!services.levelSet) throw std::runtime_error("Missing level-set port for provider interface.level-set, operation explicit.stage.execute.");
            const auto& p=std::any_cast<const std::map<std::string,double>&>(call.providerContract);
            services.levelSet->validateAdvection((int)p.at("order"),p.at("epsilon"),p.at("power"),p.at("csf")!=0,p.at("ghostFluid")!=0,p.at("sigma"),p.at("width"));
        }
        if (call.backendProvider!="interface.level-set") continue;
        if (!services.levelSet) throw std::runtime_error("Missing level-set port for provider interface.level-set, operation "+call.backendOperation+".");
        if (call.backendOperation==System::OpIds::LevelSetReference) {
            const auto p=std::any_cast<std::map<std::string,double>>(call.providerContract);
            operations.bind(System::OpIds::LevelSetReference,"interface.level-set",[&services,&realized,reference,p] {
                services.levelSet->beginReinitialization((int)p.at("order"),p.at("pseudoDt"),
                    p.at("epsilon"),p.at("power"),p.at("signFactor"));
                *reference=services.levelSet->referenceView();
                realized.bindView("phi0",System::StateViewKind::Workspace,*reference);
            });
        } else if (call.backendOperation==System::OpIds::LevelSetPseudoStage)
            operations.bind(System::OpIds::LevelSetPseudoStage,"interface.level-set",[&services] {services.levelSet->reinitializeStage();});
        else if (call.backendOperation==System::OpIds::LevelSetGeometry)
            operations.bind(System::OpIds::LevelSetGeometry,"interface.level-set",[&services,&state] {services.levelSet->publishGeometry(state.dt);});
        else throw std::runtime_error("Unimplemented level-set operation '"+call.backendOperation+"' for provider interface.level-set.");
    }
}
}
