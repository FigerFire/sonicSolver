/// @file SF_fields.cpp
/// @brief 从 active simulation state 构造 ResultWriter 的 field views。
///
/// Data flow:
///   mesh/model state references
///       -> scalar callbacks and piece ownership metadata
///       -> ResultWriter registration
///
/// Views 不复制 field storage；本文件不决定输出时机或求解行为。

#include "app/application/output/SF_fields.h"

#include "SF_MultiBlockMesh.h"
#include "SF_equationSystem.h"
#include "SF_heaviside.h"
#include "SF_interfaceModel.h"
#include "SF_multiphase.h"
#include "SF_phaseSystem.h"
#include "levelSet/SF_state.h"
#include "core/state/SF_state.h"

#include <stdexcept>

namespace SF::Application::Output {

std::vector<ResultWriter::PieceLayout> buildMultiBlockVTKPieces(
        const MultiBlockMesh& mesh) {
    std::vector<ResultWriter::PieceLayout> pieces;
    pieces.reserve(mesh.size());
    for (size_t blockId = 0; blockId < mesh.size(); ++blockId) {
        ResultWriter::PieceLayout piece;
        piece.blockId = (int)blockId;
        piece.ownerRank = mesh.block(blockId).ownerRank;
        piece.partitionId = mesh.block(blockId).ownerRank;
        pieces.push_back(piece);
    }
    return pieces;
}

bool needsMultiFieldMPI(const MultiBlockMesh& mesh, int mpiSize) {
    if ((int)mesh.size() != mpiSize) return true;
    for (size_t blockId = 0; blockId < mesh.size(); ++blockId) {
        if (mesh.block(blockId).ownerRank != (int)blockId) return true;
    }
    return false;
}

std::vector<ResultWriter::ScalarField> interfaceVTKScalars(
        Physics::InterfaceModels::Model& model,
        const State::VariableRegistry* registry) {
    namespace MP = Physics::Multiphase;
    std::vector<ResultWriter::ScalarField> scalars;
    if (!model.initialized()) return scalars;
    auto* levelSet = model.levelSetState();
    if (model.representation()
            != Physics::InterfaceModels::Representation::LevelSet
        || levelSet == nullptr) {
        throw std::runtime_error(
            "Interface VTK output currently requires a Level Set model.");
    }
    auto* modelPtr = &model;
    auto* state = levelSet;
    if (registry) {
        for (const auto& variable : registry->variables()) {
            auto* scalar = variable.value;
            const std::string name = variable.descriptor.name;
            scalars.push_back({name, [scalar](int i, int j, int k) {
                return (*scalar)(i, j, k);
            }});
        }
    } else {
        scalars.push_back({"phi", [state](int i, int j, int k) {
            return state->phi(i, j, k);
        }});
    }
    const double epsilon = model.config().levelSet.interfaceThickness;
    scalars.push_back({"LiquidFraction", [state, epsilon](int i, int j, int k) {
        const double phi = state->phi(i, j, k);
        return epsilon > 0.0 ? MP::Heaviside::regularized(phi, epsilon)
                             : MP::Heaviside::sharp(phi);
    }});
    scalars.push_back({"GasFraction", [state, epsilon](int i, int j, int k) {
        const double phi = state->phi(i, j, k);
        const double liquid = epsilon > 0.0
            ? MP::Heaviside::regularized(phi, epsilon)
            : MP::Heaviside::sharp(phi);
        return 1.0 - liquid;
    }});
    scalars.push_back({"LevelSetCurvature", [state](int i, int j, int k) {
        return state->curvature(i, j, k);
    }});
    scalars.push_back({"PhaseDensity", [state](int i, int j, int k) {
        return state->hasMaterialProperties()
            ? state->density(i, j, k) : 0.0;
    }});
    scalars.push_back({"PhaseViscosity", [state](int i, int j, int k) {
        return state->hasMaterialProperties()
            ? state->viscosity(i, j, k) : 0.0;
    }});
    scalars.push_back({"LevelSetBand", [state](int i, int j, int k) {
        return state->inNarrowBand(i, j, k) ? 1.0 : 0.0;
    }});
    scalars.push_back({"PhaseIndicatorIntegral", [modelPtr](int, int, int) {
        return modelPtr->conservation().currentPhaseIndicator;
    }});
    scalars.push_back({"PhaseIndicatorRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().phaseIndicatorRelativeError;
    }});
    scalars.push_back({"PhaseMassIntegral", [modelPtr](int, int, int) {
        return modelPtr->conservation().currentPhaseMass;
    }});
    scalars.push_back({"PhaseMassRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().phaseMassRelativeError;
    }});
    scalars.push_back({"FlowMassRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().flowMassRelativeError;
    }});
    return scalars;
}

std::vector<ResultWriter::ScalarField> multiPhaseVTKScalars(
        Physics::Multiphase::MultiPhaseModel& model,
        const State::VariableRegistry* registry) {
    namespace MP = Physics::Multiphase;
    std::vector<ResultWriter::ScalarField> scalars;
    if (!model.enabled() || !model.initialized()) return scalars;
    auto* modelPtr = &model;
    if (registry) {
        for (const auto& variable : registry->variables()) {
            auto* scalar = variable.value;
            const std::string name = variable.descriptor.name;
            const bool temperature = variable.descriptor.unit == "K";
            scalars.push_back({
                name,
                [scalar](int i, int j, int k) {
                    return (*scalar)(i, j, k);
                },
                temperature});
        }
    }
    if (model.isThermal()) {
        if (!registry && model.hasTemperature()) {
            const std::string name = model.temperature().name().empty()
                ? "T" : model.temperature().name();
            scalars.push_back({name, [modelPtr](int i, int j, int k) {
                return modelPtr->temperature()(i, j, k);
            }, true});
        }
        return scalars;
    }
    if (model.isMixture()) {
        const std::string alphaName = model.alpha().name().empty()
            ? "alpha" : model.alpha().name();
        scalars.push_back({alphaName, [modelPtr](int i, int j, int k) {
            return modelPtr->alpha()(i, j, k);
        }});
        if (!registry && model.hasTemperature()) {
            const std::string name = model.temperature().name().empty()
                ? "T" : model.temperature().name();
            scalars.push_back({name, [modelPtr](int i, int j, int k) {
                return modelPtr->temperature()(i, j, k);
            }, true});
        }
        scalars.push_back({"PhaseDensity", [modelPtr](int i, int j, int k) {
            return modelPtr->hasMixtureProperties()
                ? modelPtr->mixtureDensity(i, j, k) : 0.0;
        }});
        scalars.push_back({"PhaseViscosity", [modelPtr](int i, int j, int k) {
            return modelPtr->hasMixtureProperties()
                ? modelPtr->mixtureViscosity(i, j, k) : 0.0;
        }});
        const bool phaseChangeConfigured =
            model.config().phaseChange.enabled
            && MP::normalizeModelType(model.config().phaseChange.model) != "none";
        if (phaseChangeConfigured || model.hasPhaseChangeSource()) {
            scalars.push_back({"PhaseChangeEnergySource",
                [modelPtr](int i, int j, int k) {
                    return modelPtr->hasPhaseChangeSource()
                        ? modelPtr->phaseChangeEnergySource(i, j, k) : 0.0;
                }});
        }
    }
    scalars.push_back({"PhaseIndicatorIntegral", [modelPtr](int, int, int) {
        return modelPtr->conservation().currentPhaseIndicator;
    }});
    scalars.push_back({"PhaseIndicatorRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().phaseIndicatorRelativeError;
    }});
    scalars.push_back({"PhaseMassIntegral", [modelPtr](int, int, int) {
        return modelPtr->conservation().currentPhaseMass;
    }});
    scalars.push_back({"PhaseMassRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().phaseMassRelativeError;
    }});
    scalars.push_back({"FlowMassRelativeError", [modelPtr](int, int, int) {
        return modelPtr->conservation().flowMassRelativeError;
    }});
    return scalars;
}

std::vector<ResultWriter::ScalarField> eulerianVTKScalars(
        Physics::PhaseSystems::PhaseSystem& system) {
    std::vector<ResultWriter::ScalarField> scalars;
    auto* pressure = &system.sharedPressure();
    scalars.push_back({"p", [pressure](int i,int j,int k) {
        return (*pressure)(i,j,k);
    }});
    for (auto& phase : system.phases()) {
        auto* state = &phase;
        const size_t phaseIndex = (size_t)(&phase - system.phases().data());
        auto* sources = &system.sources();
        const std::string& name = phase.name;
        scalars.push_back({"alpha."+name,[state](int i,int j,int k){return state->primitive.alpha(i,j,k);}});
        scalars.push_back({"rho."+name,[state](int i,int j,int k){return state->primitive.density(i,j,k);}});
        scalars.push_back({"T."+name,[state](int i,int j,int k){return state->primitive.temperature(i,j,k);},true});
        scalars.push_back({"Ux."+name,[state](int i,int j,int k){return state->primitive.velocity[0](i,j,k);}});
        scalars.push_back({"Uy."+name,[state](int i,int j,int k){return state->primitive.velocity[1](i,j,k);}});
        scalars.push_back({"Uz."+name,[state](int i,int j,int k){return state->primitive.velocity[2](i,j,k);}});
        scalars.push_back({"phaseMass."+name,[state](int i,int j,int k){return state->primary.phaseMass(i,j,k);}});
        scalars.push_back({"phaseEnthalpy."+name,[state](int i,int j,int k){return state->primary.phaseEnthalpy(i,j,k);}});
        scalars.push_back({"massSource."+name,[sources,phaseIndex](int i,int j,int k){return sources->mass[phaseIndex](i,j,k);}});
        scalars.push_back({"momentumSourceX."+name,[sources,phaseIndex](int i,int j,int k){return sources->momentum[phaseIndex][0](i,j,k);}});
        scalars.push_back({"momentumSourceY."+name,[sources,phaseIndex](int i,int j,int k){return sources->momentum[phaseIndex][1](i,j,k);}});
        scalars.push_back({"momentumSourceZ."+name,[sources,phaseIndex](int i,int j,int k){return sources->momentum[phaseIndex][2](i,j,k);}});
        scalars.push_back({"energySource."+name,[sources,phaseIndex](int i,int j,int k){return sources->energy[phaseIndex](i,j,k);}});
    }
    auto* sources = &system.sources();
    scalars.push_back({"wallBoilingMdot",[sources](int i,int j,int k){return sources->wallBoilingMass(i,j,k);}});
    scalars.push_back({"wallBoilingQconv",[sources](int i,int j,int k){return sources->wallBoilingConvectiveHeat(i,j,k);}});
    scalars.push_back({"wallBoilingQquench",[sources](int i,int j,int k){return sources->wallBoilingQuenchingHeat(i,j,k);}});
    scalars.push_back({"wallBoilingQevap",[sources](int i,int j,int k){return sources->wallBoilingEvaporativeHeat(i,j,k);}});
    scalars.push_back({"wallBoilingDepartureDiameter",[sources](int i,int j,int k){return sources->wallBoilingDepartureDiameter(i,j,k);}});
    scalars.push_back({"interphaseMechanicalHeating",[sources](int i,int j,int k){return sources->interphaseMechanicalHeating(i,j,k);}});
    return scalars;
}

std::vector<ResultWriter::ScalarField> eulerianTurbulenceVTKScalars(
        const Turbulence::EquationSystem& turbulence) {
    std::vector<ResultWriter::ScalarField> scalars;
    for (const auto& state : turbulence.states()) {
        const auto* statePtr = &state;
        for (const auto* equation : turbulence.transportEquations(state)) {
            const auto* equationPtr = equation;
            scalars.push_back({equation->variable.name(),
                [equationPtr](int i,int j,int k) {
                    return equationPtr->variable(i,j,k);
                }});
        }
        scalars.push_back({state.eddyViscosity.name(),
            [statePtr](int i,int j,int k) {
                return statePtr->eddyViscosity(i,j,k);
            }});
    }
    return scalars;
}

} // namespace SF::Application::Output
