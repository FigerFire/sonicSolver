/// @file SF_interfaceCoupling.cpp
/// @brief 实现 single/multi-patch interface equation coupling 与同步声明。
///
/// Data flow:
///   Field + level-set state + ExecutionRuntime interface
///       -> boundary/halo-ready interface state
///       -> RHS、jump condition 与 committed interface state
///
/// 本文件保留既有数值顺序，不创建 communicator 或 global timestep loop。

#include "solver/equation/coupling/SF_interfaceCoupling.h"

#include "SF_MultiBlockMesh.h"
#include "SF_interfaceModel.h"
#include "SF_hjWeno.h"
#include "SF_levelSet.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::Equation::Coupling {

void registerInterfaceState(
        Physics::InterfaceModels::Model& model,
        const Field& field,
        State::VariableRegistry& registry) {
    (void)field;
    State::VariableDescriptor descriptor;
    descriptor.name = "phi";
    descriptor.conservation = State::ConservationKind::AdvectedGeometry;
    descriptor.updatePolicy = State::UpdatePolicy::HamiltonJacobi;
    descriptor.unit = "m";
    registry.add(std::move(descriptor), model.primaryScalar(),
                 model.primaryScalarRHS());
}

InterfaceEquationProvider::InterfaceEquationProvider(
        Physics::InterfaceModels::Model& model,
        State::VariableRegistry& variables,
        FDM::IExecutionRuntime& runtime)
    : model_(model), variables_(variables), runtime_(runtime) {}

void InterfaceEquationProvider::beginStep(Field& field, double dt) {
    (void)field;
    model_.beginTimeStep(dt);
}

void InterfaceEquationProvider::prepareInterfaceState(Field& field) {
    model_.applyBoundary(field);
    const int depth = Physics::Multiphase::HJWeno::requiredGhostLayers(
        model_.config().levelSet.advectionOrder);
    runtime_.prepare({"level-set advection stencil",
                      {Execution::readHalo("phi", depth)}});
    model_.applyBoundary(field);
    model_.refresh(field);
    runtime_.finalize({"level-set geometry update",
                       {Execution::writeOwned("levelSetCurvature")}});
    runtime_.prepare({"level-set curvature stencil",
                      {Execution::readHalo("levelSetCurvature", 1)}});
}

void InterfaceEquationProvider::prepareRHS(Field& field, double) {
    prepareInterfaceState(field);
    model_.assembleTransportRHS(field);
}

void InterfaceEquationProvider::assembleRHS(
        Field& field, Residual& residual, double) {
    model_.addSourceTerms(field, residual);
}

void InterfaceEquationProvider::preparePressureCorrection(Field& field) {
    prepareInterfaceState(field);
}

void InterfaceEquationProvider::commitStep(Field& field, double dt) {
    if (auto* levelSet = model_.levelSetState()) {
        const auto& config = model_.config();
        if (config.levelSet.reinitializationSteps > 0) {
            levelSet->reinitialize(
                field, config.levelSet.reinitializationSteps,
                config.levelSet.pseudoTimeStep,
                config.levelSet.reinitializationOrder,
                config.levelSet.wenoEpsilon,
                config.levelSet.wenoPower,
                config.levelSet.signSmoothingFactor,
                [&]() {
                    model_.applyBoundary(field);
                    runtime_.finalize({"level-set reinitialization stage",
                                       {Execution::writeOwned("phi")}});
                    const int depth =
                        Physics::Multiphase::HJWeno::requiredGhostLayers(
                            config.levelSet.reinitializationOrder);
                    runtime_.prepare({"level-set reinitialization stencil",
                                      {Execution::readHalo("phi", depth)}});
                    model_.applyBoundary(field);
                });
        }
    }
    model_.completeTimeStep(field, dt);
    prepareInterfaceState(field);
}

State::VariableRegistry* InterfaceEquationProvider::variables(Field&) {
    return variables_.empty() ? nullptr : &variables_;
}

const FDM::ITransportModel* InterfaceEquationProvider::transportModel(
        const Field&) const {
    return &model_;
}

const FDM::IInterfaceJumpCondition*
InterfaceEquationProvider::interfaceJumpCondition(const Field&) const {
    return &model_;
}

MultiPatchInterfaceEquationProvider::MultiPatchInterfaceEquationProvider(
        MultiBlockMesh& mesh,
        const std::vector<int>& localPatchIds,
        FDM::IExecutionRuntime& runtime,
        std::vector<std::unique_ptr<Physics::InterfaceModels::Model>>& models,
        std::vector<State::VariableRegistry>& variables)
    : mesh_(mesh), localPatchIds_(localPatchIds), runtime_(runtime)
    , models_(models), variables_(variables) {}

void MultiPatchInterfaceEquationProvider::beginStep(
        const std::vector<Field*>&, double dt) {
    for (int patchId : localPatchIds_) {
        auto& model = *models_.at((size_t)patchId);
        if (mesh_.size() > 1
            && model.config().levelSet.reinitializationSteps > 0) {
            throw std::runtime_error(
                "Multi-patch Level Set reinitialization requires a coupled "
                "pseudo-time stage scheduler. Set reinitializationSteps 0; "
                "the solver will not perform unsynchronized reinitialization.");
        }
        model.beginTimeStep(dt);
    }
}

void MultiPatchInterfaceEquationProvider::prepareRHS(
        const std::vector<Field*>&, double) {
    for (int patchId : localPatchIds_) {
        auto& field = mesh_.block((size_t)patchId).field;
        models_.at((size_t)patchId)->applyBoundary(field);
    }
    preparePhiHalo("multi-patch level-set advection stencil");
    for (int patchId : localPatchIds_) {
        auto& field = mesh_.block((size_t)patchId).field;
        auto& model = *models_.at((size_t)patchId);
        model.applyBoundary(field);
        model.refresh(field);
        model.assembleTransportRHS(field);
    }
    publishCurvature("multi-patch level-set curvature stencil");
}

void MultiPatchInterfaceEquationProvider::assembleRHS(
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals, double) {
    if (fields.size() != residuals.size()) {
        throw std::runtime_error("multi-patch interface workspace size mismatch.");
    }
    for (int patchId : localPatchIds_) {
        auto& field = mesh_.block((size_t)patchId).field;
        Residual* residual = nullptr;
        for (size_t index = 0; index < fields.size(); ++index) {
            if (fields[index] == &field) { residual = residuals[index]; break; }
        }
        if (!residual) throw std::runtime_error("missing residual for local interface patch.");
        models_.at((size_t)patchId)->addSourceTerms(field, *residual);
    }
}

void MultiPatchInterfaceEquationProvider::preparePressureCorrection(
        const std::vector<Field*>&) {
    for (int patchId : localPatchIds_) {
        models_.at((size_t)patchId)->applyBoundary(
            mesh_.block((size_t)patchId).field);
    }
    preparePhiHalo("multi-patch pressure-interface phi stencil");
    for (int patchId : localPatchIds_) {
        auto& field = mesh_.block((size_t)patchId).field;
        auto& model = *models_.at((size_t)patchId);
        model.applyBoundary(field);
        model.refresh(field);
    }
    publishCurvature("multi-patch pressure-interface curvature stencil");
}

void MultiPatchInterfaceEquationProvider::commitStep(
        const std::vector<Field*>&, double dt) {
    for (int patchId : localPatchIds_) {
        models_.at((size_t)patchId)->applyBoundary(
            mesh_.block((size_t)patchId).field);
    }
    preparePhiHalo("multi-patch committed level-set stencil");
    for (int patchId : localPatchIds_) {
        auto& field = mesh_.block((size_t)patchId).field;
        auto& model = *models_.at((size_t)patchId);
        model.completeTimeStep(field, dt);
    }
}

State::VariableRegistry* MultiPatchInterfaceEquationProvider::variables(
        Field& field) {
    auto& registry = variables_.at(findPatch(field));
    return registry.empty() ? nullptr : &registry;
}

const FDM::ITransportModel*
MultiPatchInterfaceEquationProvider::transportModel(
        const Field& field) const {
    return models_.at(findPatch(field)).get();
}

const FDM::IInterfaceJumpCondition*
MultiPatchInterfaceEquationProvider::interfaceJumpCondition(
        const Field& field) const {
    return models_.at(findPatch(field)).get();
}

int MultiPatchInterfaceEquationProvider::phiHaloDepth() const {
    int depth = 0;
    for (const auto& model : models_) {
        if (!model) continue;
        depth = std::max(
            depth,
            Physics::Multiphase::HJWeno::requiredGhostLayers(
                model->config().levelSet.advectionOrder));
    }
    if (depth <= 0) {
        throw std::runtime_error(
            "Multi-patch interface coupling has no active Level Set model.");
    }
    return depth;
}

void MultiPatchInterfaceEquationProvider::preparePhiHalo(
        const char* operation) {
    runtime_.prepare({operation,
                      {Execution::readHalo("phi", phiHaloDepth())}});
}

void MultiPatchInterfaceEquationProvider::publishCurvature(
        const char* operation) {
    runtime_.finalize({"multi-patch level-set geometry update",
                       {Execution::writeOwned("levelSetCurvature")}});
    runtime_.prepare({operation,
                      {Execution::readHalo("levelSetCurvature", 1)}});
}

size_t MultiPatchInterfaceEquationProvider::findPatch(
        const Field& field) const {
    for (size_t patchId = 0; patchId < mesh_.size(); ++patchId) {
        if (&mesh_.block(patchId).field == &field) return patchId;
    }
    throw std::runtime_error(
        "Multi-patch interface coupling received an unregistered Field.");
}

} // namespace SF::Equation::Coupling
