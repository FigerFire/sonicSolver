#include "SF_advection.h"
/// @file SF_levelSetOperations.cpp
/// @brief 实现 single/multi-patch interface equation coupling 与同步声明。
///
/// Data flow:
///   Field + level-set state + ExecutionRuntime interface
///       -> boundary/halo-ready interface state
///       -> RHS、jump condition 与 committed interface state
///
/// 本文件保留既有数值顺序，不创建 communicator 或 global timestep loop。

#include "solver/algorithm/interface/SF_levelSetOperations.h"

#include "SF_MultiBlockMesh.h"
#include "SF_interfaceModel.h"
#include "SF_hjWeno.h"
#include "SF_levelSet.h"
#include "SF_publicationTrace.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::SolverAlgorithm {

static void registerInterfaceState(
        Physics::InterfaceModels::Model& model,
        const Field& field,
        State::VariableRegistry& registry,std::vector<double>& rhs) {
    (void)field;
    State::VariableDescriptor descriptor;
    descriptor.name = "phi";
    descriptor.conservation = State::ConservationKind::AdvectedGeometry;
    descriptor.updatePolicy = State::UpdatePolicy::HamiltonJacobi;
    descriptor.unit = "m";
    registry.add(std::move(descriptor), model.primaryScalar(),
                 rhs);
}

State::DistributedFieldView interfaceNormalView(Physics::InterfaceModels::Model& model,
        Field& field,int patchId) {
    auto* state=model.levelSetState();
    if (!state) throw std::runtime_error("Interface normal view requires level-set storage.");
    State::DistributedFieldView view;
    view.name="levelSetNormal";view.geometry=&field;view.blockId=patchId;
    view.components=3;view.haloDepth=field.NG();view.exchange=State::ExchangeKind::State;
    view.read=[state](int cell,int component) {return state->normals().at((size_t)cell)[component];};
    view.write=[state](int cell,int component,double value) {state->normals().at((size_t)cell)[component]=value;};
    return view;
}

InterfaceEquationProvider::InterfaceEquationProvider(
        Physics::InterfaceModels::Model& model,
        State::VariableRegistry& variables,
        FDM::IExecutionRuntime& runtime,Field& field)
    : model_(model), variables_(variables), runtime_(runtime) {
    rhs_.assign(static_cast<size_t>(field.TotalSize()),0.0);
    registerInterfaceState(model,field,variables,rhs_);
}

void InterfaceEquationProvider::beginStep(Field& field, double dt) {
    field_=&field;
    if(!std::isfinite(dt) || dt<=0)throw std::runtime_error("Invalid interface physical dt");
    if(model_.config().levelSet.massCorrection)throw std::runtime_error("Unsupported global level-set mass correction");
}

void InterfaceEquationProvider::prepareInterfaceState(Field& field) {
    model_.applyBoundary(field);
    if (model_.levelSetState()->takeBoundaryModification())
        runtime_.finalize({"level-set physical boundary write",{Execution::writeOwned("phi")}});
    const int depth = Physics::Multiphase::HJWeno::requiredGhostLayers(
        model_.config().levelSet.advectionOrder);
    runtime_.prepare({"level-set advection stencil",{Execution::readHalo("phi",depth)}});
    Physics::Multiphase::tracePublication("C2-phi-halo",field,*model_.levelSetState());
    model_.applyBoundary(field);
    if (model_.levelSetState()->takeBoundaryModification())
        runtime_.finalize({"level-set boundary after halo",{Execution::writeOwned("phi")}});
    runtime_.prepare({"level-set boundary-ready phi",{Execution::readHalo("phi",depth)}});
    model_.computeNormals(field);
    runtime_.finalize({"level-set normal update",{Execution::writeOwned("levelSetNormal")}});
    runtime_.prepare({"level-set normal stencil",{Execution::readHalo("levelSetNormal",1)}});
    Physics::Multiphase::tracePublication("C6-published-normal",field,*model_.levelSetState());
    model_.levelSetState()->publishNormalStencil();
    model_.refreshFromNormals(field);
    runtime_.finalize({"level-set geometry update",{Execution::writeOwned("levelSetCurvature")}});
    runtime_.prepare({"level-set curvature stencil",{Execution::readHalo("levelSetCurvature",1)}});
    model_.levelSetState()->publishCurvatureStencil();
    // A single local patch may still be one partition of a distributed domain.
    const auto diagnostics=model_.conservation();
    const double indicator=runtime_.globalSum(diagnostics.currentPhaseIndicator);
    const double mass=runtime_.globalSum(diagnostics.currentFlowMass);
    model_.publishConservation(indicator,mass);
    Physics::Multiphase::tracePublication("C12-published-geometry",field,*model_.levelSetState());
}

void InterfaceEquationProvider::prepareRHS(Field& field, double) {
    prepareInterfaceState(field);
    const auto& config=model_.config().levelSet;
    AdvectOptions options{config.advectionOrder,config.wenoEpsilon,config.wenoPower,config.advectFluidCellsOnly,true};
    Advect::assembleRHS(field,*model_.levelSetState(),options,rhs_);
    Physics::Multiphase::tracePublication("C4-advection-RHS",field,*model_.levelSetState(),&rhs_);
}

void InterfaceEquationProvider::assembleRHS(
        Field& field, Residual& residual, double) {
    model_.addSourceTerms(field, residual);
}

void InterfaceEquationProvider::preparePressureCorrection(Field& field) {
    prepareInterfaceState(field);
}

void InterfaceEquationProvider::validateAdvection(int order,double epsilon,double power,bool csf,bool ghostFluid,double sigma,double width) const {
    const auto& selected=model_.config().levelSet;
    if (selected.advectionOrder!=order || selected.wenoEpsilon!=epsilon || selected.wenoPower!=power
        || (Physics::Multiphase::normalizeModelType(selected.surfaceTensionModel)=="csf")!=csf
        || (Physics::Multiphase::normalizeModelType(selected.surfaceTensionModel)=="ghostfluid")!=ghostFluid
        || selected.surfaceTension!=sigma || selected.interfaceThickness!=width)
        throw std::runtime_error("Level-set bound kernel differs from frozen advection WHICH parameters.");
}
void InterfaceEquationProvider::beginReinitialization(int order,double pseudoDt,
        double epsilon,double power,double signFactor) {
    if (!field_ || !model_.levelSetState())
        throw std::runtime_error("Level-set pseudo-time requires a bound active patch.");
    pseudoOptions_={1,pseudoDt,order,epsilon,power,signFactor,[this]() {
        model_.applyBoundary(*field_);
        runtime_.finalize({"level-set reinitialization stage",{Execution::writeOwned("phi")}});
        runtime_.prepare({"level-set reinitialization stencil",
            {Execution::readHalo("phi",Physics::Multiphase::HJWeno::requiredGhostLayers(pseudoOptions_.order))}});
        model_.applyBoundary(*field_);
    }};
    Reinit::begin(*field_,*model_.levelSetState(),pseudoOptions_,pseudoWorkspace_);
}
State::DistributedFieldView InterfaceEquationProvider::referenceView() {
    if (!field_ || pseudoWorkspace_.reference.size()!=(size_t)field_->TotalSize())
        throw std::runtime_error("Level-set reference view requires a frozen pseudo-time snapshot.");
    State::DistributedFieldView view;view.name="levelSet.phi0";view.geometry=field_;
    view.components=1;view.haloDepth=field_->NG();view.exchange=State::ExchangeKind::None;
    view.read=[this](int cell,int) {return pseudoWorkspace_.reference.at((size_t)cell);};
    view.write=[](int,int,double) {throw std::runtime_error("Frozen phi0 view is read-only during pseudo-time.");};
    return view;
}
void InterfaceEquationProvider::reinitializeStage() {
    if (!field_) throw std::runtime_error("Level-set pseudo stage before physical step begin.");
    Reinit::stage(*field_,*model_.levelSetState(),pseudoOptions_,pseudoWorkspace_);
}
void InterfaceEquationProvider::publishGeometry(double dt) {
    if (!field_) throw std::runtime_error("Level-set geometry publication before physical step begin.");
    if (!std::isfinite(dt) || dt<=0) throw std::runtime_error("Invalid level-set physical commit dt.");
    prepareInterfaceState(*field_);
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
    , models_(models), variables_(variables) {
    rhs_.resize(mesh.size());
    for(size_t id=0;id<mesh.size();++id) {
        auto& field=mesh.block(id).field;rhs_[id].assign(static_cast<size_t>(field.TotalSize()),0.0);
        registerInterfaceState(*models[id],field,variables[id],rhs_[id]);
    }
}

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
        if(!std::isfinite(dt) || dt<=0)throw std::runtime_error("Invalid multi-patch interface dt");
    }
}

void MultiPatchInterfaceEquationProvider::prepareRHS(
        const std::vector<Field*>&, double) {
    prepareInterfaceState("multi-patch level-set advection stencil");
    for (int patchId : localPatchIds_) {
        auto& model=*models_.at(static_cast<size_t>(patchId));auto& field=mesh_.block(static_cast<size_t>(patchId)).field;
        const auto& config=model.config().levelSet;
        AdvectOptions options{config.advectionOrder,config.wenoEpsilon,config.wenoPower,config.advectFluidCellsOnly,true};
        Advect::assembleRHS(field,*model.levelSetState(),options,rhs_[patchId]);
        Physics::Multiphase::tracePublication("C4-advection-RHS",field,*model.levelSetState(),&rhs_[patchId]);
    }
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
    prepareInterfaceState("multi-patch pressure-interface phi stencil");
}

void MultiPatchInterfaceEquationProvider::publishGeometry(double dt) {
    if (!std::isfinite(dt) || dt<=0) throw std::runtime_error("Invalid multi-patch level-set commit dt.");
    prepareInterfaceState("multi-patch committed level-set stencil");
}

void MultiPatchInterfaceEquationProvider::validateAdvection(int order,double epsilon,double power,bool csf,bool ghostFluid,double sigma,double width) const {
    for (int id:localPatchIds_) {
        const auto& selected=models_.at((size_t)id)->config().levelSet;
        if (selected.advectionOrder!=order || selected.wenoEpsilon!=epsilon || selected.wenoPower!=power
        || (Physics::Multiphase::normalizeModelType(selected.surfaceTensionModel)=="csf")!=csf
        || (Physics::Multiphase::normalizeModelType(selected.surfaceTensionModel)=="ghostfluid")!=ghostFluid
        || selected.surfaceTension!=sigma || selected.interfaceThickness!=width)
            throw std::runtime_error("Multi-patch level-set kernel differs from frozen WHICH parameters.");
    }
}
State::DistributedFieldView MultiPatchInterfaceEquationProvider::referenceView() {
    throw std::runtime_error("Unsupported: multi-patch level-set pseudo-time snapshot has no coupled scheduler.");
}
void MultiPatchInterfaceEquationProvider::beginReinitialization(int,double,double,double,double) {
    throw std::runtime_error("Unsupported: multi-patch level-set requires a coupled pseudo-time scheduler.");
}
void MultiPatchInterfaceEquationProvider::reinitializeStage() {
    throw std::runtime_error("Unsupported: multi-patch level-set pseudo-time stage.");
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

void MultiPatchInterfaceEquationProvider::prepareInterfaceState(const char* operation) {
    const auto boundary=[&] {
        bool modified=false;
        for (int id:localPatchIds_) {
            auto& model=*models_.at((size_t)id);
            model.applyBoundary(mesh_.block((size_t)id).field);
            modified=model.levelSetState()->takeBoundaryModification() || modified;
        }
        if (modified) runtime_.finalize({"multi-patch level-set physical boundary write",
                                        {Execution::writeOwned("phi")}});
    };
    boundary();
    preparePhiHalo(operation);
    for (int id:localPatchIds_)
        Physics::Multiphase::tracePublication("C2-phi-halo",mesh_.block((size_t)id).field,*models_[id]->levelSetState());
    boundary();
    // The second physical BC can consume a newly published tangential donor.
    // The runtime reuses the halo unless an owner actually changed.
    preparePhiHalo("multi-patch boundary-ready level-set stencil");
    for (int id:localPatchIds_) models_[id]->computeNormals(mesh_.block((size_t)id).field);
    runtime_.finalize({"multi-patch level-set normal update",{Execution::writeOwned("levelSetNormal")}});
    runtime_.prepare({"multi-patch normal stencil",{Execution::readHalo("levelSetNormal",1)}});
    for (int id:localPatchIds_) {
        auto& field=mesh_.block((size_t)id).field;
        auto& model=*models_[id];
        Physics::Multiphase::tracePublication("C6-published-normal",field,*model.levelSetState());
        model.levelSetState()->publishNormalStencil();
        model.refreshFromNormals(field);
    }
    publishCurvature("multi-patch curvature stencil");
    for (int id:localPatchIds_) models_[id]->levelSetState()->publishCurvatureStencil();
    double indicator=0.0,mass=0.0;
    for (int id:localPatchIds_) {
        indicator+=models_[id]->conservation().currentPhaseIndicator;
        mass+=models_[id]->conservation().currentFlowMass;
    }
    indicator=runtime_.globalSum(indicator);
    mass=runtime_.globalSum(mass);
    for (int id:localPatchIds_) models_[id]->publishConservation(indicator,mass);
    for (int id:localPatchIds_)
        Physics::Multiphase::tracePublication("C12-published-geometry",mesh_.block((size_t)id).field,*models_[id]->levelSetState());
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

} // namespace SF::SolverAlgorithm
