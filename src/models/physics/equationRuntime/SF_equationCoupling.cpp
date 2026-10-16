/// @file SF_equationCoupling.cpp
/// @brief 把 multiphase/model state 适配为 equation assembly contribution。
///
/// Data flow:
///   phase/model state + residual workspace
///       -> source/closure evaluation and registered-variable updates
///       -> equation-system assembly/commit contract
///
/// 本文件不选择 runner，不推进 global clock，也不实现 MPI backend。

#include "models/physics/equationRuntime/SF_mixtureEquationProvider.h"

#include "SF_MultiBlockMesh.h"
#include "SF_homogeneousPhaseChange.h"
#include "SF_multiphase.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Equation::Coupling {

void registerMixtureState(
        Physics::Multiphase::MultiPhaseModel& model,
        const Field& field,
        std::vector<double>& phaseMassRHS,
        State::VariableRegistry& registry) {
    auto descriptor = [&](const std::string& name,
                          State::ConservationKind conservation,
                          State::UpdatePolicy update,
                          const std::string& unit) {
        State::VariableDescriptor result;
        result.name = name;
        result.conservation = conservation;
        result.updatePolicy = update;
        result.unit = unit;
        return result;
    };
    if (model.isMixture()) {
        phaseMassRHS.assign((size_t)field.TotalSize(), 0.0);
        registry.add(
            descriptor(model.phaseMass().name(),
                       State::ConservationKind::Conservative,
                       State::UpdatePolicy::Explicit, "kg/m3"),
            model.phaseMass(), phaseMassRHS);
        if (model.hasTemperature()) {
            registry.add(
                descriptor(model.temperature().name(),
                           State::ConservationKind::Algorithmic,
                           State::UpdatePolicy::DerivedOnly, "K"),
                model.temperature());
        }
        return;
    }
    if (model.isThermal() && model.hasTemperature()) {
        registry.add(
            descriptor(model.temperature().name(),
                       State::ConservationKind::Algorithmic,
                       State::UpdatePolicy::ExternalExplicit, "K"),
            model.temperature());
    }
}

void CompositeTransportProvider::setTurbulence(FDM::ITransportModel* model) {
    turbulence_ = model;
}

void CompositeTransportProvider::setMultiPhase(
        Physics::Multiphase::MultiPhaseModel* model) {
    multiPhase_ = model;
}

void CompositeTransportProvider::setInterfaceModel(FDM::ITransportModel* model) {
    interfaceModel_ = model;
}

bool CompositeTransportProvider::active() const {
    return turbulence_ != nullptr || multiPhase_ != nullptr
        || interfaceModel_ != nullptr;
}

void CompositeTransportProvider::applyBoundary(const Field& field) {
    if (turbulence_) turbulence_->applyBoundary(field);
}

void CompositeTransportProvider::refreshClosure(const Field& field) {
    if (turbulence_) turbulence_->refreshClosure(field);
}
FDM::OnceStepTransportData CompositeTransportProvider::onceStepTransport() {
    return turbulence_ ? turbulence_->onceStepTransport() : FDM::OnceStepTransportData{};
}

double CompositeTransportProvider::dynamicViscosity(
        const Field& field, int i, int j, int k,
        double laminarMu) const {
    double mu = laminarMu;
    if (multiPhase_) {
        mu = multiPhase_->dynamicViscosity(field, i, j, k, mu);
    }
    if (interfaceModel_) {
        mu = interfaceModel_->dynamicViscosity(field, i, j, k, mu);
    }
    if (turbulence_) {
        mu = turbulence_->dynamicViscosity(field, i, j, k, mu);
    }
    return mu;
}

std::vector<std::string> CompositeTransportProvider::stageScalarTargets() const {
    // The current stage coefficient port has one molecular transport binding.
    // Mixture/interface combinations require their own resolved coefficient port.
    if(!turbulence_ || multiPhase_ || interfaceModel_)return {};
    return turbulence_->stageScalarTargets();
}
FDM::ITransportModel& CompositeTransportProvider::stageTransport() const {
    if(stageScalarTargets().empty())throw std::runtime_error("Unsupported composite stage coefficient port: requires one active transport model, no mixture/interface coefficient composition.");
    return *turbulence_;
}
FDM::TransportStageCoefficients CompositeTransportProvider::stageCoefficients(const Field& field,
        const FDM::TransportStageSamples& samples,int i,int j,int k) const {
    return stageTransport().stageCoefficients(field,samples,i,j,k);
}
FDM::ScalarEquationData CompositeTransportProvider::stageScalarData(const Field& field,const std::string& target) const {
    return stageTransport().stageScalarData(field,target);
}
void CompositeTransportProvider::publishStageViscosity(const Field& field,const std::vector<double>& values) {
    stageTransport().publishStageViscosity(field,values);
}

namespace {

void appendUnique(std::vector<std::string>& destination,
                  const std::vector<std::string>& source) {
    for (const auto& name : source) {
        bool duplicate = false;
        for (const auto& current : destination) duplicate |= current == name;
        if (!duplicate) destination.push_back(name);
    }
}

} // namespace

std::vector<std::string>
CompositeTransportProvider::distributedReadFields() const {
    std::vector<std::string> result;
    if (turbulence_) appendUnique(
        result, turbulence_->distributedReadFields());
    if (interfaceModel_) appendUnique(
        result, interfaceModel_->distributedReadFields());
    return result;
}

std::vector<std::string>
CompositeTransportProvider::distributedWriteFields() const {
    std::vector<std::string> result;
    if (turbulence_) appendUnique(
        result, turbulence_->distributedWriteFields());
    if (interfaceModel_) appendUnique(
        result, interfaceModel_->distributedWriteFields());
    return result;
}

int CompositeTransportProvider::distributedHaloDepth() const {
    int depth = 0;
    if (turbulence_) depth = std::max(
        depth, turbulence_->distributedHaloDepth());
    if (interfaceModel_) depth = std::max(
        depth, interfaceModel_->distributedHaloDepth());
    return depth;
}

HomogeneousPhaseChangeProvider::HomogeneousPhaseChangeProvider(
        Physics::FluidStateModel::HomogeneousPhaseChange& model)
    : model_(model) {}

void HomogeneousPhaseChangeProvider::assembleRHS(
        Field& field, Residual& residual, double dt) {
    model_.assemble(field, residual, dt);
}

void HomogeneousPhaseChangeProvider::commitStep(Field&, double) {
    Physics::PhaseChange::reportDiagnostics(model_.diagnostics());
}

MixtureEquationProvider::MixtureEquationProvider(
        Physics::Multiphase::MultiPhaseModel& model,
        std::vector<double>& rhs,
        State::VariableRegistry& variables,
        FDM::IExecutionRuntime& runtime)
    : model_(model), rhs_(rhs), variables_(variables), runtime_(runtime) {}

void MixtureEquationProvider::prepareRHS(Field& field, double) {
    model_.applyAuxiliaryBoundaryConditions(field);
    std::vector<Execution::FieldAccess> reads;
    for (const auto& variable : variables_.variables()) {
        reads.push_back(Execution::readHalo(variable.descriptor.name, 1));
    }
    if (!reads.empty()) {
        runtime_.prepare({"multiphase mixture auxiliary stencil", reads});
    }
    model_.applyAuxiliaryBoundaryConditions(field);
    model_.updateFlowCoupling(field);
}

void MixtureEquationProvider::assembleRHS(
        Field& field, Residual& residual, double dt) {
    if (model_.isMixture()) {
        model_.assembleIntegratedRHS(field, residual, rhs_, dt);
    } else {
        model_.addSourceTerms(field, residual);
    }
}

void MixtureEquationProvider::commitStep(Field& field, double) {
    model_.commitTimeLevel(field);
    std::vector<Execution::FieldAccess> writes;
    for (const auto& variable : variables_.variables()) {
        writes.push_back(Execution::writeOwned(variable.descriptor.name));
    }
    if (!writes.empty()) {
        runtime_.finalize({"multiphase mixture committed state", writes});
        std::vector<Execution::FieldAccess> reads;
        for (const auto& variable : variables_.variables()) {
            reads.push_back(Execution::readHalo(variable.descriptor.name, 1));
        }
        runtime_.prepare({"multiphase mixture committed halo", reads});
    }
    model_.applyAuxiliaryBoundaryConditions(field);
    model_.updateFlowCoupling(field);
}

State::VariableRegistry* MixtureEquationProvider::variables(Field&) {
    return variables_.empty() ? nullptr : &variables_;
}

MultiPatchMixtureEquationProvider::MultiPatchMixtureEquationProvider(
        MultiBlockMesh& mesh,
        const std::vector<int>& localPatchIds,
        FDM::IExecutionRuntime& runtime,
        std::vector<Physics::Multiphase::MultiPhaseModel>& models,
        std::vector<std::vector<double>>& rhs,
        std::vector<State::VariableRegistry>& variables,
        std::vector<CompositeTransportProvider>& transports)
    : mesh_(mesh), localPatchIds_(localPatchIds), runtime_(runtime)
    , models_(models), rhs_(rhs), variables_(variables)
    , transports_(transports) {}

void MultiPatchMixtureEquationProvider::prepareRHS(
        const std::vector<Field*>&, double) {
    for (int patchId : localPatchIds_) {
        models_.at((size_t)patchId).applyAuxiliaryBoundaryConditions(
            mesh_.block((size_t)patchId).field);
    }
    const auto halo = auxiliaryHaloContract(
        "multi-patch mixture auxiliary stencil");
    if (!halo.accesses.empty()) runtime_.prepare(halo);
    for (int patchId : localPatchIds_) {
        auto& model = models_.at((size_t)patchId);
        auto& field = mesh_.block((size_t)patchId).field;
        model.applyAuxiliaryBoundaryConditions(field);
        model.updateFlowCoupling(field);
    }
}

void MultiPatchMixtureEquationProvider::assembleRHS(
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals, double dt) {
    if (fields.size() != residuals.size()) {
        throw std::runtime_error("multi-patch source workspace size mismatch.");
    }
    for (int patchId : localPatchIds_) {
        auto& model = models_.at((size_t)patchId);
        auto& field = mesh_.block((size_t)patchId).field;
        Residual* residual = nullptr;
        for (size_t index = 0; index < fields.size(); ++index) {
            if (fields[index] == &field) { residual = residuals[index]; break; }
        }
        if (!residual) throw std::runtime_error("missing residual for local patch.");
        if (model.isMixture()) {
            model.assembleIntegratedRHS(
                field, *residual, rhs_.at((size_t)patchId), dt);
        } else {
            model.addSourceTerms(field, *residual);
        }
    }
}

void MultiPatchMixtureEquationProvider::commitStep(
        const std::vector<Field*>&, double) {
    for (int patchId : localPatchIds_) {
        models_.at((size_t)patchId).commitTimeLevel(
            mesh_.block((size_t)patchId).field);
    }
    std::vector<Execution::FieldAccess> writes;
    for (const auto& access : auxiliaryHaloContract("unused").accesses) {
        writes.push_back(Execution::writeOwned(access.field));
    }
    if (!writes.empty()) {
        runtime_.finalize({"multi-patch mixture committed state", writes});
        runtime_.prepare(auxiliaryHaloContract(
            "multi-patch mixture committed halo"));
    }
    for (int patchId : localPatchIds_) {
        auto& model = models_.at((size_t)patchId);
        auto& field = mesh_.block((size_t)patchId).field;
        model.applyAuxiliaryBoundaryConditions(field);
        model.updateFlowCoupling(field);
    }
}

State::VariableRegistry* MultiPatchMixtureEquationProvider::variables(
        Field& field) {
    auto& registry = variables_.at(findPatch(field));
    return registry.empty() ? nullptr : &registry;
}

const FDM::ITransportModel* MultiPatchMixtureEquationProvider::transportModel(
        const Field& field) const {
    return &transports_.at(findPatch(field));
}

Execution::OperatorContract
MultiPatchMixtureEquationProvider::auxiliaryHaloContract(
        const char* operation) const {
    Execution::OperatorContract contract;
    contract.name = operation;
    for (const auto& registry : variables_) {
        for (const auto& variable : registry.variables()) {
            bool duplicate = false;
            for (const auto& access : contract.accesses) {
                duplicate |= access.field == variable.descriptor.name;
            }
            if (!duplicate) {
                // 当前 legacy 标量对流为一阶迎风；模板需求由这里声明，
                // 不写入 VariableDescriptor 或 Field 注册信息。
                contract.accesses.push_back(Execution::readHalo(
                    variable.descriptor.name, 1));
            }
        }
    }
    return contract;
}

size_t MultiPatchMixtureEquationProvider::findPatch(
        const Field& field) const {
    for (size_t patchId = 0; patchId < mesh_.size(); ++patchId) {
        if (&mesh_.block(patchId).field == &field) return patchId;
    }
    throw std::runtime_error(
        "Multi-patch equation system received an unregistered Field.");
}

} // namespace SF::Equation::Coupling
