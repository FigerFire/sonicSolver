/// @file SF_IBM.cpp
/// @brief IBM 薄 facade 的生命周期和端口转发实现。

#include "SF_IBM.h"

#include "common/SF_capabilityConfig.h"
#include "common/SF_geometrySetup.h"
#include "common/SF_methodSelection.h"
#include "common/SF_methodSetup.h"
#include "common/SF_stageTime.h"

#include <stdexcept>

namespace SF::IBM {

bool IB::setup(Field& field,
               const std::string& caseDir,
               const std::vector<std::string>& stlFiles,
               const IBMRuntimeConfig& config) {
    active_ = false;
    config_ = config;
    wallField_ = &field;
    selection_ = {};
    capabilities_ = {};
    descriptor_ = {};
    fluidPorts_.clear();
    if (config_.wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic
        && (config_.method!=FDM::IBMMethod::Ghost || config_.ilwEnabled))
        throw std::runtime_error("stationaryNoSlipAdiabatic requires Ghost and explicit Euler ILW disablement.");
    if (stlFiles.empty()) return false;

    const Common::GeometrySetup geometry =
        Common::loadGeometry(geometry_, stlFiles, caseDir);
    if (!geometry.loaded) return false;

    selection_ = Common::selectMethod(config_);
    fluidPorts_ = selection_.fluidPorts;
    capabilities_ = Common::configureCapabilities(selection_);
    descriptor_ = Common::setupMethod(
        field, geometry_, config_, geometry.loadSeconds,
        ghostWeights_, forcing_);
    active_ = true;
    return true;
}

bool IB::usesGhostCells() const {
    return active_ && config_.method == FDM::IBMMethod::Ghost;
}

bool IB::usesForcing() const {
    return active_ && config_.method == FDM::IBMMethod::VariationalForcing;
}

const FDM::ImmersedMethodSelection& IB::methodSelection() const {
    return selection_;
}

const FDM::ImmersedMethodCapabilities& IB::capabilities() const {
    return capabilities_;
}

const FDM::ImmersedAlgorithmDescriptor& IB::algorithmDescriptor() const {
    return descriptor_;
}

const std::vector<FDM::ImmersedFluidPort>& IB::fluidPorts() const {
    return fluidPorts_;
}

void IB::applyGhostCells(Field& field, double time, double dt) {
    if (!usesGhostCells()) return;
    if (config_.wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic && wallField_!=&field)
        throw std::runtime_error("Viscous immersed wall is bound to a different physical patch.");
    Common::validateStageTime(time, dt);
    ghostCell_.apply(field, ghostWeights_, config_);
}

void IB::applySST(const Field& field,double laminarMu,std::vector<double>& k,
        std::vector<double>& omega,std::vector<double>& muT) const {
    if (!usesGhostCells() || config_.wallClosure!=FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic)
        throw std::runtime_error("SST requires the bound stationary no-slip immersed wall.");
    if (wallField_!=&field) throw std::runtime_error("SST immersed wall is bound to a different physical patch.");
    GhostIBM::applySSTWall(field,ghostWeights_.viscousPlans(),laminarMu,k,omega,muT);
}

FDM::ImmersedConstraintResult IB::applyConstraint(
        const std::vector<Field*>& fields, double targetTime, double dt) {
    if (!usesForcing()) return {};
    Common::validateStageTime(targetTime, dt);
    return forcing_.projectPredictedState(fields, targetTime, dt);
}

void IB::setExecutionRuntime(FDM::IExecutionRuntime* runtime) {
    // `IB` 是 application 与具体 forcing method 的唯一薄转发边界。
    // Runtime 只携带 freshness/reduction/COPY 能力，IBM 方法本身不获取
    // communicator 或 rank。
    forcing_.setExecutionRuntime(runtime);
}

const FDM::ImmersedSurfaceSystem& IB::prepareMonolithicSystem(
        Field& field, double targetTime, double dt) {
    if (!usesForcing()) {
        throw std::runtime_error("IBM forcing is not active.");
    }
    Common::validateStageTime(targetTime, dt);
    return forcing_.prepareMonolithicSystem(field, targetTime, dt);
}

FDM::ImmersedConstraintResult IB::acceptMonolithicSolution(
        Field& field, double targetTime, double dt,
        const FDM::ImmersedKKTState& state) {
    if (!usesForcing()) {
        throw std::runtime_error("IBM forcing is not active.");
    }
    Common::validateStageTime(targetTime, dt);
    return forcing_.acceptMonolithicSolution(field, targetTime, dt, state);
}

void IB::publishDiagnostics(Field& field) {
    if(usesForcing()) forcing_.publishDiagnostics(field);
}

double IB::constraintMask(
        const Field& field, int i, int j, int k) const {
    return usesForcing() ? forcing_.constraintMask(field, i, j, k) : 0.0;
}

Vector3 IB::multiplier(
        const Field& field, int i, int j, int k) const {
    return usesForcing()
        ? forcing_.multiplier(field, i, j, k) : Vector3();
}

FDM::ImmersedStorageView IB::storageView(FDM::ImmersedStorageKind kind) const {
    return usesForcing()?forcing_.storageView(kind):FDM::ImmersedStorageView{};
}

} // namespace SF::IBM
