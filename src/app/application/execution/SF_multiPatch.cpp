/// @file SF_multiPatch.cpp
/// @brief 组装 multi-patch state、equation coupling 与 distributed providers。
///
/// Data flow:
///   decomposed mesh + local patch IDs + resolved system
///       -> patch StateBundle / adapters / output views
///       -> existing stepper and flowLoop
///
/// 本文件不实现 halo/COPY/SUM backend，也不改写数值 lifecycle。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.08.04-----------*/

#include "app/application/execution/SF_execution.h"
#include "app/application/execution/SF_flowLoop.h"

#include "SF_resultWriter.h"
#include "SF_MultiBlockMesh.h"
#include "SF_boundaryGeometry.h"
#include "app/application/output/SF_report.h"
#include "core/interfaces/SF_log.h"
#include "SF_compositeIBM.h"
#include "solver/equation/coupling/SF_equationCoupling.h"
#include "solver/equation/coupling/SF_interfaceCoupling.h"
#include "SF_interfaceModel.h"
#include "models/physics/interfaceModel/levelSet/SF_state.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include "SF_numericsPolicy.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
#include "SF_multiphase.h"
#include "app/application/output/SF_fields.h"
#include "SF_pipeline.h"
#include "app/application/adapters/SF_ibmAdapters.h"
#include "solver/algorithm/time/SF_time.h"
#include "core/state/SF_state.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

namespace SF::Application::Execution {

int executeMulti(
        MultiBlockMesh& mesh,
        ResultWriter& writer,
        Parallel::ParallelContext& parallel,
        IBM::CompositeIB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig,
        bool ibmEnabled,
        bool initialOutputOnly) {
    using Equation::Coupling::CompositeTransportProvider;
    using Equation::Coupling::MultiPatchMixtureEquationProvider;
    using Equation::Coupling::MultiPatchInterfaceEquationProvider;
    using Equation::Coupling::registerInterfaceState;
    using Equation::Coupling::registerMixtureState;
    using Output::multiPhaseVTKScalars;
    using Output::interfaceVTKScalars;
    using Report::broadcastSolverConfig;
    using Report::formatRunControl;
    using Report::formatTimeStepStatus;
    using Report::multiPhaseSummary;
    using Adapters::MultiPatchIBAdapter;

    // This is capability validation against the resolved system, not a
    // second workflow plan.  It stays next to the multi-patch adapter until
    // those numerical providers are operation-level implementations.
    if (System::requiresProvider(system,"thermodynamics.homogeneous")
        || System::requiresProvider(system,"flow.eulerian-pressure")) {
        broadcast("Fatal execution capability: ",
            "the selected executable system has no multi-patch provider.");
        return -1;
    }
    if (System::requiresProvider(system,"equation.turbulence-transport")) {
        broadcast("Fatal execution capability: ",
            "multi-patch transported turbulence state is not implemented.");
        return -1;
    }
    if (System::requiresProvider(system,"equation.legacy-mixture")
        && System::requiresCapability(system, "CanonicalScalarInterfaceFlux")
        && !mesh.haloExchangePlan().interfaces.empty()) {
        broadcast("Fatal execution capability: ",
            "legacy transported alpha has no canonical scalar interface flux.");
        return -1;
    }
    for (auto& block : mesh.blocks()) {
        Boundary::Geometry::configureEmptyDimensions(
            block.field,
            solverConfig.boundaries.density,
            solverConfig.boundaries.velocity,
            solverConfig.boundaries.energyFromPressure,
            solverConfig.boundaries.thermal,
            solverConfig.turbulence.scalars.kBoundary,
            solverConfig.turbulence.scalars.epsilonBoundary,
            solverConfig.turbulence.scalars.omegaBoundary);
    }

    const bool canonicalFlux = solverConfig.numerics.interfaceFlux
        == FDM::InterfaceFluxPolicy::SharedInterfaceFlux;
    parallel.configureBlocks(
        &mesh.haloExchangePlan(), &mesh.blocks(), canonicalFlux);
    auto& coordinator = parallel.coordinator();
    std::vector<Field*> localFields;
    std::vector<int> localPatchIds;
    const MeshPartition& partition =
        mesh.partition((size_t)parallel.rank());
    for (int patchId : partition.patchIds) {
        if (patchId < 0 || patchId >= (int)mesh.size()) {
            broadcast("Fatal: ",
                      "partition contains invalid structured patch id.");
            return -1;
        }
        localPatchIds.push_back(patchId);
        localFields.push_back(&mesh.block((size_t)patchId).field);
    }

    const bool boundaryIBMRequired =
        System::requiresProvider(system,"ibm.boundary");
    if (boundaryIBMRequired != ibmEnabled) {
        throw std::runtime_error(
            "Multi-patch IBM resource does not match the compiled boundary provider requirement.");
    }
    if (boundaryIBMRequired
        && !ibm.setupLocalPatches(mesh, localPatchIds, parallel.rank())) {
        return -1;
    }
    parallel.exchangeInitial(mesh.blocks());
    if (boundaryIBMRequired) {
        ibm.applyLocalPatches(
            mesh, localPatchIds, caseConfig.time.startTime, 0.0);
        parallel.exchangeInitial(mesh.blocks());
    }

    std::vector<Physics::Multiphase::MultiPhaseModel> phaseModels;
    std::vector<std::unique_ptr<Physics::InterfaceModels::Model>>
        interfaceModels;
    std::vector<CompositeTransportProvider> transportModels;
    std::vector<std::vector<double>> phaseRHS;
    std::vector<State::VariableRegistry> variables;
    const bool legacyMultiPhaseActive =
        System::requiresProvider(system,"equation.legacy-mixture");
    const bool interfaceActive =
        System::requiresProvider(system,"equation.level-set");
    const bool multiPhaseActive =
        legacyMultiPhaseActive || interfaceActive;
    if (interfaceActive && !solverConfig.numerics.timeRecipeDeclared) {
        throw std::runtime_error(
            "OneFluidInterface requires an explicit "
            "system/fvSchemes ddtSchemes/default entry.");
    }
    if (legacyMultiPhaseActive) {
        phaseModels.resize(mesh.size());
        transportModels.resize(mesh.size());
        phaseRHS.resize(mesh.size());
        variables.resize(mesh.size());
        for (size_t blockId = 0; blockId < mesh.size(); ++blockId) {
            auto& model = phaseModels[blockId];
            model.configure(caseConfig.multiPhase);
            model.configureBoundaryNumerics(
                solverConfig.boundaries.ilwEnabled,
                solverConfig.boundaries.ilwOrder);
            model.initialize(mesh.block(blockId).field);
            model.initializeConservedFromPrimitive(
                mesh.block(blockId).field);
            transportModels[blockId].setMultiPhase(&model);
            registerMixtureState(
                model, mesh.block(blockId).field,
                phaseRHS[blockId], variables[blockId]);
        }
    } else if (interfaceActive) {
        interfaceModels.resize(mesh.size());
        variables.resize(mesh.size());
        for (size_t blockId = 0; blockId < mesh.size(); ++blockId) {
            auto& field = mesh.block(blockId).field;
            interfaceModels[blockId] =
                Physics::InterfaceModels::makeInterfaceModel(
                    caseConfig.multiPhase);
            interfaceModels[blockId]->initialize(field);
            registerInterfaceState(
                *interfaceModels[blockId], field, variables[blockId]);
        }
    }
    if (multiPhaseActive) {
        broadcast("MultiPhase model  : ",
                  multiPhaseSummary(
                      caseConfig.multiPhase,
                      caseConfig.compatFlowLabel)
                  + " (blocks=" + std::to_string(mesh.size())
                  + ", local patches="
                  + std::to_string(localFields.size()) + ")");
    }

    ::SF::Execution::Runtime executionRuntime(&coordinator);
    FDM::FunctionObserver observer([](const FDM::SolverMessage& message) {
        if (message.kind == FDM::SolverMessageKind::TimeStep) {
            broadcast("Step time: ",
                      formatTimeStepStatus(message.time, message.dt));
        } else if (message.kind
                   == FDM::SolverMessageKind::StateClosure) {
            broadcast("State closure: ", message.detail);
        } else if (message.kind == FDM::SolverMessageKind::Setup
                   && message.topic == "Convection contract") {
            broadcast("Convection contract: ", message.detail);
        }
    });

    std::unique_ptr<MultiPatchIBAdapter> ibmAdapter;
    if (boundaryIBMRequired) {
        ibmAdapter = std::make_unique<MultiPatchIBAdapter>(
                &ibm, &mesh, &localPatchIds);
    }
    Boundary::Applicator physicalBoundary(solverConfig.boundaries);
    const int conservativeHaloDepth = std::max(
        system.numericalSystem.requiredHaloWidth,
        FDM::requiredGhostLayersForILW(solverConfig.numerics.ilwOrder));
    Boundary::Pipeline boundaryPipeline({
        &physicalBoundary,
        &executionRuntime,
        ibmAdapter.get(),
        ibmAdapter != nullptr,
        conservativeHaloDepth});
    std::unique_ptr<MultiPatchMixtureEquationProvider> equations;
    std::unique_ptr<MultiPatchInterfaceEquationProvider>
        interfaceEquations;
    if (legacyMultiPhaseActive) {
        equations = std::make_unique<MultiPatchMixtureEquationProvider>(
            mesh, localPatchIds, executionRuntime, phaseModels, phaseRHS,
            variables, transportModels);
    } else if (interfaceActive) {
        interfaceEquations =
            std::make_unique<MultiPatchInterfaceEquationProvider>(
                mesh, localPatchIds, executionRuntime,
                interfaceModels, variables);
    }

    broadcastSolverConfig(solverConfig);
    broadcast("Run control: ", formatRunControl(caseConfig));
    broadcast("MPI partition ownership: ",
              "rank " + std::to_string(parallel.rank()) + " owns "
              + std::to_string(localFields.size())
              + " internal structured patch(es)");

    std::vector<ResultWriter::FieldPiece> pieces;
    pieces.reserve(mesh.size());
    for (size_t patchId = 0; patchId < mesh.size(); ++patchId) {
        pieces.push_back({
            &mesh.block(patchId).field,
            &mesh.block(patchId).globalPointIds,
            &mesh.block(patchId).topology,
            (int)patchId,
            legacyMultiPhaseActive
                ? multiPhaseVTKScalars(
                    phaseModels[patchId], &variables[patchId])
                : interfaceActive
                    ? interfaceVTKScalars(
                        *interfaceModels[patchId], &variables[patchId])
                    : std::vector<ResultWriter::ScalarField>{}
        });
    }
    auto saveStep = [&](int value) { writer.savePieces(pieces, value); };
    auto saveTime = [&](double value) { writer.savePieces(pieces, value); };
    if (initialOutputOnly) {
        if (caseConfig.time.writeByStep) saveStep(0);
        else saveTime(caseConfig.time.startTime);
        return 0;
    }
    if (caseConfig.writeInitial) {
        if (caseConfig.time.writeByStep) saveStep(0);
        else saveTime(caseConfig.time.startTime);
    }

    State::StateBundle bundle;
    bundle.patches = localFields;
    if (System::requiresProvider(system,"thermodynamics.single-fluid")) {
        bundle.stateModel = Physics::FluidStateModel::makeSingleFluidPerfectGas(
            solverConfig.numerics.idealGasGamma,
            solverConfig.numerics.idealGasConstant,
            solverConfig.numerics.dynamicViscosity,
            solverConfig.numerics.prandtl);
        for (Field* field : bundle.patches) {
            if (field) field->setStateModel(bundle.stateModel);
        }
    }
    bundle.time = caseConfig.time.startTime;
    for (size_t blockId = 0; blockId < mesh.size(); ++blockId) {
        auto& geometry = mesh.block(blockId).field;
        bundle.distributed.add(State::conservativeView(
            "conservative", (int)blockId, geometry));
        if (blockId < variables.size()) {
            variables[blockId].registerDistributed(
                bundle.distributed, (int)blockId, geometry);
        }
        if (interfaceActive) {
            auto* levelSet = interfaceModels[blockId]->levelSetState();
            if (levelSet) {
                bundle.distributed.add(State::workspaceView(
                    "levelSetCurvature", (int)blockId, geometry,
                    levelSet->curvatures(), 1));
            }
        }
    }
    FDM::IEquationSystemCoupling* equationProvider = nullptr;
    if (equations) equationProvider = equations.get();
    if (interfaceEquations) equationProvider = interfaceEquations.get();
    return Detail::runConservative(
        solverConfig,system,plan,bundle,boundaryPipeline,executionRuntime,
        observer,
        ibmAdapter
            ? FDM::ImmersedCouplingPorts{ibmAdapter.get(),nullptr,nullptr}
            : FDM::ImmersedCouplingPorts{},
        nullptr,equationProvider,caseConfig.time,"Multi-patch density",
        saveStep,
        saveTime,
        [&](bool finished) { return coordinator.allRanksAgree(finished); });
}

} // namespace SF::Application::Execution
