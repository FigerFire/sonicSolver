/// @file SF_singleFluid.cpp
/// @brief 组装 conservative single-field state、providers 与 solver interfaces。
///
/// Data flow:
///   Field + resolved system + case models
///       -> StateBundle / services / equation coupling
///       -> compiled-plan-bound stepper and flowLoop
///
/// 本文件不定义 governing equations、RK 系数或 MPI ownership 规则。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.08.04-----------*/

#include "app/application/execution/SF_execution.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/execution/SF_flowLoop.h"

#include "SF_resultWriter.h"
#include "SF_IBM.h"
#include "app/application/output/SF_report.h"
#include "core/interfaces/SF_log.h"
#include "solver/equation/coupling/SF_equationCoupling.h"
#include "SF_equationSystem.h"
#include "SF_factory.h"
#include "SF_field.h"
#include "SF_homogeneousPhaseChange.h"
#include "solver/equation/coupling/SF_interfaceCoupling.h"
#include "SF_interfaceModel.h"
#include "models/physics/interfaceModel/levelSet/SF_state.h"
#include "SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include "SF_numericsPolicy.h"
#include "SF_multiphase.h"
#include "app/application/output/SF_fields.h"
#include "SF_pipeline.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
#include "app/application/adapters/SF_ibmAdapters.h"
#include "solver/algorithm/time/SF_time.h"
#include "SF_turbulence.h"
#include "core/state/SF_state.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

namespace SF::Application::Execution::Detail {

int executeConservativeEquations(
        Field& field,
        ResultWriter& writer,
        Parallel::ParallelContext& parallel,
        IBM::IB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig,
        bool ibmEnabled,
        bool initialOutputOnly,
        int localBlockId) {
    using Report::broadcastSolverConfig;
    using Report::formatRunControl;
    using Report::formatTimeStepStatus;
    using Report::multiPhaseSummary;
    using Equation::Coupling::CompositeTransportProvider;
    using Equation::Coupling::HomogeneousPhaseChangeProvider;
    using Equation::Coupling::InterfaceEquationProvider;
    using Equation::Coupling::MixtureEquationProvider;
    using Equation::Coupling::registerInterfaceState;
    using Equation::Coupling::registerMixtureState;
    using Output::multiPhaseVTKScalars;
    using Output::interfaceVTKScalars;
    using Adapters::IBBoundaryAdapter;
    using Adapters::IBConstraintAdapter;

    const bool ghostIBMActive =
        System::requiresProvider(system,"ibm.boundary");
    const bool forcingIBMActive =
        System::requiresProvider(system,"ibm.constraint");
    if (ghostIBMActive != (ibmEnabled && ibm.usesGhostCells())
        || forcingIBMActive != (ibmEnabled && ibm.usesForcing())) {
        throw std::runtime_error(
            "Resolved IBM provider requirements do not match initialized IBM resources.");
    }

    SF::Physics::Multiphase::MultiPhaseModel multiPhaseModel;
    std::unique_ptr<SF::Physics::InterfaceModels::Model> interfaceModel;
    const bool legacyMultiPhaseActive =
        System::requiresProvider(system,"equation.legacy-mixture");
    const bool interfaceActive =
        System::requiresProvider(system,"equation.level-set");
    const bool usesHomogeneousThermodynamics =
        System::requiresProvider(system,"thermodynamics.homogeneous");
    const bool multiPhaseActive = legacyMultiPhaseActive || interfaceActive;
    if (interfaceActive && !solverConfig.numerics.timeRecipeDeclared) {
        throw std::runtime_error(
            "OneFluidInterface requires an explicit "
            "system/fvSchemes ddtSchemes/default entry.");
    }
    std::shared_ptr<const SF::Physics::FluidStateModel::HomogeneousMultiphaseStateModel>
        homogeneousEquations;
    std::shared_ptr<const SF::Physics::FluidStateModel::SingleFluidStateModel>
        singleFluidEquations;
    if (usesHomogeneousThermodynamics) {
        homogeneousEquations = SF::Physics::FluidStateModel::makeHomogeneous(
            caseConfig.multiPhase);
        SF::Physics::FluidStateModel::initializeHomogeneousField(
            field, caseConfig.multiPhase, homogeneousEquations);
        SF::broadcast("FluidStateModel       : ",
                      "homogeneousMultiphase (generic thermodynamics, primary="
                      + std::to_string(homogeneousEquations->variableCount())
                      + ")");
    }
    if ((!legacyMultiPhaseActive && !usesHomogeneousThermodynamics)
        || interfaceActive) {
        singleFluidEquations=
            SF::Physics::FluidStateModel::makeSingleFluidPerfectGas(
                solverConfig.numerics.idealGasGamma,
                solverConfig.numerics.idealGasConstant,
                solverConfig.numerics.dynamicViscosity,
                solverConfig.numerics.prandtl);
        if (interfaceActive) {
            // OneFluidInterface still advances the conservative density
            // system.  Its level-set coupling is auxiliary state, so bind the
            // existing PerfectGas EOS without reinitializing conservative Q.
            field.setStateModel(singleFluidEquations);
            SF::broadcast("FluidStateModel       : ",
                          "singleFluid/perfectGas bound to one-fluid interface");
        } else if (ghostIBMActive || solverConfig.boundaries.ilwEnabled) {
            // IBM/ILW setup owns its existing conservative initialization.  The
            // FluidStateModel is still the density solver's one thermodynamic
            // binding and must be attached before the first timestep.
            field.setStateModel(singleFluidEquations);
            SF::broadcast("FluidStateModel       : ",
                          "singleFluid/perfectGas attached to IBM/ILW state");
        } else if (caseConfig.time.startTime > 0.0) {
            if (solverConfig.initial.energy.empty()) {
                SF::broadcast(
                    "Fatal restart: ",
                    "startTime>0 requires an explicit rhoE field in the "
                    "selected time directory.");
                return -1;
            }
            SF::Physics::FluidStateModel::initializeSingleFluidRestart(
                field,singleFluidEquations);
            SF::broadcast(
                "FluidStateModel       : ",
                "singleFluid/perfectGas restart Q closure");
        } else {
            SF::Physics::FluidStateModel::initializeSingleFluidField(
                field, singleFluidEquations,
                solverConfig.initial.pressure,
                solverConfig.initial.temperature);
            SF::broadcast(
                "FluidStateModel       : ",
                "singleFluid/perfectGas cold start");
        }
    }
    if (legacyMultiPhaseActive) {
        multiPhaseModel.configure(caseConfig.multiPhase);
        multiPhaseModel.configureBoundaryNumerics(
            solverConfig.boundaries.ilwEnabled,
            solverConfig.boundaries.ilwOrder);
        multiPhaseModel.initialize(field);
        multiPhaseModel.initializeConservedFromPrimitive(field);
    } else if (interfaceActive) {
        interfaceModel = SF::Physics::InterfaceModels::makeInterfaceModel(
            caseConfig.multiPhase);
        interfaceModel->initialize(field);
    }
    if (multiPhaseActive) {
        SF::broadcast("MultiPhase model  : ",
                      multiPhaseSummary(
                          caseConfig.multiPhase,
                          caseConfig.compatFlowLabel));
    }
    SF::State::StateBundle stateBundle;
    stateBundle.patches = {&field};
    stateBundle.stateModel = field.stateModel();
    stateBundle.time = caseConfig.time.startTime;

    std::unique_ptr<SF::Physics::FluidStateModel::HomogeneousPhaseChange>
        homogeneousPhaseChange;
    std::unique_ptr<HomogeneousPhaseChangeProvider>
        homogeneousCoupling;
    if (homogeneousEquations && caseConfig.multiPhase.phaseChange.enabled) {
        homogeneousPhaseChange = std::make_unique<
            SF::Physics::FluidStateModel::HomogeneousPhaseChange>(
                caseConfig.multiPhase, homogeneousEquations);
        homogeneousPhaseChange->initialize(field);
        homogeneousCoupling =
            std::make_unique<HomogeneousPhaseChangeProvider>(
                *homogeneousPhaseChange);
        SF::broadcast("Phase change      : ",
                      caseConfig.multiPhase.phaseChange.model
                      + " -> partialDensity, rhoE source=0");
    }
    SF::Turbulence::Manager turbulenceManager(solverConfig.turbulence);
    CompositeTransportProvider coupledTransport;
    auto& parallelCoordinator = parallel.coordinator();
    // A pressure-based serial case initializes an MPI-capable backend for
    // HYPRE, but one process still has serial ownership semantics.
    ::SF::Execution::Runtime executionRuntime(
        parallel.active() && parallel.size() > 1
            ? &parallelCoordinator : nullptr);
    IBBoundaryAdapter immersedBoundary(&ibm);
    IBConstraintAdapter immersedConstraint(&ibm);
    Boundary::Applicator physicalBoundary(solverConfig.boundaries);
    const int conservativeHaloDepth = std::max(
        system.numericalSystem.requiredHaloWidth,
        FDM::requiredGhostLayersForILW(solverConfig.numerics.ilwOrder));
    Boundary::Pipeline boundaryPipeline({
        &physicalBoundary,
        &executionRuntime,
        ghostIBMActive ? &immersedBoundary : nullptr,
        ghostIBMActive && parallel.active(),
        conservativeHaloDepth});
    SF::FDM::FunctionObserver solverObserver([](const SF::FDM::SolverMessage& message) {
        if (message.kind == SF::FDM::SolverMessageKind::TimeStep) {
            SF::broadcast("Step time: ",
                          formatTimeStepStatus(message.time, message.dt));
        } else if (message.kind == SF::FDM::SolverMessageKind::StateClosure) {
            SF::broadcast("State closure: ", message.detail);
        } else if (message.kind == SF::FDM::SolverMessageKind::Setup
                   && message.topic == "Convection contract") {
            SF::broadcast("Convection contract: ", message.detail);
        }
    });

    FDM::ImmersedCouplingPorts immersedPorts;
    if (ghostIBMActive) immersedPorts.boundary = &immersedBoundary;
    if (forcingIBMActive) {
        immersedPorts.constraint = &immersedConstraint;
        immersedPorts.system = &ibm;
    }
    std::vector<double> multiPhaseAlphaRHS;
    std::unique_ptr<MixtureEquationProvider>
        legacyMultiphaseCoupling;
    std::unique_ptr<InterfaceEquationProvider> interfaceCoupling;
    if (legacyMultiPhaseActive) {
        coupledTransport.setMultiPhase(&multiPhaseModel);
        registerMixtureState(
            multiPhaseModel, field,
            multiPhaseAlphaRHS, stateBundle.transported);
        legacyMultiphaseCoupling =
            std::make_unique<MixtureEquationProvider>(
                multiPhaseModel, multiPhaseAlphaRHS,
                stateBundle.transported, executionRuntime);
    } else if (interfaceActive) {
        registerInterfaceState(
            *interfaceModel, field, stateBundle.transported);
        interfaceCoupling = std::make_unique<InterfaceEquationProvider>(
            *interfaceModel, stateBundle.transported, executionRuntime);
        coupledTransport.setInterfaceModel(interfaceModel.get());
    }
    const bool turbulenceActive = turbulenceManager.initialize(field);
    const bool transportedTurbulence = turbulenceActive
        && solverConfig.turbulence.family == FDM::TurbulenceFamily::RAS
        && (solverConfig.turbulence.model
                == FDM::TurbulenceModelKind::kEpsilon
            || solverConfig.turbulence.model
                == FDM::TurbulenceModelKind::kOmegaSST);
    if (transportedTurbulence
        != System::requiresProvider(system,"equation.turbulence-transport")) {
        throw std::runtime_error(
            "Resolved turbulence-equation provider does not match the active density "
            "turbulence equations.");
    }
    if (turbulenceActive) {
        coupledTransport.setTurbulence(&turbulenceManager);
        SF::broadcast("Turbulence active : ", turbulenceManager.description());
    } else if (solverConfig.turbulence.enabled) {
        SF::broadcast("Turbulence warning: ", "enabled but no valid turbulence model was created.");
    }
    broadcastSolverConfig(solverConfig);
    SF::broadcast("Run control: ", formatRunControl(caseConfig));

    auto turbulenceVTKScalars = [&]() {
        std::vector<SF::ResultWriter::ScalarField> scalars;
        if (!turbulenceManager.active()) return scalars;
        if (!solverConfig.turbulence.enabled
            || solverConfig.turbulence.family == SF::FDM::TurbulenceFamily::None
            || solverConfig.turbulence.family == SF::FDM::TurbulenceFamily::DNS) {
            return scalars;
        }

        const auto& state = turbulenceManager.scalarFields();
        if (solverConfig.turbulence.family == SF::FDM::TurbulenceFamily::RAS) {
            scalars.push_back({"TurbulenceK", [&](int i, int j, int k) {
                return state.K(i, j, k);
            }});
            if (solverConfig.turbulence.model
                == SF::FDM::TurbulenceModelKind::kOmegaSST) {
                scalars.push_back({"TurbulenceOmega", [&](int i, int j, int k) {
                    return state.Omega(i, j, k);
                }});
            } else if (solverConfig.turbulence.model
                       == SF::FDM::TurbulenceModelKind::kEpsilon) {
                scalars.push_back({"TurbulenceEpsilon", [&](int i, int j, int k) {
                    return state.Epsilon(i, j, k);
                }});
            }
        }
        scalars.push_back({"TurbulentViscosity", [&](int i, int j, int k) {
            return state.EddyMu(i, j, k);
        }});
        return scalars;
    };

    auto solverVTKScalars = [&]() {
        std::vector<SF::ResultWriter::ScalarField> scalars = turbulenceVTKScalars();
        if (forcingIBMActive) {
            scalars.push_back({"IBMConstraintMask", [&](int i, int j, int k) {
                return ibm.constraintMask(field, i, j, k);
            }});
            scalars.push_back({"IBMMultiplierX", [&](int i, int j, int k) {
                return ibm.multiplier(field, i, j, k).x;
            }});
            scalars.push_back({"IBMMultiplierY", [&](int i, int j, int k) {
                return ibm.multiplier(field, i, j, k).y;
            }});
            scalars.push_back({"IBMMultiplierZ", [&](int i, int j, int k) {
                return ibm.multiplier(field, i, j, k).z;
            }});
        }
        if (singleFluidEquations) {
            scalars.push_back({"MomentumX", [&](int i, int j, int k) {
                return field(i,j,k,RU);
            }});
            scalars.push_back({"MomentumY", [&](int i, int j, int k) {
                return field(i,j,k,RV);
            }});
            scalars.push_back({"MomentumZ", [&](int i, int j, int k) {
                return field(i,j,k,RW);
            }});
            scalars.push_back({"TotalEnergyDensity", [&](int i, int j, int k) {
                return field(i,j,k,E);
            }});
            scalars.push_back({"Temperature", [&](int i, int j, int k) {
                return field.thermodynamicState(i,j,k).temperature;
            }});
            scalars.push_back({"SoundSpeed", [&](int i, int j, int k) {
                return field.thermodynamicState(i,j,k).soundSpeed;
            }});
        }
        std::vector<SF::ResultWriter::ScalarField> phaseScalars;
        if (legacyMultiPhaseActive) {
            phaseScalars = multiPhaseVTKScalars(
                multiPhaseModel, &stateBundle.transported);
        } else if (interfaceActive) {
            phaseScalars = interfaceVTKScalars(
                *interfaceModel, &stateBundle.transported);
        }
        scalars.insert(scalars.end(), phaseScalars.begin(), phaseScalars.end());
        if (homogeneousEquations) {
            const auto& variables = homogeneousEquations->primaryVariables();
            for (int v = 0; v < homogeneousEquations->densityVariableCount(); ++v) {
                scalars.push_back({variables[(size_t)v].name,
                                   [&, v](int i, int j, int k) { return field(i,j,k,v); }});
            }
            scalars.push_back({"Temperature", [&](int i, int j, int k) {
                return field.thermodynamicState(i,j,k).temperature;
            }});
            scalars.push_back({"SoundSpeed", [&](int i, int j, int k) {
                return field.thermodynamicState(i,j,k).soundSpeed;
            }});
            for (int phase = 0; phase < homogeneousEquations->phaseCount(); ++phase) {
                const std::string phaseName = homogeneousEquations->phaseName(phase);
                scalars.push_back({"phaseMass." + phaseName,
                                   [&, phase](int i, int j, int k) {
                    return field.thermodynamicState(i,j,k).phaseMass[(size_t)phase];
                }});
                scalars.push_back({"alpha." + phaseName,
                                   [&, phase](int i, int j, int k) {
                    return field.thermodynamicState(i,j,k).volumeFraction[(size_t)phase];
                }});
            }
        }
        return scalars;
    };

    auto saveStepVTK = [&](int outStep) {
        auto scalars = solverVTKScalars();
        if (parallel.active()) writer.save(field, outStep, localBlockId, scalars);
        else               writer.save(field, outStep, scalars);
    };
    auto saveTimeVTK = [&](double outTime) {
        auto scalars = solverVTKScalars();
        if (parallel.active()) writer.save(field, outTime, localBlockId, scalars);
        else               writer.save(field, outTime, scalars);
    };

    stateBundle.distributed.add(State::conservativeView(
        "conservative", localBlockId, field));
    if (!stateBundle.transported.empty()) {
        stateBundle.transported.registerDistributed(
            stateBundle.distributed, localBlockId, field);
    }
    if (interfaceActive) {
        auto* levelSet = interfaceModel->levelSetState();
        if (levelSet) {
            stateBundle.distributed.add(State::workspaceView(
                "levelSetCurvature", localBlockId, field,
                levelSet->curvatures(), 1));
        }
    }
    if (turbulenceManager.active()) {
        auto& turbulence = turbulenceManager.scalarFields();
        const auto names = turbulenceManager.distributedWriteFields();
        for (const auto& name : names) {
            const Turbulence::ScalarSlot slot = name == "k"
                ? Turbulence::ScalarSlot::K
                : name == "epsilon" ? Turbulence::ScalarSlot::Epsilon
                : name == "omega" ? Turbulence::ScalarSlot::Omega
                : Turbulence::ScalarSlot::EddyMu;
            stateBundle.distributed.add(State::workspaceView(
                name, localBlockId, field,
                turbulence.values(slot), 1));
        }
    }

    // 初值先由各 rank 独立闭合为物理保守量，再通过 Runtime 声明其 stencil
    // 需求。不能在 application 装配网格后立刻交换，否则 canonical owner
    // 仍持有 Field 的零初始化保守量，MPI-4 等更细分区会显式暴露该顺序错误。
    executionRuntime.attachState(stateBundle);
    executionRuntime.prepare({
        "initialized conservative-state halo",
        {::SF::Execution::readHalo("conservative", conservativeHaloDepth)}});

    if (initialOutputOnly) {
        if (caseConfig.time.writeByStep) saveStepVTK(0);
        else            saveTimeVTK(caseConfig.time.startTime);
        return 0;
    }

    if (caseConfig.writeInitial) {
        if (caseConfig.time.writeByStep) saveStepVTK(0);
        else            saveTimeVTK(caseConfig.time.startTime);
    }

    FDM::IEquationSystemCoupling* equationProvider = nullptr;
    if (homogeneousCoupling) equationProvider = homogeneousCoupling.get();
    if (legacyMultiphaseCoupling) {
        equationProvider = legacyMultiphaseCoupling.get();
    }
    if (interfaceCoupling) equationProvider = interfaceCoupling.get();
    return runConservative(
        solverConfig,system,plan,stateBundle,boundaryPipeline,
        executionRuntime,solverObserver,immersedPorts,
        coupledTransport.active() ? &coupledTransport : nullptr,
        equationProvider,caseConfig.time,
        "Unified single-field",
        saveStepVTK,
        saveTimeVTK,
        [&](bool finished) {
            return parallelCoordinator.allRanksAgree(finished);
        });
}

} // namespace SF::Application::Execution::Detail
