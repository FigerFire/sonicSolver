/// @file SF_densityBasedRHS.cpp
/// @brief Density-based stage RHS 的共享编排；不拥有 RK 或 timestep commit。

#include "solver/algorithm/SF_conservativeRHS.h"
#include "solver/algorithm/SF_highOrderTrace.h"

#include "SF_equationState.h"
#include "SF_numericsPolicy.h"
#include "methods/numerics/structured/SF_structured.h"
#include "solver/system/SF_numericalCompiler.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SF::SolverAlgorithm::ConservativeRHS {
namespace {

const FDM::ITransportModel* transportFor(
        const Field& field, const FDM::SolverServices& services) {
    if (services.transportModel) return services.transportModel;
    return services.equationSystem
        ? services.equationSystem->transportModel(field) : nullptr;
}

void addWriteOnce(std::vector<Execution::FieldAccess>& writes,
                  const std::string& field) {
    for (const auto& access : writes) {
        if (access.field == field) return;
    }
    writes.push_back(Execution::writeOwned(field));
}

} // namespace

void validateTimestepState(
        const State::StateBundle& state,
        const FDM::SolverConfig& config,
        const System::CompiledNumericalSystem& numericalSystem) {
    state.validateDensityEquationBinding();
    const auto& convection = System::NumericalCompiler::requireUniqueRecipe(
        numericalSystem,FDM::TermRole::Convection);
    const auto gamma = state.stateModel->perfectGasGamma();
    if (!gamma || state.stateModel->variableCount() != 5) {
        std::ostringstream message;
        message
            << "Density-based reconstructed convection '"
            << FDM::toString(convection.id())
            << "' requires a five-variable PerfectGas FluidStateModel; selected "
            << "FluidStateModel family="
            << static_cast<int>(state.stateModel->family())
            << ", variables=" << state.stateModel->variableCount()
            << ". No whole-field Rusanov bypass is permitted.";
        throw std::runtime_error(
            message.str());
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(*gamma),
                    std::abs(config.numerics.idealGasGamma)});
    if (std::abs(*gamma - config.numerics.idealGasGamma) > tolerance) {
        std::ostringstream message;
        message << "Density-based reconstructed convection '"
                << FDM::toString(convection.id())
                << "' received PerfectGas FluidStateModel gamma=" << *gamma
                << " but configured characteristic flux gamma="
                << config.numerics.idealGasGamma
                << "; no fixed-gamma fallback is permitted.";
        throw std::runtime_error(message.str());
    }
}

void prepareBoundaryState(
        const std::vector<Field*>& fields,
        int conservativeHaloDepth,
        Boundary::Applicator& boundaryApplicator,
        FDM::SolverServices& services,
        double time,
        double dt) {
    if (services.boundaryPipeline) {
        services.boundaryPipeline->prepare(fields, time, dt);
        return;
    }
    for (Field* field : fields) {
        if (!field) continue;
        field->invalidateThermodynamicCache();
        boundaryApplicator.apply(*field);
        field->invalidateThermodynamicCache();
    }
    const Execution::OperatorContract halo{
        "boundary-state halo",
        {Execution::readHalo("conservative", conservativeHaloDepth)}};
    if (services.executionRuntime) services.executionRuntime->prepare(halo);
    if (services.immersed.boundary) services.immersed.boundary->apply(fields, time, dt);
    if (services.immersed.boundary && services.executionRuntime) {
        services.executionRuntime->finalize({
            "IBM ghost-state update", {Execution::writeOwned("conservative")}});
        services.executionRuntime->prepare(halo);
    }
}

void assembleAllPatches(
        const std::vector<Field*>& fields,
        std::vector<PatchWorkspace>& workspaces,
        const FDM::SolverConfig& config,
        const System::CompiledNumericalSystem& numericalSystem,
        Equation::Compressible::System& equations,
        State::StateBundle& state,
        FDM::SolverServices& services,
        Boundary::Applicator& boundaryApplicator,
        double stageTime) {
    if (fields.empty() || fields.size() != workspaces.size()) {
        throw std::runtime_error(
            "Density-based RHS requires one workspace for every patch.");
    }
    const bool trace = HighOrderTrace::activeFor(state.step);
    std::vector<PatchWorkspace*> traceWorkspaces;
    if (trace) {
        HighOrderTrace::beginStage(state.step, stageTime);
        traceWorkspaces.reserve(workspaces.size());
        for (std::size_t index = 0; index < workspaces.size(); ++index) {
            if (!fields[index]) {
                throw std::runtime_error(
                    "high-order trace found a null patch with a workspace slot.");
            }
            traceWorkspaces.push_back(&workspaces[index]);
            std::cout << "[SF TRACE] workspace patch=" << index
                      << " field=" << static_cast<const void*>(fields[index])
                      << " workspace=" << static_cast<const void*>(&workspaces[index])
                      << '\n';
        }
    }
    prepareBoundaryState(
        fields,numericalSystem.requiredHaloWidth,boundaryApplicator,services,
                         stageTime, state.dt);
    if (services.equationSystem) services.equationSystem->prepareRHS(fields, state.dt);

    const Execution::OperatorContract convectionContract{
        "compressible convection",
        {Execution::readHalo(
            "conservative",numericalSystem.requiredHaloWidth)}};
    if (services.executionRuntime) services.executionRuntime->prepare(convectionContract);

    const auto& convectionRecipe =
        System::NumericalCompiler::requireUniqueRecipe(
            numericalSystem,FDM::TermRole::Convection);
    const auto* diffusionRecipe =
        System::NumericalCompiler::findUniqueRecipe(
            numericalSystem,FDM::TermRole::Diffusion);
    const auto boundSources =
        System::NumericalCompiler::sourceKinds(numericalSystem);
    const Equation::Compressible::AssemblyContext assembly{
        &convectionRecipe,diffusionRecipe,&boundSources,&config.sources,
        state.dt,config.numerics.ibmBoundary,config.numerics.ilwOrder,
        config.numerics.dynamicViscosity,config.numerics.prandtl,
        config.numerics.idealGasGamma,config.numerics.idealGasConstant,
        nullptr,state.stateModel.get()};

    for (size_t index = 0; index < fields.size(); ++index) {
        Field* field = fields[index];
        if (!field) continue;
        auto& workspace = workspaces[index];
        workspace.ensureFor(*field);
        equations.begin(workspace.convectiveFlux, workspace.residual);
        if (trace) HighOrderTrace::workspaceAfterClear(*field, workspace, index);
        auto patchAssembly = assembly;
        patchAssembly.transport = transportFor(*field,services);
        equations.convection(
            *field,workspace.convectiveFlux,workspace.residual,patchAssembly);
    }
    if (trace) {
        HighOrderTrace::flux("candidate face flux before canonical COPY",
                             fields, traceWorkspaces);
        HighOrderTrace::residual("directional residual before canonical COPY",
                                 fields, traceWorkspaces, false, false);
    }

    // 每个 patch 必须先产生候选 face flux，才能进入唯一的 canonical-face
    // barrier。Runtime 在此处 COPY canonical F* 并做 +/- 残差装配；不能把
    // finalize 放进 patch 循环，否则跨 patch face 会失去同一个物理通量。
    if (services.executionRuntime) {
        std::vector<FluxField*> fluxes;
        fluxes.reserve(workspaces.size());
        for (auto& workspace : workspaces) fluxes.push_back(&workspace.convectiveFlux);
        std::vector<Residual*> residuals;
        residuals.reserve(workspaces.size());
        for (auto& workspace : workspaces) residuals.push_back(&workspace.residual);
        services.executionRuntime->synchronizeCanonicalFaceFluxes(fields, fluxes, residuals);
    }
    if (trace) {
        HighOrderTrace::flux("face flux after canonical COPY", fields, traceWorkspaces);
        HighOrderTrace::residual("directional residual after canonical COPY",
                                 fields, traceWorkspaces, false, false);
    }

    for (size_t index = 0; index < fields.size(); ++index) {
        Field* field = fields[index];
        if (!field) continue;
        auto patchAssembly = assembly;
        patchAssembly.transport = transportFor(*field,services);
        equations.diffusionAndSources(
            *field,workspaces[index].residual,patchAssembly);
    }
    if (services.equationSystem) {
        std::vector<Residual*> residuals;
        residuals.reserve(workspaces.size());
        for (auto& workspace : workspaces) residuals.push_back(&workspace.residual);
        services.equationSystem->assembleRHS(fields, residuals, state.dt);
    }
    if (services.executionRuntime) {
        for (size_t index = 0; index < fields.size(); ++index) {
            if (fields[index]) Math::stageLocalSpatialResidual(
                *fields[index], workspaces[index].residual);
        }
        if (trace) {
            HighOrderTrace::residual("local residual before GlobalDof SUM",
                                     fields, traceWorkspaces, true, false);
        }
        std::vector<Residual*> residuals;
        residuals.reserve(workspaces.size());
        for (auto& workspace : workspaces) residuals.push_back(&workspace.residual);
        services.executionRuntime->accumulateGlobalDofResiduals(fields, residuals);
        if (trace) {
            HighOrderTrace::residual("global residual after GlobalDof SUM",
                                     fields, traceWorkspaces, true, true);
        }
    }
}

void publishIntegratedState(
        const std::vector<Field*>& fields,
        const State::StateBundle& state,
        FDM::SolverServices& services) {
    if (!services.executionRuntime) return;
    std::vector<Execution::FieldAccess> writes{Execution::writeOwned("conservative")};
    for (const auto& variable : state.transported.variables()) {
        if (variable.rhs) addWriteOnce(writes, variable.descriptor.name);
    }
    if (services.equationSystem) {
        for (Field* field : fields) {
            if (!field) continue;
            const auto* variables = services.equationSystem->variables(*field);
            if (!variables) continue;
            for (const auto& variable : variables->variables()) {
                if (variable.rhs) addWriteOnce(writes, variable.descriptor.name);
            }
        }
    }
    services.executionRuntime->finalize({"explicit integrated-state update", writes});
}

void publishCorrectedConservativeState(
        bool wroteConservativeState, FDM::SolverServices& services) {
    if (!wroteConservativeState || !services.executionRuntime) return;
    services.executionRuntime->finalize({
        "immersed forcing corrected state", {Execution::writeOwned("conservative")}});
}

void validateStateClosure(
        const std::vector<Field*>& fields,
        const State::StateBundle& state,
        FDM::SolverServices& services,
        const char* stage,
        bool labelPatches) {
    for (size_t index = 0; index < fields.size(); ++index) {
        Field* field = fields[index];
        if (!field) continue;
        std::ostringstream label;
        label << stage;
        if (labelPatches) label << ", block=" << index;
        const std::string context = label.str();
        const auto summary = Numerics::validateEquationState(
            *field, *state.stateModel, context.c_str());
        if (!services.observer) continue;
        std::ostringstream detail;
        detail << context << ": checked=" << summary.checkedCells
               << ", min(rho)=" << summary.minDensity
               << ", min(p)=" << summary.minPressure
               << ", min(T)=" << summary.minTemperature
               << ", min(c)=" << summary.minSoundSpeed;
        services.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
            "EOS state closure", detail.str(), state.step, state.time, state.dt});
    }
}

} // namespace SF::SolverAlgorithm::ConservativeRHS
