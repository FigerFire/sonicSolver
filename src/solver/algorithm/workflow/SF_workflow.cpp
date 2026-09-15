/// @file SF_workflow.cpp
/// @brief 从求解器、相系统和界面选择构造顶层工作流。

#include "SF_workflow.h"

#include <stdexcept>

namespace SF::Workflow {
namespace {

/// @brief 由 SolveBlock 的 strategy 推断粗粒度阶段类别。
/// @details 该映射是临时迁移层：runner 真正消费 Stage 后，应由
///          WorkflowBuilder 直接产出精确阶段，而不是从字符串 strategy 反推。
StageKind classifyStage(
        const System::SolveBlock& block,
        const std::string& timeIntegrator) {
    const std::string& strategy = block.strategy;
    if (strategy == timeIntegrator) {
        return StageKind::Predictor;
    }
    if (strategy.find("correction") != std::string::npos
        || strategy.find("PIMPLE") != std::string::npos
        || strategy.find("pimple") != std::string::npos) {
        return StageKind::Correction;
    }
    if (strategy.find("predictor") != std::string::npos) {
        return StageKind::Predictor;
    }
    if (strategy.find("constraint") != std::string::npos
        || strategy.find("KKT") != std::string::npos
        || strategy.find("monolithic") != std::string::npos) {
        return StageKind::Constraint;
    }
    return StageKind::Commit;
}

std::vector<Stage> buildStages(const System::ResolvedSimulationSystem& system) {
    std::vector<Stage> stages;
    stages.reserve(system.solveBlocks.size());
    for (const auto& block : system.solveBlocks) {
        Stage stage;
        stage.id = block.id;
        stage.kind = classifyStage(block, system.timeIntegrator);
        stage.equations = block.equations;
        stage.constraints = block.constraints;
        stages.push_back(std::move(stage));
    }
    return stages;
}

} // namespace

std::string Plan::name() const {
    if (stages.empty()) return "empty";
    std::string result;
    for (const auto& stage : stages) {
        if (!result.empty()) result += " -> ";
        result += stage.id;
    }
    return result;
}

Plan makePlan(
        const System::ResolvedSimulationSystem& system,
        const FDM::SolverConfig& config) {
    (void)config;
    Plan plan;
    plan.stages = buildStages(system);
    plan.timeIntegrator = system.timeIntegrator;
    return plan;
}

void validateMeshExecution(
        const System::ResolvedSimulationSystem& system,
        const FDM::SolverConfig& config,
        bool multiField,
        bool hasCoupledInterfaces) {
    if (!multiField) return;
    if (System::hasEquation(system, "E_PHASE_MASS")) {
        throw std::runtime_error(
            "homogeneous EquationSet has not migrated to the multi-field "
            "MPI stepper.");
    }
    if (System::hasEquationPrefix(system, "E_CONTINUITY.")) {
        throw std::runtime_error(
            "eulerianEulerian requires one structured Field per MPI rank; "
            "multi-patch-per-rank execution is not implemented.");
    }
    if (config.turbulence.enabled
        && config.turbulence.family != FDM::TurbulenceFamily::None
        && config.turbulence.family != FDM::TurbulenceFamily::DNS) {
        throw std::runtime_error(
            "Multi-field MPI has not migrated transported turbulence scalar "
            "state and halo plans.");
    }
    if (System::hasEquation(system, "E_LEGACY_ALPHA")
        && System::requiresCapability(system, "CanonicalScalarInterfaceFlux")
        && hasCoupledInterfaces) {
        throw std::runtime_error(
            "Legacy transported alpha on coupled patches has no canonical "
            "scalar interface flux.");
    }
}

} // namespace SF::Workflow
