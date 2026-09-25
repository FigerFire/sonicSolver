#!/usr/bin/env python3
"""Check quoted internal includes against SonicSolver's dependency baseline."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"]+)"')
SOURCE_SUFFIXES = {".h", ".hpp", ".cpp", ".cc", ".cxx", ".inl"}
HEADER_SUFFIXES = {".h", ".hpp"}


def relative(root: Path, path: Path) -> str:
    return path.relative_to(root).as_posix()


def module(path: str) -> str:
    parts = Path(path).parts
    if len(parts) < 3 or parts[0] != "src":
        return ""
    if parts[1] == "solver" and len(parts) >= 3:
        return "/".join(parts[1:3])
    return parts[1]


def dependency_rule(source: str, target: str) -> str | None:
    source_module, target_module = module(source), module(target)
    if source_module == "core" and target_module.startswith("solver/"):
        return "core -> solver"
    if source_module == "core" and target_module == "models":
        return "core -> models"
    if source_module == "methods" and target_module.startswith("solver/"):
        return "methods -> solver"
    if source_module == "solver/linearAlgebra" and target_module in {"solver/equation", "solver/algorithm"}:
        return "solver/linearAlgebra -> upper solver layer"
    if source_module in {"solver/discretization", "solver/boundary"} and target_module in {"solver/equation", "solver/algorithm"}:
        return f"{source_module} -> upper solver layer"
    if source_module == "models" and target_module in {"solver/discretization", "solver/boundary"}:
        return f"models -> {target_module}"
    if source_module == "infrastructure" and target_module.startswith("solver/"):
        return "infrastructure -> solver"
    return None


def resolve_include(root: Path, source: Path, include: str, headers: dict[str, list[Path]]) -> Path | None:
    source_root = root / "src"
    if "/" in include:
        candidate = source_root / include
        if candidate.is_file():
            return candidate
        suffix_matches = [path for paths in headers.values() for path in paths
                          if path.as_posix().endswith(include)]
        return suffix_matches[0] if len(suffix_matches) == 1 else None
    local = source.parent / include
    if local.is_file():
        return local
    candidates = headers.get(Path(include).name, [])
    return candidates[0] if len(candidates) == 1 else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    source_root = root / "src"
    allowlist_path = root / "tools" / "architecture_allowlist.json"
    baseline = json.loads(allowlist_path.read_text())
    allowed = {(item["source"], item["target"], item["rule"]): item["reason"]
               for item in baseline["dependency_allowlist"]}

    # AppleDouble sidecar files can be created by macOS on external volumes.
    # They are filesystem metadata, not C/C++ translation units or headers.
    files = [path for path in source_root.rglob("*")
             if path.is_file() and not path.name.startswith("._")
             and path.suffix in SOURCE_SUFFIXES]
    headers: dict[str, list[Path]] = defaultdict(list)
    for path in files:
        if path.suffix in HEADER_SUFFIXES:
            headers[path.name].append(path)

    duplicate_headers = {name: sorted(relative(root, path) for path in paths)
                         for name, paths in headers.items() if len(paths) > 1}
    potential_bare: dict[str, list[str]] = defaultdict(list)
    known_debt, violations = [], []
    for source in files:
        source_name = relative(root, source)
        for line_number, line in enumerate(source.read_text(errors="replace").splitlines(), 1):
            match = INCLUDE.match(line)
            if not match:
                continue
            include = match.group(1)
            if "/" not in include and len(headers.get(Path(include).name, [])) > 1:
                potential_bare[include].append(f"{source_name}:{line_number}")
            target = resolve_include(root, source, include, headers)
            if target is None or not target.is_relative_to(source_root):
                continue
            target_name = relative(root, target)
            rule = dependency_rule(source_name, target_name)
            if not rule:
                continue
            key = (source_name, target_name, rule)
            if key in allowed:
                known_debt.append((key, allowed[key]))
            else:
                violations.append(key)

    duplicate_errors = []
    baseline_duplicates = baseline["duplicate_baseline"]
    for name, paths in duplicate_headers.items():
        expected = baseline_duplicates.get(name)
        if expected is None:
            duplicate_errors.append(f"new duplicate basename {name}: {paths}")
        elif any(path not in expected for path in paths):
            duplicate_errors.append(f"expanded duplicate basename {name}: {paths}")

    authority_errors = []
    authority_checks = {
        "src/solver/equation/compressible/SF_compressible.cpp": [
            ("definition_.add", "compressible adapter recreates equation definitions"),
            ("System::System()", "compressible adapter has an unbound default authority"),
        ],
        "src/solver/discretization/source/SF_sourceTerm.h": [
            ("switch (kind)", "density source execution uses a central SourceKind switch"),
        ],
        "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp": [
            ("switch (kind)", "Eulerian source composition uses a central SourceKind switch"),
        ],
    }
    for name, forbidden in authority_checks.items():
        content = (root / name).read_text(errors="replace")
        for token, detail in forbidden:
            if token in content:
                authority_errors.append(f"{name}: {detail}")
    resolved_header = (root / "src/solver/system/SF_resolvedSimulationSystem.h").read_text()
    equation_header = (root / "src/core/system/SF_equationIR.h").read_text(errors="replace")
    solve_program_header = (root / "src/core/system/SF_solveProgram.h").read_text(errors="replace")
    numerical_header = (root / "src/solver/system/SF_numericalSystem.h").read_text(errors="replace")
    runtime_header = (root / "src/solver/system/SF_runtimeRequirements.h").read_text(errors="replace")
    if "Equation::System equationDefinitions" not in equation_header:
        authority_errors.append(
            "Executable equation authority does not own equation definitions")
    run_flow = (root / "src/app/application/execution/SF_flowLoop.cpp").read_text()
    if "stepper.bindSolvePlan(plan)" not in run_flow:
        authority_errors.append(
            "Generic flowLoop does not bind the compiled plan to runtime execution")
    stepper_interface = (
        root / "src/core/interfaces/SF_solverStepper.h"
    ).read_text()
    if "bindSolvePlan(const System::CompiledSolvePlan& plan)" not in stepper_interface:
        authority_errors.append(
            "INavierStokesStepper has no stable compiled-plan binding")
    if "virtual void prepare(SolverState& state)" not in stepper_interface:
        authority_errors.append(
            "Transient executor has no pre-timestep state realization hook")
    if "stepper.prepare(state)" not in run_flow:
        authority_errors.append(
            "Generic flowLoop does not realize state before Time::Driver")
    validator = (root / "src/solver/system/SF_systemValidator.cpp").read_text()
    if "AssemblyPlanRegistry plans" not in validator:
        authority_errors.append(
            "System validator does not build AssemblyPlans for every equation")
    state_realizer = root / "src/solver/system/SF_stateRealizer.cpp"
    if not state_realizer.is_file() or "workspaceRequirements" not in state_realizer.read_text():
        authority_errors.append(
            "Resolved unknown/workspace requirements are not realized before run")
    pressure_stepper = (
        root / "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp"
    ).read_text()
    if "writeCanonicalFace(prefix" in pressure_stepper:
        authority_errors.append(
            "Eulerian face workspace still uses implicit finalize instead of borrowed storage")
    plan_executor = (root / "src/solver/run/SF_planExecutor.cpp").read_text()
    plan_executor_header = (
        root / "src/solver/run/SF_planExecutor.h"
    ).read_text(errors="replace")
    if "PlanNodeKind::StageLoop" not in plan_executor \
            or "execution.stageIndex = stage" not in plan_executor:
        authority_errors.append(
            "PlanExecutor does not implement StageLoop execution context")
    for token in ["Euler", "SSPRK3", "RK4", "density", "IBM"]:
        if token in plan_executor:
            authority_errors.append(
                "PlanExecutor inspects explicit/solver-specific token '"
                + token + "'")
    if "executeOperation" in plan_executor \
            or "executeOperation" in plan_executor_header:
        authority_errors.append(
            "PlanExecutor transitional executeOperation entry still exists")
    if "ExecutionBackendKind" in plan_executor or "GenericPisoPlan" in plan_executor:
        authority_errors.append(
            "PlanExecutor contains solver-family runtime dispatch")
    if "using OpId = std::string" not in solve_program_header:
        authority_errors.append(
            "Compiled solve plan still uses a closed operation enum")
    if "legacyBlocks" in resolved_header or "ExecutionBackendKind" in resolved_header:
        authority_errors.append(
            "Resolved plan still contains legacy block/backend authority")
    if "SolveStage" in stepper_interface or "SolveStageKind" in stepper_interface:
        authority_errors.append(
            "Dead Workflow-era SolveStage types remain in core interfaces")
    for token in ["Eulerian", "IBM", "turbulence", "SolverAlgorithm"]:
        if token in plan_executor:
            authority_errors.append(
                "PlanExecutor contains physics/solver-specific dispatch token '"
                + token + "'")

    solve_plan_source = (
        root / "src/solver/system/SF_solvePlan.cpp"
    ).read_text(errors="replace")
    # Numerical providers are declared by formulation/model contributors;
    # the generic planner only lowers fragment/control-flow values.
    if any(token in solve_plan_source for token in
           ["compileEulerianPimple", "PIMPLE", "PISO", "SIMPLE",
            "ee.pressure", "ee.momentum"]):
        authority_errors.append("SolvePlanner still dispatches a pressure/Eulerian algorithm")
    if "requireAllPoliciesConsumed" not in solve_plan_source:
        authority_errors.append(
            "SolvePlanner has no consume-or-fail validation for active policies")
    if 'operations.bind("ee.pressure.publish",[&] {})' in pressure_stepper:
        authority_errors.append(
            "Eulerian pressure publish is still represented by an empty callback")
    advance_start = pressure_stepper.find("FDM::StepResult EulerianStepper::advance")
    register_start = pressure_stepper.find("void EulerianStepper::registerOperations")
    advance_body = pressure_stepper[advance_start:register_start]
    for token in ["stableTimeStep(", "state_->time+=", "++state_->step"]:
        if token in advance_body:
            authority_errors.append(
                "EulerianStepper::advance still owns timestep lifecycle token '"
                + token + "'")
    builder_source = (
        root / "src/solver/system/SF_systemBuilder.cpp"
    ).read_text(errors="replace")
    provider_resolver = (
        root / "src/solver/system/SF_providerResolver.cpp"
    ).read_text(errors="replace")
    if "SolvePlanner::requiredOperations(plan)" not in provider_resolver:
        authority_errors.append(
            "Provider resolver does not derive operations from CompiledSolvePlan")
    if "phase-wise IBM fluid-port assembly is unavailable" not in builder_source:
        authority_errors.append(
            "Eulerian variational IBM has no explicit unsupported diagnostic")
    single_fluid_preset = (
        root / "src/solver/system/SF_singleFluidPreset.cpp"
    )
    preset_text = single_fluid_preset.read_text(errors="replace") \
        if single_fluid_preset.is_file() else ""
    if not single_fluid_preset.is_file() \
            or "builtin.singleFluidNavierStokes" not in preset_text:
        authority_errors.append(
            "density single-fluid core equations have no authoritative preset")
    for token in ['Equation::ddt({"rhoU"})']:
        if token in builder_source:
            authority_errors.append(
                "SF_systemBuilder.cpp recreates density core equation DSL: "
                + token)
    if "addDensityBasedFluid" in builder_source:
        authority_errors.append(
            "legacy addDensityBasedFluid equation authority remains")
    if "bool fluidDiffusion" in resolved_header:
        authority_errors.append(
            "BuildRequest still owns the legacy fluidDiffusion feature flag")
    for token in ["struct BuildRequest", "SF_turbulenceSystemContribution.h",
                  "SF_levelSetSystemContribution.h"]:
        if token in resolved_header:
            authority_errors.append(
                "resolved runtime contract still exposes composition-only token '"
                + token + "'")
    for token, owner in [
            ("addSourceContributions", "physical source models"),
            ("addTurbulenceEquations", "turbulence model"),
            ('Equation::ddt({"phi"})', "LevelSet model"),
            ("void addImmersed(", "IBM model")]:
        if token in builder_source:
            authority_errors.append(
                "SF_systemBuilder.cpp still owns " + owner
                + " contribution token '" + token + "'")
    for source in files:
        if source.name == "SF_interfaces.h":
            continue
        content = source.read_text(errors="replace")
        if re.search(r'#\s*include\s*[<\"](?:core/interfaces/)?SF_interfaces\.h[>\"]',
                     content):
            authority_errors.append(
                relative(root,source)
                + ": production source includes the compatibility interface umbrella")
    runtime_contract = (
        root / "src/core/interfaces/SF_executionRuntime.h"
    ).read_text(errors="replace")
    if "LocalExecutionRuntime" in runtime_contract or "LocalRuntime final" in runtime_contract:
        authority_errors.append(
            "core execution contract still defines a concrete local runtime")

    # Composition -> Transformation -> Planning authority guards.
    if "struct CompiledNumericalSystem" not in numerical_header:
        authority_errors.append("CompiledNumericalSystem type is missing")
    for token, detail, header in [
        ("struct RawEquationSystem", "RawEquationSystem type is missing", equation_header),
        ("struct ExecutableEquationSystem", "ExecutableEquationSystem type is missing", equation_header),
        ("struct CompiledSolvePlan", "CompiledSolvePlan type is missing", solve_program_header),
        ("class IEquationSystemTransformer",
         "IEquationSystemTransformer contract is missing",
         (root / "src/solver/system/SF_transformation.h").read_text(errors="replace")),
    ]:
        if token not in header:
            authority_errors.append(detail)
    if "RawEquationSystem rawSystem" not in resolved_header:
        authority_errors.append("Resolved system has no raw composition authority")
    if "ExecutableEquationSystem executableSystem" not in resolved_header:
        authority_errors.append("Resolved system has no executable equation authority")
    if "CompiledSolvePlan solvePlan" not in resolved_header:
        authority_errors.append("Resolved system has no compiled solve-plan authority")

    builder_source = (
        root / "src/solver/system/SF_systemBuilder.cpp"
    ).read_text(errors="replace")
    raw_pressure_slice = builder_source[
        builder_source.find("void addPressureBasedFluid("):
        builder_source.find("void addPhaseEquationPack(")]
    if 'addUnknown(system,"pPrime"' in raw_pressure_slice:
        authority_errors.append(
            "pPrime is still composed as a raw pressure-system unknown")
    if 'addEquation(system,{"E_PRESSURE"' in raw_pressure_slice:
        authority_errors.append(
            "Pressure correction is still composed as a builtin raw equation")

    transformation = (
        root / "src/solver/system/SF_transformation.cpp"
    ).read_text(errors="replace")
    if "MPI_" in transformation or "mpi.h" in transformation:
        authority_errors.append(
            "Equation-system transformer directly accesses MPI implementation")
    if "SolverAlgorithm::DensityBased" in transformation \
            or "SolverAlgorithm::PressureBased" in transformation:
        authority_errors.append(
            "Transformer applicability depends on density/pressure solver identity")
    if "StorageBinding::TransientWorkspace" not in transformation:
        authority_errors.append(
            "Generated pressure correction is not declared as transient workspace")
    solve_planner = (root / "src/solver/system/SF_solvePlan.cpp").read_text(
        errors="replace")
    if "struct ExecutionCapabilitySignature" not in runtime_header \
            or "compileExecutionCapabilities(" not in (root / "src/solver/system/SF_providerResolver.cpp").read_text(errors="replace"):
        authority_errors.append(
            "Generic pressure routing has no execution capability signature")
    if "isMinimalGenericPiso" in solve_planner:
        authority_errors.append(
            "Generic pressure routing still depends on an exact equation-id set")
    generic_executor = (
        root / "src/solver/run/SF_planExecutor.cpp"
    ).read_text(errors="replace")
    for token, detail in [
        ("LegacyPressureExecutionAdapter",
         "Generic PISO executor calls the legacy pressure adapter"),
        ("turbulence", "Generic PISO executor queries turbulence identity"),
        ("phase", "Generic PISO executor queries phase identity"),
        ("IBM", "Generic PISO executor queries IBM identity"),
        ("AssembleMomentumPredictor",
         "Generic executor duplicates PISO operation order"),
        ("SolvePressureCorrection",
         "Generic executor duplicates PISO operation order"),
    ]:
        if token in generic_executor:
            authority_errors.append(detail)
    if "is not registered" not in transformation \
            or "cannot apply" not in transformation:
        authority_errors.append(
            "Unknown/inapplicable transformations are not fail-visible")

    execution_sources = "\n".join(
        path.read_text(errors="replace")
        for path in (root / "src/app/application/execution").glob("*.cpp")
    )
    if "templateOrigin" in execution_sources or "PhysicsTemplateKind" in execution_sources:
        authority_errors.append(
            "Application execution dispatches lifecycle from physics template identity")

    # Application execution/ is orchestration-only. Concrete composition,
    # adapters, output field views, and equation coupling have explicit owners.
    execution_dir = root / "src/app/application/execution"
    allowed_execution_files = {
        "SF_execution.h", "SF_execution.cpp",
        "SF_flowLoop.h", "SF_flowLoop.cpp",
        "SF_singleFluid.cpp", "SF_multiPatch.cpp", "SF_eulerian.cpp",
    }
    unexpected_execution_files = sorted(
        path.name for path in execution_dir.iterdir()
        if path.is_file() and not path.name.startswith("._")
        and path.suffix in SOURCE_SUFFIXES
        and path.name not in allowed_execution_files)
    if unexpected_execution_files:
        authority_errors.append(
            "application/execution contains non-orchestration files: "
            + ", ".join(unexpected_execution_files))

    obsolete_paths = [
        "src/solver/algorithm/SF_densityBasedTime.h",
        "src/solver/algorithm/SF_densityBasedTime.cpp",
        "src/core/interfaces/SF_flowAlgorithm.h",
        "src/solver/algorithm/SF_solverAlgorithm.h",
        "src/solver/algorithm/SF_solverAlgorithm.cpp",
        "src/solver/algorithm/densityBased/SF_correction.h",
        "src/solver/algorithm/densityBased/SF_correction.cpp",
        "src/solver/algorithm/pressureBased/SF_adapter.h",
        "src/solver/algorithm/pressureBased/SF_adapter.cpp",
        "src/app/application/run/SF_equationCoupling.h",
        "src/app/application/run/SF_equationCoupling.cpp",
        "src/app/application/run/SF_interfaceCoupling.h",
        "src/app/application/run/SF_interfaceCoupling.cpp",
        "src/app/application/run/SF_services.h",
        "src/app/application/run/SF_services.cpp",
        "src/app/application/run/SF_runtime.h",
        "src/app/application/run/SF_runtime.cpp",
        "src/app/application/run/SF_output.h",
        "src/app/application/run/SF_output.cpp",
        "src/app/application/run/SF_run.h",
        "src/app/application/run/SF_run.cpp",
        "src/app/application/run/SF_runFlow.h",
        "src/app/application/run/SF_runFlow.cpp",
        "src/app/application/SF_preflight.h",
        "src/app/application/SF_preflight.cpp",
        "src/app/application/execution/SF_executionBuilder.cpp",
        "src/app/application/execution/SF_executionAssemblers.h",
        "src/app/application/execution/SF_runners.h",
        "src/app/application/system/SF_inspection.h",
        "src/app/application/system/SF_inspection.cpp",
    ]
    for obsolete in obsolete_paths:
        if (root / obsolete).exists():
            authority_errors.append(
                "obsolete responsibility path still exists: " + obsolete)

    explicit_time = root / "src/solver/algorithm/time/SF_explicit.cpp"
    explicit_header = root / "src/solver/algorithm/time/SF_explicit.h"
    if not explicit_time.is_file():
        authority_errors.append("generic explicit time integrator is missing")
    else:
        explicit_source = explicit_time.read_text(errors="replace")
        explicit_contract = explicit_header.read_text(errors="replace") \
            if explicit_header.is_file() else ""
        if "app/application/run" in explicit_source:
            authority_errors.append(
                "explicit time integration depends on application/run")
        if "void advance(" in explicit_source \
                or "void advance(" in explicit_contract:
            authority_errors.append(
                "Time::Explicit still owns a full-step advance entry")
        if "for (int stage" in explicit_source:
            authority_errors.append(
                "Time::Explicit still owns a global stage loop")

    # Time recipe authority and mathematical Equation::Term purity.
    expression_header = (
        root / "src/core/system/SF_expression.h"
    ).read_text(errors="replace")
    for token in ["EvaluationMode", "TimeLevel", "ExplicitSource",
                  "ImplicitSource", "implicitSource", "TimeRecipe",
                  "TimeScheme", "StageLoop", "Newton", "HYPRE"]:
        if token in expression_header:
            authority_errors.append(
                "Equation::Term contains execution authority token '"
                + token + "'")

    physics_slice = builder_source[
        builder_source.find("void addDensityBasedFluid"):
        builder_source.find("struct SourceEquationContribution")]
    for token in ["TimeRecipe", "TimeScheme", "timeIntegrator",
                  "ExecutionPolicyKind::ExplicitStages"]:
        if token in physics_slice:
            authority_errors.append(
                "physics composition contains time authority token '"
                + token + "'")

    for token in ["parseTimeScheme", "resolveTimeRecipe", "explicitStageCount",
                  '"Euler"', '"RK4"', '"SSPRK3"',
                  "TimeRecipeId::ForwardEuler", "TimeRecipeId::SSPRK3",
                  "TimeRecipeId::ClassicalRK4"]:
        if token in solve_plan_source:
            authority_errors.append(
                "SolvePlanner contains concrete/parsed time authority token '"
                + token + "'")

    all_source_text = "\n".join(
        path.read_text(errors="replace") for path in files)
    for token in ["explicitStageCount(", "ddtDispatch(",
                  "class ExplicitSolver", "class ImplicitSolver",
                  "class SemiImplicitSolver", "class RK4Solver",
                  "class SDIRKSolver"]:
        if token in all_source_text:
            authority_errors.append(
                "obsolete/global time authority symbol remains: " + token)
    for obsolete in [
        "src/solver/discretization/time/SF_time.h",
        "src/methods/numerics/time/SF_time.h",
        "src/methods/numerics/time/SF_Euler.h",
        "src/methods/numerics/time/SF_Euler.cpp",
        "src/methods/numerics/time/SF_RK4.h",
        "src/methods/numerics/time/SF_RK4.cpp",
    ]:
        if (root / obsolete).exists():
            authority_errors.append(
                "obsolete time-dispatch path still exists: " + obsolete)

    # Mathematical terms and compiled recipes are the only spatial-method
    # authority beyond the IO compatibility boundary.
    if "requiredGhostLayersForConvection" in all_source_text:
        authority_errors.append(
            "scheme-name halo authority remains: requiredGhostLayersForConvection")
    numerical_compiler = (
        root / "src/solver/system/SF_numericalCompiler.cpp"
    ).read_text(errors="replace")
    if "equations.equations" not in numerical_compiler \
            or "result.terms.emplace_back" not in numerical_compiler:
        authority_errors.append(
            "NumericalCompiler does not enumerate and bind executable terms")
    if "result.numericalSystem = NumericalCompiler::compile" \
            not in builder_source:
        authority_errors.append(
            "Resolved system does not compile the production numerical system")
    compressible_equation = (
        root / "src/solver/equation/compressible/SF_compressible.cpp"
    ).read_text(errors="replace")
    for token in ["viscousEnabled", "numerics.convection",
                  "numerics.flux", "parseConvectionScheme",
                  "parseFluxSplitter"]:
        if token in compressible_equation:
            authority_errors.append(
                "compressible equation execution reads raw numerical authority '"
                + token + "'")
    source_provider = (
        root / "src/solver/discretization/source/SF_sourceTerm.h"
    ).read_text(errors="replace")
    if "config.enabled" in source_provider:
        authority_errors.append(
            "source provider reselects mathematical source presence")
    solve_plan_text = solve_plan_source.lower()
    for token in ["weno", "teno", "rusanov", "stegerwarming"]:
        if token in solve_plan_text:
            authority_errors.append(
                "SolvePlanner contains spatial recipe routing token '"
                + token + "'")

    compressible_source = (
        root / "src/solver/algorithm/SF_singleFluidStepper.cpp"
    ).read_text(errors="replace")
    advance_start = compressible_source.find(
        "FDM::StepResult SingleFluidStepper::advance")
    advance_body = compressible_source[advance_start:]
    for token in ["switch (config_.numerics.time", "for (int stage",
                  "Time::Explicit::advance",
                  "SolverAlgorithm::DensityBased",
                  "SolverAlgorithm::PressureBased"]:
        if token in advance_body:
            authority_errors.append(
                "SingleFluidStepper::advance contains forbidden runtime "
                "routing token '" + token + "'")
    flow_authority_text = "\n".join(
        path.read_text(errors="replace") for path in files)
    for token in ["IFlowAlgorithm", "FlowAlgorithmContext",
                  "FlowAlgorithmResult", "makeFlowAlgorithm",
                  "flowAlgorithm_", "DensityBased::Algorithm",
                  "PressureBased::Algorithm",
                  "SingleFluidStepper::stepPressure",
                  "LegacyAdapterRequired", "requiredAdapters"]:
        if token in flow_authority_text:
            authority_errors.append(
                "removed flow-algorithm authority symbol remains: " + token)
    if "PlanExecutor::validateBindings" not in compressible_source:
        authority_errors.append(
            "single-fluid runtime does not validate Plan operation coverage")
    if "validateBindings" not in plan_executor:
        authority_errors.append(
            "PlanExecutor has no pre-execution operation coverage validation")
    solve_plan_source = (
        root / "src/solver/system/SF_solvePlan.cpp"
    ).read_text(errors="replace")
    for token in ['strategyName == "PISO"',
                  'strategyName == "SIMPLE"',
                  'strategyName == "PIMPLE"']:
        if token in solve_plan_source:
            authority_errors.append(
                "pressure Plan lowering dispatches on display name: " + token)
    if ('"pressure.schedule."+node.id' not in (
            root / "src/solver/system/SF_pressureCoupling.cpp"
        ).read_text(errors="replace")):
        authority_errors.append(
            "unsupported pressure fixed-point schedule has no dedicated "
            "predictor operation")
    validator_source = (
        root / "src/solver/system/SF_systemValidator.cpp"
    ).read_text(errors="replace")
    if "has no explicit runtime operation ID" not in validator_source:
        authority_errors.append(
            "compiled Plan leaves are not required to carry an OpId")
    # 压力 schedule 映射的唯一 authority 是 coupling contribution。
    schedule_source = (
        root / "src/solver/system/SF_pressureCoupling.cpp"
    ).read_text(errors="replace")
    for token in ["correction.repeatCount = request.outerCorrectors",
                  "correction.nestedRepeatCount = request.pressureCorrectors",
                  "correction.innerRepeatCount = request.nonOrthogonalCorrectors+1"]:
        if token not in schedule_source:
            authority_errors.append(
                "pressure policy lost canonical schedule mapping: " + token)
    system_builder_source = (
        root / "src/solver/system/SF_systemBuilder.cpp"
    ).read_text(errors="replace")
    for token in ["S_PRESSURE", "S_EE_PIMPLE"]:
        if token in system_builder_source:
            authority_errors.append(
                "system builder still fabricates coupling schedule '" + token
                + "' instead of consuming the coupling contribution")
    coupling_dir = root / "src/solver/equation/coupling"
    for path in coupling_dir.glob("*.*"):
        content = path.read_text(errors="replace")
        if "app/application/execution" in content \
                or "app/application/run" in content:
            authority_errors.append(
                "equation coupling depends on application execution: "
                + relative(root, path))

    # Execution composition consumes compiler-derived provider requirements.
    # It must not rediscover physics by inspecting equation IDs, strategy
    # labels, or solver-family enums.
    execution_dir = root / "src/app/application/execution"
    execution_sources = "\n".join(
        path.read_text(errors="replace")
        for path in execution_dir.glob("*.cpp")
    )
    for token, detail in [
        ("System::hasEquation(",
         "application execution selects providers from equation IDs"),
        ("System::hasSolveStrategy(",
         "application execution selects a runner from solve strategy"),
    ]:
        if token in execution_sources:
            authority_errors.append(detail)
    execution_builder = (
        execution_dir / "SF_execution.cpp"
    ).read_text(errors="replace")
    for token in ["PressureVelocityCoupling", "SolverAlgorithm::",
                  "PhysicsTemplateKind"]:
        if token in execution_builder:
            authority_errors.append(
                "execution builder contains solver-family routing token '"
                + token + "'")
    environment_source = (
        root / "src/app/application/SF_environment.cpp"
    ).read_text(errors="replace")
    if "solverConfig.numerics.solver" in environment_source:
        authority_errors.append(
            "execution environment infers runtime services from solver identity")
    for token in ["NoopImmersedBoundary", "NoopBoundaryPipeline",
                  "LocalParallelCoordinator", "NoopTransportModel"]:
        if token in flow_authority_text:
            authority_errors.append(
                "zero-caller compatibility type remains: " + token)
    for token in ["LegacyMultiphaseEquationCoupling",
                  "MultiPatchLegacyEquationCoupling",
                  "registerLegacyState", "InterfaceEquationCoupling"]:
        if token in flow_authority_text:
            authority_errors.append(
                "removed legacy/wiring provider symbol remains: " + token)
    # RUNTIME authority is RuntimeRequirements (SF_runtimeRequirements.h).
    # The aggregate must hold exactly one such member and must not re-declare
    # the individual requirement vectors as loose duplicate fields.
    if "providerRequirements" not in runtime_header \
            or "runtimeServiceRequirements" not in runtime_header:
        authority_errors.append(
            "runtime authority has no explicit provider/runtime-service requirements")
    if "RuntimeRequirements runtime" not in resolved_header:
        authority_errors.append(
            "resolved system does not own a single RuntimeRequirements aggregate")
    for token in ["providerRequirements", "runtimeServiceRequirements",
                  "workspaceRequirements"]:
        if token in resolved_header:
            authority_errors.append(
                "resolved system duplicates runtime requirement field '" + token
                + "' instead of delegating to RuntimeRequirements")
    if "deriveExecutionComposition(result,request)" not in builder_source:
        authority_errors.append(
            "execution-provider requirements are not compiler-derived")

    source_text = "\n".join(
        path.read_text(errors="replace") for path in files)
    if "DensityBasedTime" in source_text or "SF_densityBasedTime" in source_text:
        authority_errors.append(
            "DensityBasedTime naming still owns explicit integration")

    composition = (
        root / "src/solver/system/SF_equationContribution.cpp"
    ).read_text(errors="replace")
    if "SystemModification::Add must use a typed add operation" not in composition \
            or "lowering is not implemented" not in composition:
        authority_errors.append(
            "Incomplete user modifications can be silently ignored")

    # §20：solver-family 命名与 foreign cpp ownership 必须保持消失。
    production_text = "\n".join(
        path.read_text(errors="replace")
        for path in list((root / "src").rglob("*.cpp"))
        + list((root / "src").rglob("*.h")))
    for token in ["CompressibleAlgorithm", "DensityBasedRHS",
                  "DensityBasedSolver", "PressureBasedSolver",
                  "runSimpleSolver", "runPisoSolver",
                  "PisoStepper", "SimpleStepper"]:
        if token in production_text:
            authority_errors.append(
                "removed solver-family symbol remains in production source: "
                + token)
    for cmake in (root / "src").rglob("CMakeLists.txt"):
        text = cmake.read_text(errors="replace")
        if cmake.parent.name == "system":
            for line in text.splitlines():
                stripped = line.strip()
                if stripped.endswith(".cpp") and ".." in stripped:
                    authority_errors.append(
                        "solver/system CMake compiles foreign implementation: "
                        + stripped)
    # 压耦合 preset 必须由 resolved equation/constraint structure 匹配，而不是
    # 由 density/pressure 标签或 solver enum 选择。
    coupling_source = (
        root / "src/solver/system/SF_pressureCoupling.cpp"
    ).read_text(errors="replace")
    for token in ["C_INCOMPRESSIBILITY", "C_SHARED_PRESSURE"]:
        if token not in coupling_source:
            authority_errors.append(
                "pressure-coupling matching ignores constraint '" + token + "'")
    if "SolverAlgorithm" in coupling_source:
        authority_errors.append(
            "pressure-coupling matching dispatches on solver family enum")
    resolved_text = (
        root / "src/solver/system/SF_resolvedSimulationSystem.h"
    ).read_text(errors="replace")
    if "FDM::SolverAlgorithm" in resolved_text:
        authority_errors.append(
            "resolved system still carries a solver-family formulation field")

    # §P25：executable operation authority 必须唯一。
    # 1) coupling layer 不得自带 operation/provider 名单。
    coupling_text = (
        root / "src/solver/system/SF_pressureCoupling.cpp"
    ).read_text(errors="replace")
    for token in ["couplingOperations", "couplingProviders",
                  '"pressure.assemble"', '"pressure.solve"',
                  '"pressure.prepare"', '"velocity.correct"',
                  '"flux.correct"', '"momentum.solve"', '"momentum.assemble"']:
        if token in coupling_text:
            authority_errors.append(
                "coupling layer duplicates operation/provider authority: "
                + token)
    # 2) SolvePlanner 不得知道 pressure 数学或 provider 名单。
    planner_text = (
        root / "src/solver/system/SF_solvePlan.cpp"
    ).read_text(errors="replace")
    for token in ["compileGenericPiso", "compileUnavailablePressureSchedule",
                  "compileEulerianPimple", "PIMPLE", "PISO", "SIMPLE",
                  "ee.pressure", "ee.momentum",
                  '"pressure.assemble"', '"pressure.solve"',
                  '"pressure.prepare"', '"pressure.boundary.prepare"',
                  '"pressure.update.prepare"', '"velocity.correct"',
                  '"flux.correct"', '"pressure.correction.commit"',
                  '"pressure.step.commit"']:
        if token in planner_text:
            authority_errors.append(
                "SolvePlanner still owns pressure operation semantics: "
                + token)
    # 3) formulation 必须声明 executable operations（唯一 authority）。
    transformation_text = (
        root / "src/solver/system/SF_transformation.cpp"
    ).read_text(errors="replace")
    for token in ["addExecutableOperation", "OperationStage::PressureSolve",
                  "OperationStage::VelocityCorrect"]:
        if token not in transformation_text:
            authority_errors.append(
                "pressure formulation does not declare executable operation: "
                + token)

    # §P25B: descriptions stay provider-neutral, models consume only core IR,
    # and unused recipes cannot be exempted by a system-wide constraint flag.
    operation_slice = equation_header[
        equation_header.find("struct ExecutableOperation {"):
        equation_header.find("enum class ResourceAccessMode")]
    if "std::string provider" in operation_slice:
        authority_errors.append(
            "ExecutableOperation still chooses a concrete provider")
    resolver_source = (
        root / "src/solver/system/SF_providerResolver.cpp"
    ).read_text(errors="replace")
    if "resolveOperationBindings" not in resolver_source \
            or "reportOperationBindings" not in resolver_source:
        authority_errors.append(
            "compiled provider resolution/report authority is missing")
    numerical_compiler = (
        root / "src/solver/system/SF_numericalCompiler.cpp"
    ).read_text(errors="replace")
    if "strictCoreTermCoverage" in numerical_compiler \
            or "recipeBindings" not in numerical_compiler:
        authority_errors.append(
            "numerical recipe validation still uses a coarse system exemption")
    for model_file in (root / "src/models").rglob("*"):
        if model_file.suffix not in {".h", ".hpp", ".cpp"}:
            continue
        if '#include "solver/system/' in model_file.read_text(errors="replace"):
            authority_errors.append(
                f"model depends on solver/system implementation: {model_file}")
    for model_cmake in (root / "src/models").rglob("CMakeLists.txt"):
        if "SF_solverSystem" in model_cmake.read_text(errors="replace"):
            authority_errors.append(
                f"model target links solver/system implementation: {model_cmake}")

    if not args.quiet:
        print(f"Architecture dependency check: {len(files)} source files")
        print(f"Known dependency debt: {len(known_debt)} allowlisted entries")
        for (source, target, rule), reason in known_debt:
            print(f"  ALLOW {rule}: {source} -> {target} ({reason})")
        print(f"Duplicate header basenames: {len(duplicate_headers)}")
        for name in sorted(duplicate_headers):
            print(f"  {name}: {duplicate_headers[name]}")
            print("    exposure: src is exported by SF_headers; use qualified includes across modules")
            if potential_bare.get(name):
                print(f"    bare callsites: {potential_bare[name]}")
        if potential_bare:
            print("Potential duplicate-basename bare includes are reported, not treated as a compiler-resolution failure.")
        for source, target, rule in violations:
            print(f"ERROR new dependency violation ({rule}): {source} -> {target}", file=sys.stderr)
        for error in duplicate_errors:
            print(f"ERROR {error}", file=sys.stderr)
        for error in authority_errors:
            print(f"ERROR executable authority: {error}", file=sys.stderr)

    return 1 if violations or duplicate_errors or authority_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
