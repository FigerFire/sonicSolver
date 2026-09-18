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

    files = [path for path in source_root.rglob("*")
             if path.is_file() and path.suffix in SOURCE_SUFFIXES]
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
        "src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp": [
            ("switch (kind)", "Eulerian source composition uses a central SourceKind switch"),
        ],
    }
    for name, forbidden in authority_checks.items():
        content = (root / name).read_text(errors="replace")
        for token, detail in forbidden:
            if token in content:
                authority_errors.append(f"{name}: {detail}")
    resolved_header = (root / "src/solver/system/SF_resolvedSimulationSystem.h").read_text()
    if "Equation::System equationDefinitions" not in resolved_header:
        authority_errors.append(
            "ResolvedSimulationSystem does not own executable equation definitions")
    run_flow = (root / "src/app/application/run/SF_runFlow.cpp").read_text()
    if "stepper.bindSolvePlan(plan)" not in run_flow:
        authority_errors.append(
            "Generic runFlow does not bind the compiled plan to runtime execution")
    stepper_interface = (root / "src/core/interfaces/SF_interfaces.h").read_text()
    if "bindSolvePlan(const System::CompiledSolvePlan& plan)" not in stepper_interface:
        authority_errors.append(
            "INavierStokesStepper has no stable compiled-plan binding")
    if "virtual void prepare(SolverState& state)" not in stepper_interface:
        authority_errors.append(
            "Transient executor has no pre-timestep state realization hook")
    if "stepper.prepare(state)" not in run_flow:
        authority_errors.append(
            "Generic runFlow does not realize state before Time::Driver")
    validator = (root / "src/solver/system/SF_systemValidator.cpp").read_text()
    if "AssemblyPlanRegistry plans" not in validator:
        authority_errors.append(
            "System validator does not build AssemblyPlans for every equation")
    state_realizer = root / "src/solver/system/SF_stateRealizer.cpp"
    if not state_realizer.is_file() or "workspaceRequirements" not in state_realizer.read_text():
        authority_errors.append(
            "Resolved unknown/workspace requirements are not realized before run")
    pressure_stepper = (
        root / "src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp"
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
    if "using OpId = std::string" not in resolved_header:
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
    for operation in ["ee.dt.compute", "ee.time.commit",
                      "ibm.constraint.project", "ibm.kkt.solve"]:
        if operation not in solve_plan_source:
            authority_errors.append(
                "Compiled Plan is missing required operation '" + operation + "'")
    if "requireAllPoliciesConsumed" not in solve_plan_source:
        authority_errors.append(
            "SolvePlanner has no consume-or-fail validation for active policies")
    if 'operations.bind("ee.pressure.publish",[&] {})' in pressure_stepper:
        authority_errors.append(
            "Eulerian pressure publish is still represented by an empty callback")
    advance_start = pressure_stepper.find("FDM::StepResult PressureStepper::advance")
    register_start = pressure_stepper.find("void PressureStepper::registerOperations")
    advance_body = pressure_stepper[advance_start:register_start]
    for token in ["stableTimeStep(", "state_->time+=", "++state_->step"]:
        if token in advance_body:
            authority_errors.append(
                "PressureStepper::advance still owns timestep lifecycle token '"
                + token + "'")
    builder_source = (
        root / "src/solver/system/SF_systemBuilder.cpp"
    ).read_text(errors="replace")
    if "SolvePlanner::requiredOperations(result.solvePlan)" not in builder_source:
        authority_errors.append(
            "Runtime operation requirements are not derived from CompiledSolvePlan")
    if "phase-wise IBM fluid-port assembly is unavailable" not in builder_source:
        authority_errors.append(
            "Eulerian variational IBM has no explicit unsupported diagnostic")

    # Composition -> Transformation -> Planning authority guards.
    for token, detail in [
        ("struct RawEquationSystem", "RawEquationSystem type is missing"),
        ("struct ExecutableEquationSystem", "ExecutableEquationSystem type is missing"),
        ("struct CompiledSolvePlan", "CompiledSolvePlan type is missing"),
        ("class IEquationSystemTransformer",
         "IEquationSystemTransformer contract is missing"),
    ]:
        haystack = resolved_header
        if token == "class IEquationSystemTransformer":
            transform_header = root / "src/solver/system/SF_transformation.h"
            haystack = transform_header.read_text() if transform_header.is_file() else ""
        if token not in haystack:
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
    if "struct ExecutionCapabilitySignature" not in resolved_header \
            or "capabilitySignature(" not in solve_planner:
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

    run_sources = "\n".join(
        path.read_text(errors="replace")
        for path in (root / "src/app/application/run").glob("*.cpp")
    )
    if "templateOrigin" in run_sources or "PhysicsTemplateKind" in run_sources:
        authority_errors.append(
            "Application runner dispatches lifecycle from physics template identity")

    # Application run/ is orchestration-only. Concrete composition, adapters,
    # output field views, and equation coupling have explicit owners.
    run_dir = root / "src/app/application/run"
    allowed_run_files = {
        "SF_run.h", "SF_run.cpp", "SF_runFlow.h", "SF_runFlow.cpp",
    }
    unexpected_run_files = sorted(
        path.name for path in run_dir.iterdir()
        if path.is_file() and path.suffix in SOURCE_SUFFIXES
        and path.name not in allowed_run_files)
    if unexpected_run_files:
        authority_errors.append(
            "application/run contains non-orchestration files: "
            + ", ".join(unexpected_run_files))

    obsolete_paths = [
        "src/solver/algorithm/SF_densityBasedTime.h",
        "src/solver/algorithm/SF_densityBasedTime.cpp",
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

    compressible_source = (
        root / "src/solver/algorithm/SF_compressible.cpp"
    ).read_text(errors="replace")
    density_start = compressible_source.find(
        "void CompressibleAlgorithm::stepDensity")
    advance_start = compressible_source.find(
        "FDM::StepResult CompressibleAlgorithm::advance", density_start)
    density_step = compressible_source[density_start:advance_start]
    for token in ["switch (config_.numerics.time", "for (int stage",
                  "Time::Explicit::advance"]:
        if token in density_step:
            authority_errors.append(
                "CompressibleAlgorithm::stepDensity still owns explicit "
                "stage control token '" + token + "'")
    coupling_dir = root / "src/solver/equation/coupling"
    for path in coupling_dir.glob("*.*"):
        content = path.read_text(errors="replace")
        if "app/application/execution" in content \
                or "app/application/run" in content:
            authority_errors.append(
                "equation coupling depends on application execution: "
                + relative(root, path))

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
