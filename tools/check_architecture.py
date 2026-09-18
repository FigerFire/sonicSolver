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
    if "stepper.bindSolveStages(workflow.stages)" not in run_flow:
        authority_errors.append(
            "Generic runFlow does not bind Workflow stages to runtime execution")
    stepper_interface = (root / "src/core/interfaces/SF_interfaces.h").read_text()
    if "bindSolveStages(const std::vector<SolveStage>& stages)" not in stepper_interface:
        authority_errors.append(
            "INavierStokesStepper has no stable executable solve-stage binding")
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
    workflow = (root / "src/solver/algorithm/workflow/SF_workflow.cpp").read_text()
    if "strategy.find(" in workflow or "strategy.rfind(" in workflow:
        authority_errors.append(
            "Workflow still infers execution semantics from strategy text")
    if "SolveStrategyKind strategyKind" not in resolved_header:
        authority_errors.append(
            "Resolved solve blocks do not own a typed solve strategy")
    if "stage.strategy = block.strategyKind" not in workflow:
        authority_errors.append(
            "Workflow does not propagate the typed solve strategy to runtime stages")


    # application/run is an orchestration boundary. Concrete execution
    # composition and numerical services must have their own owners.
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
    if not explicit_time.is_file():
        authority_errors.append("generic explicit time integrator is missing")
    else:
        explicit_source = explicit_time.read_text(errors="replace")
        if "app/application/run" in explicit_source:
            authority_errors.append(
                "explicit time integration depends on application/run")
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
