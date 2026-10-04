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
    formula_header = (root / "src/core/system/SF_formula.h").read_text(
        errors="replace")
    formula_source = (root / "src/core/system/SF_formula.cpp").read_text(
        errors="replace")
    for name, content in (("SF_formula.h", formula_header),
                          ("SF_formula.cpp", formula_source)):
        if re.search(r'^\s*#\s*include\s*[<\"](?:mpi|HYPRE|solver/)',
                     content, re.MULTILINE):
            authority_errors.append(
                f"{name}: Formula WHAT includes MPI/HYPRE/solver implementation")
    for token in ("FormulaMode", "solveTarget", "LegacyExecutionPolicyKind"):
        if token in formula_header:
            authority_errors.append(
                f"SF_formula.h: Formula WHAT owns execution concept {token}")
    formula_compiler = (root / "src/solver/system/SF_formulaCompiler.cpp").read_text(
        errors="replace")
    if re.search(r'\b(?:if|switch)\s*\([^)]*(?:PISO|SIMPLE|kEpsilon|Momentum)',
                 formula_compiler):
        authority_errors.append(
            "FormulaCompiler branches on a named preset/equation")
    resolved_header = (root / "src/solver/system/SF_resolvedSimulationSystem.h").read_text()
    equation_header = (root / "src/core/system/SF_equationIR.h").read_text(errors="replace")
    solve_program_header = (root / "src/core/system/SF_solveProgram.h").read_text(errors="replace")
    program_step = re.search(r"struct EquationCall\s*\{([^}]*)\};", solve_program_header)
    if not program_step or "Target target" not in program_step.group(1):
        authority_errors.append("EquationCall must declare its typed Target")
    elif re.search(r"FormulaMode|TimeRecipe|TemporalMethod|EquationMethod|WENO|HYPRE",
                   program_step.group(1)):
        authority_errors.append("EquationCall owns numerical method or mode")
    method_contract = (root / "src/solver/system/SF_methodObjects.h").read_text(
        errors="replace")
    for name in ("ITemporalMethod", "IProvider"):
        if name not in method_contract:
            authority_errors.append("Method-object contract lacks " + name)
    if "CompiledExecutionProgram" not in solve_program_header:
        authority_errors.append("SolvePlan has no compiled method-step result")
    if "ExecutionKind" not in solve_program_header or "ExecutionScope root" not in solve_program_header:
        authority_errors.append("ExecutionProgram lacks structured HOW authority")
    if "workspaceRequires" not in solve_program_header or "backendOperation" not in solve_program_header:
        authority_errors.append("Compiled method step lacks workspace/backend contract")
    method_source = (root / "src/solver/system/SF_methodObjects.cpp").read_text(
        errors="replace")
    if "compileFragment" not in method_source or "hasTemporalRoot=true" not in method_source:
        authority_errors.append("TemporalMethod does not compile a production fragment")
    numerical_header = (root / "src/solver/system/SF_numericalSystem.h").read_text(errors="replace")
    runtime_header = (root / "src/solver/system/SF_runtimeRequirements.h").read_text(errors="replace")
    if "Equation::System legacyDefinitions" not in equation_header:
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
    if "executable.legacyDefinitions.at(equation)" not in validator or "!equationIds.empty()" not in validator:
        authority_errors.append(
            "System validator lacks isolated legacy equation ownership validation")
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

    pressure_ops = "\n".join((root / ("src/solver/algorithm/pressure/SF_pressureOperators." + suffix)).read_text() for suffix in ("h", "cpp"))
    for forbidden in ['MPI_', 'SparseSystem', 'strategyName', 'CouplingStrategyKind']:
        if forbidden in pressure_ops:
            authority_errors.append("Pressure operator owns backend/preset authority: " + forbidden)
    for required in ['GlobalDofSystem', 'globalMinimum', 'globalMaximum(velocityScale)',
                     'globalMaximum(pressureScale)', 'ExchangeKind::CanonicalFaceFlux']:
        if required not in pressure_ops:
            authority_errors.append("Distributed pressure contract missing: " + required)
    for forbidden in ['SourceConfig', 'GravityConfig', 'MRFConfig',
                      'TurbulenceManager']:
        if forbidden in pressure_ops:
            authority_errors.append(
                "Pressure operator reads uncompiled model configuration: " + forbidden)
    system_root = root / "src/solver/system"
    for source in system_root.rglob("*"):
        if source.is_file() and source.suffix in SOURCE_SUFFIXES \
                and not source.name.startswith("._"):
            if re.search(r'^\s*#\s*include\s*[<\"]models/',
                         source.read_text(errors="replace"), re.MULTILINE):
                authority_errors.append(
                    relative(root, source) + ": solver/system includes a concrete model")
    model_catalog = (system_root / "SF_termProviderCatalog.cpp").read_text()
    for forbidden in ['gravity', 'mrf', 'wallHeat', 'turbulence']:
        if re.search(r'\b' + re.escape(forbidden) + r'\b',
                     model_catalog, re.IGNORECASE):
            authority_errors.append(
                "Generic TermProviderCatalog enumerates model " + forbidden)
    source_executor = (
        root / "src/solver/discretization/source/SF_sourceTerm.h"
    ).read_text()
    for forbidden in ['SourceKind::Gravity', 'SourceKind::MRF',
                      'SourceKind::WallHeat', 'SourceConfig']:
        if forbidden in source_executor:
            authority_errors.append(
                "Generic conservative source executor enumerates model " + forbidden)
    system_cmake = (system_root / "CMakeLists.txt").read_text()
    for forbidden in ['    SF_physics\n', '    SF_ibm\n', '    SF_turbulence\n']:
        if forbidden in system_cmake:
            authority_errors.append(
                "solver/system target links concrete model target " + forbidden.strip())
    for forbidden in ['SourceContribution::contribute',
                      'Turbulence::contribute',
                      'IBM::SystemContribution::contribute',
                      'LevelSetContribution::contribute']:
        if forbidden in (system_root / "SF_systemBuilder.cpp").read_text():
            authority_errors.append(
                "SystemBuilder calls concrete model contributor " + forbidden)
    if 'FDM::SourceConfig sources' in numerical_header:
        authority_errors.append(
            "PressureNumericalConfig retains raw SourceConfig")
    for name in [
            "src/solver/system/SF_numericalCompiler.cpp",
            "src/solver/system/SF_providerResolver.cpp",
            "src/solver/algorithm/SF_singleFluidStepper.cpp",
            "src/solver/algorithm/pressure/SF_pressureOperators.cpp"]:
        content = (root / name).read_text(errors="replace")
        if re.search(r'\b(?:if|switch)\s*\([^\n]*\b(?:gravity|MRF|wallHeat)\b', content):
            authority_errors.append(name + ": model-name source dispatch returned")
        if re.search(r'(?:rfind|starts_with)\s*\(\s*"E_MOMENTUM', content):
            authority_errors.append(name + ": equation role inferred from ID prefix")
    if "RawEquationSystem& raw_" in (
            root / "src/core/system/SF_systemContribution.h").read_text():
        authority_errors.append("SystemContribution retains raw-system lifetime authority")
    builtin_operations = re.findall(
        r'inline constexpr const char\* \w+ = "([^"]+)";',
        (root / "src/core/system/SF_operationIds.h").read_text())
    for name in [
            "src/solver/system/SF_transformation.cpp",
            "src/solver/system/SF_pressureCoupling.cpp",
            "src/solver/system/SF_solvePlan.cpp",
            "src/solver/system/SF_systemBuilder.cpp",
            "src/solver/algorithm/SF_singleFluidStepper.cpp",
            "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp"]:
        content = (root / name).read_text(errors="replace")
        literals = [operation for operation in builtin_operations
                    if '"' + operation + '"' in content]
        if literals:
            authority_errors.append(
                name + ": built-in OpId literals bypass constants: "
                + ", ".join(sorted(literals)))
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
    if "phase-wise IBM fluid-port assembly is unavailable" not in (
            root / "src/solver/system/SF_immersedMethods.cpp").read_text():
        authority_errors.append(
            "Compiled IBM capability has no explicit Eulerian fluid-port unsupported diagnostic")
    immersed_contribution = (root / "src/models/ibm/SF_ibmSystemContribution.cpp").read_text()
    for token in ["addLegacyExecution", "LegacyExecutionPolicy", "EquationDescriptor",
                  "Equation::Definition", "C_IBM", "strategy"]:
        if token in immersed_contribution:
            authority_errors.append("IBM contribution retains legacy authority: " + token)
    if "ImmersedConstraintTransformer" in (
            root / "src/solver/system/SF_transformation.cpp").read_text():
        authority_errors.append("Empty IBM constraint transformer has returned")
    for token in ["C_IBM", "E_IBM", "ImmersedSolveBlockDescriptor", "strategy("]:
        if token in (root / "src/models/ibm/descriptor/SF_algorithmDescriptor.cpp").read_text():
            authority_errors.append("IBM port recreates a duplicate mathematical/schedule inventory: " + token)
    for token in ["IBM", "ibm.", "Immersed"]:
        if token in (root / "src/solver/system/SF_methodObjects.cpp").read_text():
            authority_errors.append("Generic method compiler dispatches an IBM implementation: " + token)
    if "native turbulence + IBM requires" not in provider_resolver:
        authority_errors.append("IBM/turbulence lacks an explicit compiled coupling capability guard")
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
    if re.search(r'addState\s*\([^;]*pPrime', transformation):
        authority_errors.append("Pressure correction was registered as a base STATE symbol")
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
    pressure_leaf_sources = [
        root / "src/solver/algorithm/SF_singleFluidStepper.cpp",
        root / "src/solver/algorithm/pressure/SF_pressureOperators.cpp",
        root / "src/solver/algorithm/pressure/SF_pressureOperators.h",
        root / "src/solver/algorithm/pressure/SF_fixedTimeMath.h",
        root / "src/solver/discretization/pressure/SF_rhieChow.h",
        root / "src/solver/run/SF_planExecutor.cpp",
    ]
    forbidden_family = re.compile(
        r"\b(?:PisoSolver|SimpleSolver|PimpleSolver|"
        r"PisoStepper|SimpleStepper|PimpleStepper|"
        r"PressureBasedSolver|ConstantDensityPisoProvider|"
        r"ConstantDensitySimpleProvider|ConstantDensityPimpleProvider)\b")
    forbidden_branch = re.compile(
        r"\bif\s*\([^)]*\b(?:PISO|SIMPLE|PIMPLE)\b")
    for source in pressure_leaf_sources:
        content = source.read_text(errors="replace")
        if (forbidden_family.search(content) or forbidden_branch.search(content)
                or "PressureCouplingPreset::" in content):
            authority_errors.append(
                relative(root, source)
                + ": pressure numerical leaf reintroduces coupling-family dispatch")
        if source.name in {"SF_pressureOperators.cpp", "SF_rhieChow.h"} \
                and re.search(r"\bMPI_|[<\"]mpi\.h[>\"]", content):
            authority_errors.append(
                relative(root, source)
                + ": pressure numerical kernel directly accesses MPI")
    pressure_provider = (
        root / "src/solver/algorithm/pressure/SF_pressureOperators.h"
    ).read_text(errors="replace")
    if re.search(r"\bSolverConfig\s*[&*]", pressure_provider):
        authority_errors.append(
            "pressure numerical provider consumes the whole SolverConfig")
    if "energyFromPressure" in (
        root / "src/core/config/SF_configTypes.h"
    ).read_text(errors="replace"):
        authority_errors.append(
            "typed multiplier pressure boundary retains legacy energy name")
    stepper_text = (
        root / "src/solver/algorithm/SF_singleFluidStepper.cpp"
    ).read_text(errors="replace")
    if re.search(r"binding\s*->\s*provider\s*==\s*\"flow\.pressure-operators\"",
                 stepper_text):
        authority_errors.append(
            "SingleFluidStepper constructs pressure providers by name")
    binder_header = (
        root / "src/solver/algorithm/SF_pressureProviderBinding.h"
    ).read_text(errors="replace")
    binder_source = (
        root / "src/solver/algorithm/SF_pressureProviderBinding.cpp"
    ).read_text(errors="replace")
    if "SolverConfig" in binder_header or "SolverConfig" in binder_source:
        authority_errors.append(
            "pressure provider binder consumes raw SolverConfig")
    coupling_plan_text = (
        root / "src/solver/system/SF_pressureCoupling.cpp"
    ).read_text(errors="replace")
    if '"PISO.outerCorrectors"' in coupling_plan_text:
        authority_errors.append("PISO retains a fixed-point outer loop")
    if "terminationSignal" not in generic_executor \
            or "LoopSignalState" not in generic_executor \
            or "pressure" in generic_executor.lower():
        authority_errors.append(
            "PlanExecutor does not have pressure-neutral loop termination")
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
                  "LegacyExecutionPolicyKind::ExplicitStages"]:
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
    if "selectedFormulas=selectFormulas(program)" not in numerical_compiler \
            or "result.operators.emplace_back" not in numerical_compiler:
        authority_errors.append(
            "NumericalCompiler does not bind selected Formula occurrences")
    if "NumericalCompiler::compileSystem(" \
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
    if "fusedTargets" in compressible_source or re.search(
            r'call\.equation\s*==\s*"E_(?:MASS|MOMENTUM|ENERGY)"',
            compressible_source):
        authority_errors.append(
            "fused backend selects Formula identity instead of consuming compiled HOW")
    if "FormulaMode::" in compressible_source:
        authority_errors.append("production flow stepper dispatches on compatibility FormulaMode")
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
    pressure_schedule = (
        root / "src/solver/system/SF_pressureCoupling.cpp"
    ).read_text(errors="replace")
    for token in ("OperationStage::FixedTimeStepBegin",
                  "OperationStage::IterationBegin",
                  "OperationStage::RelaxationApply",
                  "OperationStage::FluxConsistencyRestore",
                  "OperationStage::ConvergenceEvaluate"):
        if token not in (root / "src/solver/system/SF_builtinProviders.cpp").read_text():
            authority_errors.append(
                "selected pressure methods omit lifecycle capability " + token)
    validator_source = (
        root / "src/solver/system/SF_systemValidator.cpp"
    ).read_text(errors="replace")
    if "has no explicit runtime operation ID" not in validator_source:
        authority_errors.append(
            "compiled Plan leaves are not required to carry an OpId")
    # Coupling count authority is now the native source HOW in both backends.
    eulerian_source = (root / "src/solver/system/SF_eulerianCoupling.cpp").read_text(errors="replace")
    for token in ['scope(ExecutionKind::Loop,"EE.outer",request.outerCorrectors)',
                  'scope(ExecutionKind::Loop,"EE.pressure",request.pressureCorrectors)',
                  'scope(ExecutionKind::Loop,"EE.nonOrthogonal",request.nonOrthogonalCorrectors+1)']:
        if token not in eulerian_source:
            authority_errors.append("Eulerian native HOW lost explicit count mapping: " + token)
    if "OpIds::" in eulerian_source or "LegacyExecutionPolicy" in eulerian_source:
        authority_errors.append("Eulerian source HOW contains numerical operations or a duplicate schedule policy")
    for path in ["src/solver/system/SF_executionComposition.cpp", "src/solver/system/SF_pressureCoupling.cpp",
                 "src/solver/system/SF_transformation.cpp"]:
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", (root / path).read_text(errors="replace"), flags=re.S)
        if "sharedPressurePlanFragment" in code or "kSharedPressureScheduleId" in code:
            authority_errors.append("obsolete shared-pressure schedule authority remains: " + path)
    transformer = (root / "src/solver/system/SF_transformation.cpp").read_text(errors="replace")
    shared = transformer[transformer.index("class SharedPressureTransformer"):transformer.index("bool hasEquationId",transformer.index("class SharedPressureTransformer"))]
    if "OpIds::Ee" in shared or "addExecutableOperation" in shared:
        authority_errors.append("shared-pressure transformer declares runtime lifecycle")
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
    for token in ("FormulaMode::", "FormulaCall{"):
        if token in coupling_text:
            authority_errors.append(
                "pressure control skeleton owns Formula execution semantics: "
                + token)
    if "void applyPressureExecution(" not in coupling_text \
            or "bindConstantDensityFormulas" in coupling_text:
        authority_errors.append(
            "constant-density pressure HOW is missing or retains LegacyPlanFragment branching")
    method_source = (
        root / "src/solver/system/SF_methodObjects.cpp"
    ).read_text(errors="replace")
    if "compileLegacyExplicitStep(" not in solve_plan_source:
        authority_errors.append("legacy explicit target inference is not isolated")
    if "selectFormulaBinding(" not in (root / "src/solver/system/SF_termProviderCatalog.cpp").read_text():
        authority_errors.append("production spatial binding bypasses shared Formula precedence")
    fused_equation=(root / "src/solver/equation/compressible/SF_compressible.cpp").read_text()
    if "contains(TermKind" in fused_equation or "makeAssemblyPlan" in fused_equation:
        authority_errors.append("fused production backend still gates terms using legacy definitions")
    if "FormulaMode" in method_source:
        authority_errors.append(
            "migrated EquationMethods retain a FormulaMode authority")
    if "legacyDefinitions.at" in numerical_compiler \
            or "result.terms.emplace_back" in numerical_compiler:
        authority_errors.append(
            "migrated spatial compilation still consumes legacy Equation terms")
    for migrated in (
        "src/solver/algorithm/SF_conservativeRHS.cpp",
        "src/solver/algorithm/pressure/SF_pressureOperators.cpp",
        "src/solver/system/SF_providerResolver.cpp",
    ):
        content=(root / migrated).read_text(errors="replace")
        if "BoundTerm" in content or "numerics.terms" in content:
            authority_errors.append(
                f"migrated numerical consumer retains BoundTerm: {migrated}")
    if "formulaFromEquation(" in (
            root / "src/solver/algorithm/SF_singleFluidStepper.cpp"
            ).read_text(errors="replace"):
        authority_errors.append(
            "single-fluid stage still compares Formula to legacy Equation math")
    for token in ("PressureOperators::", "RhieChow::", "HYPRE_", "MPI_"):
        if token in method_source:
            authority_errors.append(
                "EquationMethod object executes pressure arithmetic/backend: "
                + token)
    builder_source = (
        root / "src/solver/system/SF_systemBuilder.cpp"
    ).read_text(errors="replace")
    if re.search(r'if\s*\(\s*(?:step\.)?subject\s*==\s*"E_PRESSURE"',
                 builder_source):
        authority_errors.append(
            "central composition infers the pressure method from Formula ID")
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
    builtin_provider_text = (root / "src/solver/system/SF_builtinProviders.cpp").read_text()
    for token in ["pressureOperations", "OperationStage::PressureSolve", "OperationStage::VelocityCorrect"]:
        if token not in builtin_provider_text:
            authority_errors.append("pressure provider does not declare compiled operation: " + token)

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
    if "compileOperationBindings" not in resolver_source \
            or "reportOperationBindings" not in resolver_source:
        authority_errors.append("compile-time provider binding/report output is missing")
    if "ProviderCatalog::builtIn()" not in resolver_source:
        authority_errors.append("compiled operation validation bypasses ProviderCatalog")
    explicit_source = (root / "src/solver/algorithm/time/SF_explicit.cpp").read_text()
    for old_table in ("stageFraction[]", "baseWeight[]", "eulerWeight[]"):
        if old_table in explicit_source:
            authority_errors.append(f"explicit provider owns a second tableau: {old_table}")

    # Validate the CMake target graph as well as source includes. An INTERFACE
    # target can reintroduce a link cycle without adding a quoted include.
    target_links: dict[str, set[str]] = defaultdict(set)
    for cmake_file in (root / "src").rglob("CMakeLists.txt"):
        cmake = cmake_file.read_text(errors="replace")
        for target, body in re.findall(
                r"target_link_libraries\s*\(\s*(SF_\w+)\s+([^)]*)\)",
                cmake, flags=re.S):
            target_links[target].update(re.findall(r"\bSF_\w+\b", body))
    visiting: set[str] = set()
    visited: set[str] = set()
    def visit_target(target: str, path: list[str]) -> None:
        if target in visiting:
            authority_errors.append("CMake target cycle: "
                                    + " -> ".join(path + [target]))
            return
        if target in visited:
            return
        visiting.add(target)
        for dependency in sorted(target_links[target]):
            if dependency in target_links:
                visit_target(dependency, path + [target])
        visiting.remove(target)
        visited.add(target)
    for target in sorted(target_links):
        visit_target(target, [])
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

    # WHAT/HOW/WHICH: migrated paths must not reacquire legacy authority.
    native_equation = re.search(r"struct Equation\s*\{(.*?)\n};", formula_header, re.S)
    if not native_equation or re.search(r"\b(?:target|order|unknown)\s*[;=]", native_equation.group(1)):
        authority_errors.append("Equation definition owns HOW storage/order")
    for removed in ("FormulaGroup", "FormulaRegistry", "ProgramStep", "OutputRef", "ProgramControl"):
        if removed in formula_header or removed in solve_program_header:
            authority_errors.append("retired source authority remains: " + removed)
    if "FormulaMode" in solve_program_header:
        authority_errors.append("public HOW exposes legacy numerical realization mode")
    if "std::stable_sort" not in solve_program_header or "a.order < b.order" not in solve_program_header:
        authority_errors.append("execution scopes lack deterministic local ordering")
    if "TargetKind" not in solve_program_header or "node.target=method.target" not in method_source:
        authority_errors.append("lowered equation occurrences lack typed targets")
    density_preset=(root / "src/solver/system/SF_singleFluidPreset.cpp").read_text()
    for retired in ("EquationRole", "Equation::named", "FormulaGroup"):
        if retired in density_preset:
            authority_errors.append("native NS still constructs legacy authority: " + retired)
    pressure_transform=(root / "src/solver/system/SF_transformation.cpp").read_text()
    if "E_MOMENTUM_PREDICTOR" in pressure_transform or 'raw.legacyDefinitions.at("momentum")' in pressure_transform:
        authority_errors.append("pressure coupling duplicates momentum mathematics")
    if "EquationRegistry registry" not in equation_header:
        authority_errors.append("resolved systems lack the authoritative equation registry")
    if "NumericalCompiler::compileSystem" not in builder_source or "builtinProviders()" in builder_source:
        authority_errors.append("SystemBuilder still owns numerical provider compilation")
    source_contribution=(root / "src/models/physics/SF_sourceContribution.cpp").read_text()
    if "extendMathematics" not in source_contribution or "LegacyExecutionPolicy" in source_contribution:
        authority_errors.append("ordinary sources do not extend equation mathematics only")
    pressure_kernel=(root / "src/solver/algorithm/pressure/SF_pressureOperators.cpp").read_text()
    if "workingVelocityView_" not in pressure_kernel or "publishVelocity()" not in pressure_kernel:
        authority_errors.append("working predictor target is not backed by numerical workspace")

    # Strict source/compiled ownership: inspect declarations/function bodies,
    # not incidental comments in files that also contain legacy adapters.
    def cpp_body(content, marker):
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", content, flags=re.S)
        at = code.find(marker)
        if at < 0:
            return ""
        start = code.find("{", at)
        if start < 0:
            return ""
        depth = 0
        for index in range(start, len(code)):
            if code[index] == "{":
                depth += 1
            elif code[index] == "}":
                depth -= 1
                if depth == 0:
                    return code[start + 1:index]
        return ""

    state_registry_header = (root / "src/core/system/SF_stateRegistry.h").read_text()
    state_symbol = cpp_body(state_registry_header, "struct StateSymbol {")
    if not state_symbol or re.search(r"\b(?:EquationRef|NumericalBinding|order|repetitions|terminationSignal)\b", state_symbol):
        authority_errors.append("STATE base symbol owns equation scheduling/numerical binding")
    if "std::vector<UnknownDescriptor> unknowns" in equation_header:
        authority_errors.append("WHAT still owns the retired unknown registration vector")
    if "StateRegistry state" not in equation_header or "CompiledStateView" not in solve_program_header:
        authority_errors.append("Four-module compiler lacks STATE specification/view output")
    method_compiler = (root / "src/solver/system/SF_methodObjects.cpp").read_text()
    if "state.at(step->target.symbol)" not in method_compiler or "realizeTarget(state,compiled.target)" not in method_compiler:
        authority_errors.append("Compiler does not resolve HOW targets through independent STATE input")
    source_target = cpp_body(solve_program_header, "struct Target {")
    source_scope = cpp_body(solve_program_header, "struct ExecutionScope {")
    source_program = cpp_body(solve_program_header, "struct ExecutionProgram {")
    for name, body, forbidden in [
        ("Target", source_target, r"\b(?:workspace|resources|CompiledResourceBinding|storage|offset)\b"),
        ("ExecutionScope", source_scope, r"\b(?:OpId|before|after|NumericalBinding|TemporalMethodBinding)\b"),
        ("ExecutionProgram", source_program, r"\b(?:temporal|TemporalMethodBinding|NumericalBinding|TimeRecipe|consumedPolicies)\b"),
    ]:
        if not body or re.search(forbidden, body):
            authority_errors.append(name + " leaks numerical/runtime authority into source HOW")
    pressure_how = cpp_body(coupling_text, "void applyPressureExecution(")
    pressure_signature = re.search(r"void applyPressureExecution\s*\((.*?)\)\s*\{", coupling_text, re.S)
    if not pressure_signature or "NumericalBinding" in pressure_signature.group(1) \
            or re.search(r"\b(?:NumericalBinding|PressureMomentum|PressureCorrection|OpIds)\b", pressure_how):
        authority_errors.append("pressure HOW mutates WHICH or names runtime operations")
    composition_code = re.sub(r"/\*.*?\*/|//[^\n]*", "",
        (root / "src/solver/system/SF_executionComposition.cpp").read_text(), flags=re.S)
    legacy_numerics_code = re.sub(r"/\*.*?\*/|//[^\n]*", "",
        (root / "src/solver/system/SF_legacyNumerics.cpp").read_text(), flags=re.S)
    if "lowerPressure" in composition_code or "legacyCouplingPlanFragment" in coupling_text:
        authority_errors.append("single-fluid pressure retained legacy execution authority")
    if "conservativePressureScheduleSupported" in legacy_numerics_code:
        authority_errors.append("legacy selector still owns conservative pressure scheduling")
    native_methods = cpp_body(
        (root / "src/solver/system/SF_builtinProviders.cpp").read_text(), "class ConservativePressureMethod")
    if "legacyAdapter" in native_methods or 'return "flow.conservative"' not in native_methods:
        authority_errors.append("conservative pressure methods do not freeze native ownership")
    if "conservativePressureRelations()" not in native_methods or "canonicalFormula" not in native_methods:
        authority_errors.append("conservative pressure relation validation bypasses authoritative WHAT")
    native_transform = cpp_body(pressure_transform, "class PressureConstraintTransformer")
    if re.search(r"\b(?:CompiledResourceBinding|addCompiledEquation|addExecutableOperation|OpIds)\b", native_transform):
        authority_errors.append("native mathematical pressure transformer constructs compiled runtime data")
    executor_body = (root / "src/solver/run/SF_planExecutor.cpp").read_text()
    if re.search(r"\b(?:PISO|SIMPLE|PIMPLE)\b", re.sub(r"/\*.*?\*/|//[^\n]*", "", executor_body, flags=re.S)):
        authority_errors.append("generic executor branches on pressure preset")
    for name, content in [("NumericalCompiler", numerical_compiler), ("execution compiler", method_source)]:
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", content, flags=re.S)
        if re.search(r"\b(?:pressureMultiplier|conservativeTransportedMass|constantDensity|PISO|SIMPLE|PIMPLE|RhieChow)\b", code) \
                or re.search(r'"(?:momentum|continuity|energy|U|rhoU|rhoE|p)"', code):
            authority_errors.append(name + " recreates domain/provider authority")
    if "SF_legacyNumerics.h" in numerical_compiler:
        authority_errors.append("generic spatial compiler infers legacy equation selection")
    composition_source = (root / "src/solver/system/SF_executionComposition.cpp").read_text()
    if "realization." in composition_source:
        authority_errors.append("StateRealization chooses source execution topology")

    # Selected native methods own implementations before the frozen plan is published.
    catalog_header=(root / "src/solver/system/SF_providerCatalog.h").read_text()
    catalog_code=re.sub(r"/\*.*?\*/|//[^\n]*", "", catalog_header, flags=re.S)
    if re.search(r"\b(?:ProviderMatchContext|conservativeState|phaseState|constantPressureSchedule|conservativePressureSchedule)\b", catalog_code):
        authority_errors.append("provider catalog routes from solver-family/state flags")
    binding_body=cpp_body(resolver_source,"compileOperationBindings(")
    if re.search(r"\b(?:conservativeTransportedMass|phaseTransportedState|pressureMultiplier|constantPressureSchedule)\b",binding_body):
        authority_errors.append("generic compiled binding nominates providers from state shape")
    if "leaf.legacyAdapter" not in binding_body or "leaf.provider" not in binding_body:
        authority_errors.append("legacy provider adapter is not isolated from frozen native ownership")
    builder_body=cpp_body(builder_source,"ResolvedSimulationSystem build(")
    if re.search(r"\b(?:compileOperationBindings|resolveOperationBindings|ProviderCatalog|ProviderMatchContext)\b",builder_body):
        authority_errors.append("SystemBuilder owns a second numerical provider selection")
    compiled_leaf=cpp_body(solve_program_header,"struct SolvePlanNode {")
    if not re.search(r"std::string\s+provider\s*;",compiled_leaf):
        authority_errors.append("compiled operation lacks final provider ownership")
    executor_validate=cpp_body(executor_body,"void validateNode(")
    executor_invoke=cpp_body(executor_body,"void executeNode(")
    if "node.provider" not in executor_validate or "node.provider" not in executor_invoke:
        authority_errors.append("PlanExecutor does not validate/invoke frozen operation owners")
    state_add=cpp_body(state_registry_header,"void add(StateSymbol")
    if "find_first_of" not in state_add:
        authority_errors.append("source STATE permits HOW-qualified symbols")
    builtin_state=(root / "src/solver/system/SF_builtinState.cpp").read_text()
    catalog_ctor=cpp_body(builtin_state,"BuiltinStateCatalog::BuiltinStateCatalog()")
    if re.search(r"\b(?:Field|StateBundle|DistributedFieldView|workspaceView)\b|\.add\(",catalog_ctor):
        authority_errors.append("builtin catalog activates/allocates physical state automatically")
    if "installBuiltinState(" in builder_body or "requireTargetStates" not in builder_body:
        authority_errors.append("composition does not activate only requested base STATE")
    stepper=(root / "src/solver/algorithm/SF_singleFluidStepper.cpp").read_text()
    native_runtime=cpp_body(stepper,"SingleFluidStepper::SingleFluidStepper(")+cpp_body(stepper,"SingleFluidStepper::advance(")
    if re.search(r"capabilities\.(?:constantDensity|conservativeState)|\b(?:ProviderMatchContext|resolveOperationBindings)\b",native_runtime):
        authority_errors.append("native runtime chooses an implementation from state flags")

    contribution=(root / "src/models/turbulence/SF_turbulenceSystemContribution.cpp").read_text()
    native_branch=cpp_body(contribution,"if (!spec.eulerian) {")
    if not native_branch or re.search(r"\baddLegacyExecution\s*\(",native_branch):
        authority_errors.append("single-fluid RAS still depends on legacy HOW")
    if not all(re.search(r"\b"+name+r"\s*\(",native_branch) for name in ("addEquation","addState","addExecution","bindNumerics")):
        authority_errors.append("single-fluid RAS must contribute all four native channels")
    # Inspect the registered lambda body, so a comment cannot satisfy the guard.
    begin_body=cpp_body(stepper,"operations.bind(System::OpIds::FlowStepBegin,")
    if not begin_body or re.search(r"\bcorrectTransportModel\s*\(|transportModel\s*->\s*correct",begin_body):
        authority_errors.append("FlowStepBegin still schedules transport correction")
    advance_body=cpp_body(stepper,"SingleFluidStepper::advance(")
    if not re.search(r'bind\(System::OpIds::TurbulenceAdvance,\s*"flow.turbulence"',advance_body):
        authority_errors.append("native turbulence leaf has no frozen runtime owner binding")
    temporal=(root / "src/solver/system/SF_methodObjects.cpp").read_text()
    temporal_body=cpp_body(temporal,"CompiledExecutionProgram compileExecutionProgram(")
    if re.search(r'"(?:k|omega|epsilon|TurbulenceTransport|kOmegaSST|kEpsilon|flow.turbulence)"',temporal_body):
        authority_errors.append("generic mixed temporal compiler branches on turbulence identity")
    if "physicalStepPrefix" not in (root / "src/solver/system/SF_methodObjects.h").read_text():
        authority_errors.append("temporal method lacks an explicit compiled physical-step prefix contract")
    app=(root / "src/app/application/execution/SF_singleFluid.cpp").read_text()
    if re.search(r"\btransportedTurbulence\b",re.sub(r"/\*.*?\*/|//[^\n]*","",app,flags=re.S)):
        authority_errors.append("application independently selects transported turbulence execution")

    phase_pack=cpp_body((root / "src/solver/system/SF_presets.cpp").read_text(),"void addPhaseEquationPack(")
    if "legacyExecution" in phase_pack or not all(token in phase_pack for token in ("addExecution(","bindNumerics(","eulerianPhaseRelations(")):
        authority_errors.append("Eulerian phase pack is not native WHAT/HOW/WHICH")
    eulerian_runtime=(root / "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp").read_text()
    native_bind=cpp_body(eulerian_runtime,"void EulerianStepper::bindSolvePlan(")
    native_advance=cpp_body(eulerian_runtime,"EulerianStepper::advance(")
    if "policyKind" in native_bind or "PlanExecutor::execute" not in native_advance:
        authority_errors.append("Eulerian runtime still selects an independent coupling schedule")
    for name,content in [("execution compiler",method_source),("numerical compiler",numerical_compiler)]:
        code=re.sub(r"/\*.*?\*/|//[^\n]*","",content,flags=re.S)
        if re.search(r'"(?:C_SHARED_PRESSURE|flow\.eulerian-pressure|momentum\.)"|\bEulerianEulerian\b',code):
            authority_errors.append(name+" branches on Eulerian identities")
    methods=(root / "src/solver/system/SF_eulerianMethods.cpp").read_text()
    if not all(token in methods for token in ('"flow.eulerian-pressure"',"fusionMembers=required","StateViewOwner::NumericalProvider","pressureCorrection")):
        authority_errors.append("Eulerian methods lack explicit group/owner/storage contracts")

    # Turbulence has native WHAT/STATE/HOW/WHICH; the flow provider cannot insert its solves.
    for name in ("src/solver/system/SF_eulerianMethods.cpp",
                 "src/solver/algorithm/eulerian/SF_eulerianStepper.cpp",
                 "src/solver/system/SF_legacyNumerics.cpp"):
        content=(root/name).read_text()
        if re.search(r"E_TURB_|legacyTurbulence|legacyClosure|eulerianTurbulenceInputs",content):
            authority_errors.append(name+" retained Eulerian turbulence string/compatibility authority")
    eulerian_contribution=cpp_body(contribution,"void addEulerian(")
    if not eulerian_contribution or "addLegacyExecution" in contribution \
            or not all(token in eulerian_contribution for token in
                ("eulerianTransportMathematics(","eulerianClosureMathematics(","addExecution(","bindNumerics(")):
        authority_errors.append("Eulerian turbulence lacks native contribution channels")
    if "EeTurbulenceSolve" in methods or "EeTurbulencePrepare" in methods:
        authority_errors.append("flow Eulerian provider still inserts turbulence numerical updates")
    native_turbulence=(root/"src/solver/system/SF_eulerianTurbulence.cpp").read_text()
    if "legacyAdapter" in native_turbulence or not all(token in native_turbulence for token in
        ('"flow.eulerian-turbulence"',"fusionMembers=members","canonicalFormula", "validateEulerianTurbulenceBindings")):
        authority_errors.append("Eulerian turbulence lacks frozen native mathematical/group/provider binding")
    for name,content in [("execution compiler",method_source),("numerical compiler",numerical_compiler)]:
        if re.search(r"EulerianTurbulence|eulerianTurbulence|flow\.eulerian-turbulence",content):
            authority_errors.append(name+" interprets domain-specific Eulerian turbulence contracts")
    if not re.search(r'bind\(System::OpIds::EeTurbulenceSolve,\s*"flow.eulerian-turbulence"',eulerian_runtime) \
            or "std::to_string(state.phaseIndex)" not in eulerian_runtime:
        authority_errors.append("Eulerian turbulence runtime lost selected owner or actual phase-slot alias")

    # Native Eulerian assembly is AST/provider-owned; flat DSL remains scoped compatibility.
    native_eulerian=list((root / "src/solver/algorithm/eulerian").rglob("*.h")) \
        +list((root / "src/solver/algorithm/eulerian").rglob("*.cpp")) \
        +list((root / "src/solver/system").glob("SF_eulerian*.h")) \
        +list((root / "src/solver/system").glob("SF_eulerian*.cpp"))
    for path in native_eulerian:
        if path.name.startswith("._"): continue
        code=re.sub(r"/\*.*?\*/|//[^\n]*", "",path.read_text(),flags=re.S)
        if re.search(r"\b(?:AssemblyPlan|AssemblyPlanRegistry|TermKind|legacyDefinitions|legacyEquations|equationDefinition|formulaFromEquation)\b|SF_assemblyPlan\.h",code):
            authority_errors.append("native Eulerian depends on the legacy equation DSL: "+str(path.relative_to(root)))
    for path in (root / "src").rglob("*"):
        if path.suffix not in (".cpp",".h") or path.name.startswith("._"): continue
        code=re.sub(r"/\*.*?\*/|//[^\n]*", "",path.read_text(),flags=re.S)
        if "projectEulerianBackendDefinitions" in code or "eulerianBackendDefinition" in code:
            authority_errors.append("removed Eulerian mathematical projection remains: "+str(path.relative_to(root)))
    phase_pack=cpp_body((root / "src/solver/system/SF_presets.cpp").read_text(),"void addPhaseEquationPack(")
    shared=cpp_body((root / "src/solver/system/SF_transformation.cpp").read_text(),"class SharedPressureTransformer")
    if "EquationDescriptor" in phase_pack or "legacy" in shared or "EquationDescriptor" in shared:
        authority_errors.append("native Eulerian composition creates duplicate legacy mathematics")
    assembly=(root / "src/solver/system/SF_eulerianAssembly.h").read_text()
    if "CompiledEulerianAssemblyContract" not in methods or "canonicalFormula(lhsA)" not in methods \
            or "sourceExtensions(actual.rhs" not in methods or "shared_ptr<const" not in assembly:
        authority_errors.append("Eulerian provider lacks strict AST-derived immutable assembly contracts")
    if "validateEulerianAssemblyBindings" not in eulerian_runtime \
            or any(token in eulerian_runtime for token in ("canonicalFormula(","FormulaExpr","eulerianPhaseRelations(")):
        authority_errors.append("Eulerian runtime interprets mathematical formulas or lacks frozen contract binding")
    if (root / "src/solver/equation/SF_assemblyPlan.h").exists() or (root / "src/solver/equation/SF_assemblyPlan.cpp").exists():
        authority_errors.append("dead AssemblyPlan infrastructure was reintroduced")

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
