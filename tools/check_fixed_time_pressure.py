#!/usr/bin/env python3
"""Run the shared constant-density pressure operations under three plans."""

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def invoke(solver, action, case, trace=False):
    env = os.environ.copy()
    if trace:
        env["SF_PLAN_TRACE"] = "1"
    result = subprocess.run([str(solver), action, str(case)], env=env,
                            capture_output=True, text=True, timeout=120)
    output = result.stdout + result.stderr
    if result.returncode:
        raise AssertionError(f"{action} {case.name} failed:\n{output[-5000:]}")
    return output


def prepare(solver, source, parent, name, preset, outer, corrections,
            relaxation=1.0, steps=None, relative_tolerance=None,
            absolute_tolerance=None):
    case = parent / name
    shutil.copytree(source, case, ignore=shutil.ignore_patterns(
        "result", "mesh.sfm", "mesh.vtm", "mesh.vts", "._*", ".DS_Store"))
    algorithm = case / "solvers/algorithm.yaml"
    contents = (algorithm.read_text()
                .replace("algorithm: PISO", f"algorithm: {preset}")
                .replace("PISO:", f"{preset}:")
                .replace("outerCorrectors: 1", f"outerCorrectors: {outer}")
                .replace("pressureCorrectors: 2",
                         f"pressureCorrectors: {corrections}")
                .replace("momentumRelaxation: 1.0",
                         f"momentumRelaxation: {relaxation}")
                .replace("pressureRelaxation: 1.0",
                         f"pressureRelaxation: {relaxation}"))
    if relative_tolerance is not None or absolute_tolerance is not None:
        lines = []
        if relative_tolerance is not None:
            lines.append(f"  relativeTolerance: {relative_tolerance}")
        if absolute_tolerance is not None:
            lines.append(f"  absoluteTolerance: {absolute_tolerance}")
        anchor = f"  pressureRelaxation: {relaxation}"
        if anchor not in contents:
            raise AssertionError("pressure convergence controls have no preset block")
        contents = contents.replace(anchor, anchor + "\n" + "\n".join(lines), 1)
    algorithm.write_text(contents)
    semantic = case / "algorithms/algorithms.yaml"
    if semantic.exists():
        semantic.write_text("SonicFile:\n  object: algorithms\n  type: registry\n"
                            + f"{preset}: {{}}\n")
    runtime = case / "solvers/runtime.yaml"
    if steps is not None:
        runtime.write_text(runtime.read_text()
                           .replace("endStep: 2200", f"endStep: {steps + 1}")
                           .replace("writeInterval: 1000",
                                    f"writeInterval: {steps}"))
    invoke(solver, "run", case)
    if not (case / "mesh/mesh.sfm").is_file():
        raise AssertionError(f"mesh generation did not complete for {name}")
    runtime.write_text(runtime.read_text().replace("createMesh: true",
                                                   "createMesh: false"))
    explain = invoke(solver, "explain", case)
    for required in ("C_INCOMPRESSIBILITY", "RUNTIME STATUS\n  runnable",
                     "provider : flow.pressure-operators",
                     "provider : flow.rhie-chow"):
        if required not in explain:
            raise AssertionError(f"{name} explain lacks {required}")
    return case


def field(case, job, step, name):
    root = ET.parse(case / "result" / f"{job}_{step:06d}.vts").getroot()
    node = root.find(f'.//PointData/DataArray[@Name="{name}"]')
    if node is None:
        raise AssertionError(f"{case.name} lacks {name}")
    values = [float(value) for value in node.text.split()]
    if not values or not all(math.isfinite(value) for value in values):
        raise AssertionError(f"{case.name} has non-finite {name}")
    return values


def clock(log):
    return [tuple(map(float, match)) for match in re.findall(
        r"Step time: time=([\deE+.-]+), dt=([\deE+.-]+)", log)]


def check_cavity(solver, source, directory):
    cases = {
        "piso": prepare(solver, source, directory, "piso", "PISO", 1, 2),
        "pimple1": prepare(solver, source, directory, "pimple1", "PIMPLE", 1, 2),
        "simple2": prepare(solver, source, directory, "simple2", "SIMPLE", 2, 1,
                           0.8),
        "pimple2": prepare(solver, source, directory, "pimple2", "PIMPLE", 2, 2),
        "simpleEarly": prepare(solver, source, directory, "simpleEarly",
                               "SIMPLE", 4, 1, 0.8,
                               relative_tolerance=0.004,
                               absolute_tolerance=5e-8),
        "pimpleEarly": prepare(solver, source, directory, "pimpleEarly",
                               "PIMPLE", 4, 2, 0.8,
                               relative_tolerance=0.004,
                               absolute_tolerance=5e-8),
    }
    logs = {name: invoke(solver, "run", case, trace=True)
            for name, case in cases.items()}
    reference = clock(logs["piso"])
    if len(reference) != 2:
        raise AssertionError("PISO did not commit two physical steps")
    for name, log in logs.items():
        if clock(log) != reference:
            raise AssertionError(f"{name} advanced a different physical clock")
        expected_outer = 2 if name in ("simple2", "pimple2",
                                       "simpleEarly", "pimpleEarly") else 1
        if name != "piso":
            if log.count("op=pressure.step.begin") != 2:
                raise AssertionError(f"{name} did not freeze one base per step")
            if log.count("op=pressure.iteration.begin") != 2 * expected_outer:
                raise AssertionError(f"{name} did not execute its outer loop")
            if log.count("op=pressure.step.commit") != 2:
                raise AssertionError(f"{name} committed within its outer loop")
            corrections = 1 if name.startswith("simple") else 2
            if log.count("op=pressure.solve") != 2*expected_outer*corrections:
                raise AssertionError(f"{name} changed its inner correction count")
            operations = re.findall(r"op=(pressure\.[a-z.]+)", log)
            iteration_ends = [i for i, op in enumerate(operations)
                              if op == "pressure.iteration.end"]
            if len(iteration_ends) != 2*expected_outer:
                raise AssertionError(f"{name} changed its iteration end count")
            for end in iteration_ends:
                if operations[end-3:end] != [
                        "pressure.relaxation.apply",
                        "pressure.flux.consistency.restore",
                        "pressure.convergence.evaluate"]:
                    raise AssertionError(
                        f"{name} lost relaxation/flux/convergence order: "
                        f"{operations[end-3:end]}")
            final = [float(value) for value in re.findall(
                r"maxDelta=[^\n]*?maxDiv=([\deE+.-]+)", log)]
            candidate = [float(value) for value in re.findall(
                r"maxDelta=[^\n]*?candidateDiv=([\deE+.-]+)", log)]
            flux_delta = [float(value) for value in re.findall(
                r"maxDelta=[^\n]*?fluxDelta=([\deE+.-]+)", log)]
            converged = re.findall(
                r"maxDelta=[^\n]*?converged=(true|false)", log)
            continuity = [float(value) for value in re.findall(
                r"State closure: iterations=[^\n]*?maxDiv\(before/after\)="
                r"[\deE+.-]+/([\deE+.-]+)", log)]
            if len(final) != 2*expected_outer or len(candidate) != len(final) \
                    or len(flux_delta) != len(final) \
                    or any(not math.isfinite(value) for value in flux_delta) \
                    or len(continuity) != len(final) \
                    or len(converged) != len(final) \
                    or any(abs(a-b) > 1e-12 for a,b in zip(final,continuity)):
                raise AssertionError(f"{name} reports candidate, not final flux continuity")
            if name in ("simpleEarly", "pimpleEarly") and any(
                    status != "true" for status in converged[1::2]):
                raise AssertionError(f"{name} exited without reporting convergence")
            if name.startswith("simple") or name == "pimpleEarly":
                if not any(abs(a-b) > 1e-14 for a,b in zip(final,candidate)):
                    raise AssertionError(f"{name} retained the unrelaxed candidate flux")
        defects = [tuple(map(float, pair)) for pair in re.findall(
            r"maxDiv\(before/after\)=([\deE+.-]+)/([\deE+.-]+)", log)]
        if not defects or any(not math.isfinite(after) or after > 1e-6
                              for _, after in defects):
            raise AssertionError(f"{name} continuity did not remain finite: {defects[-4:]}")
    for name in ("U", "p", "rho"):
        piso = field(cases["piso"], "ConstantDensityPiso", 2, name)
        pimple = field(cases["pimple1"], "ConstantDensityPiso", 2, name)
        delta = max(abs(a - b) for a, b in zip(piso, pimple))
        if len(piso) != len(pimple) or delta > 1e-12:
            raise AssertionError(f"PIMPLE outer=1 differs from PISO in {name}: {delta}")
    if abs(field(cases["simple2"], "ConstantDensityPiso", 2, "p")[82]
           - 101325.0) > 1e-8:
        raise AssertionError("SIMPLE pressure relaxation changed its gauge")
    zero = prepare(solver, source, directory, "zero-gauge", "PISO", 1, 2)
    algorithm = zero / "solvers/algorithm.yaml"
    algorithm.write_text(algorithm.read_text().replace(
        "referencePressure: 101325.0", "referencePressure: 0.0"))
    invoke(solver, "run", zero)
    if abs(field(zero, "ConstantDensityPiso", 2, "p")[82]) > 1e-12:
        raise AssertionError("zero multiplier-pressure gauge was not preserved")
    print(f"fixed-time cavity: same clock {reference}; "
          "PIMPLE outer=1 matches PISO; SIMPLE/PIMPLE max=4 exit at 2 per step; "
          "relaxed flux continuity checked; zero gauge legal")


def check_simple_poiseuille(solver, source, directory):
    case = prepare(solver, source, directory, "simple-poiseuille",
                   "SIMPLE", 2, 1, 0.8, steps=512)
    log = invoke(solver, "run", case)
    deltas = [float(value) for value in re.findall(
        r"maxDelta=([\deE+.-]+), velocityDelta=[\deE+.-]+, "
        r"pressureDelta=[\deE+.-]+, maxDiv=", log)]
    convergence = re.findall(
        r"maxDelta=[^\n]*?converged=(true|false)", log)
    if len(deltas) != 2 * 513 or any(deltas[i+1] >= deltas[i]
                                       for i in range(0, len(deltas), 2)):
        raise AssertionError("SIMPLE outer fixed-point delta did not decrease")
    if len(convergence) != len(deltas):
        raise AssertionError("SIMPLE did not report every convergence decision")
    defects = [float(value) for value in re.findall(
        r"maxDiv\(before/after\)=[\deE+.-]+/([\deE+.-]+)", log)]
    if not defects or defects[-1] > 1e-6:
        raise AssertionError("SIMPLE continuity residual did not converge")
    velocity = field(case, "Poiseuille", 512, "U")
    pressure = field(case, "Poiseuille", 512, "p")
    n = 17
    errors = []
    for j in range(1, 16):
        y = -0.5 + j / 16
        errors.append(velocity[3 * (j * n + 8)] - 0.5 * (0.25 - y*y))
    l2 = math.sqrt(sum(value*value for value in errors) / len(errors))
    linf = max(abs(value) for value in errors)
    pressure_drop = pressure[8*n] - pressure[8*n+16]
    if l2 >= 0.02 or abs(pressure_drop-0.1) > 1e-7:
        raise AssertionError(f"SIMPLE Poiseuille profile/drop: {l2}, {pressure_drop}")
    print(f"SIMPLE Poiseuille: L2={l2:.9g}, Linf={linf:.9g}, "
          f"maxDiv={defects[-1]:.9g}, pressureDrop={pressure_drop:.9g}, "
          f"outerIterations=2/step, finalMaxDelta={deltas[-1]:.9g}, "
          f"finalConverged={convergence[-1]}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--cavity", type=Path, required=True)
    parser.add_argument("--poiseuille", type=Path, required=True)
    args = parser.parse_args()
    solver = args.solver.resolve()
    with tempfile.TemporaryDirectory(prefix="sonic-fixed-time-") as temp:
        directory = Path(temp)
        check_cavity(solver, args.cavity, directory)
        check_simple_poiseuille(solver, args.poiseuille, directory)


if __name__ == "__main__":
    main()
