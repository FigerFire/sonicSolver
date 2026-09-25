#!/usr/bin/env python3
"""Run the minimal compiled-plan PISO case against its frozen baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


def close(actual: float, expected: float, tolerance: float = 5.0e-12) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def require_close(name: str, actual: float, expected: float) -> None:
    if not close(actual, expected):
        raise RuntimeError(f"{name}: expected {expected:.17g}, got {actual:.17g}")


def match_float(log: str, pattern: str, name: str) -> float:
    match = re.search(pattern, log)
    if not match:
        raise RuntimeError(f"missing {name} diagnostic")
    return float(match.group(1))


def arrays(path: Path) -> dict[str, tuple[float, float, float]]:
    root = ET.parse(path).getroot()
    result: dict[str, tuple[float, float, float]] = {}
    for array in root.findall(".//PointData/DataArray"):
        name = array.attrib.get("Name")
        if not name:
            continue
        values = [float(value) for value in (array.text or "").split()]
        if values:
            result[name] = (
                min(values), max(values),
                math.sqrt(sum(value * value for value in values)))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()

    baseline = json.loads((args.case / "baseline.json").read_text())
    result_dir = args.case / "result"
    shutil.rmtree(result_dir, ignore_errors=True)
    run = subprocess.run(
        [str(args.solver), "run", str(args.case)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False)
    if run.returncode != 0:
        sys.stdout.write(run.stdout)
        raise RuntimeError(f"sonicSolver exited with status {run.returncode}")
    log = run.stdout
    if "RUNTIME STATUS\n  runnable" not in log:
        raise RuntimeError("minimal PISO did not compile to a runnable plan")
    adapter_section = log.split("RUNTIME ADAPTERS", 1)[-1].split(
        "EXECUTION REQUIREMENTS", 1)[0]
    if "LegacyPressureExecutionAdapter" in adapter_section:
        raise RuntimeError("minimal PISO still reports LegacyPressureExecutionAdapter")

    step = baseline["step"]
    require_close("time", match_float(log, r"Step time: time=([^,]+)", "time"), step["time"])
    require_close("dt", match_float(log, r"Step time: time=[^,]+, dt=([^\n]+)", "dt"), step["dt"])
    require_close("pressure residual", match_float(log, r"iterations=\d+, residual=([^,]+)", "pressure residual"), step["pressure_residual"])
    require_close("maxDiv before", match_float(log, r"maxDiv\(before\)=([^,]+)", "maxDiv before"), step["max_divergence_before"])
    require_close("maxDiv after", match_float(log, r"maxDiv\(after\)=([^,]+)", "maxDiv after"), step["max_divergence_after"])
    iterations = int(match_float(log, r"iterations=(\d+)", "iterations"))
    if iterations != step["iterations"]:
        raise RuntimeError(f"iterations: expected {step['iterations']}, got {iterations}")
    hypre = re.search(r"HYPRE\(rebuilds/solves\)=(\d+)/(\d+)", log)
    if not hypre or [int(hypre.group(1)), int(hypre.group(2))] != [
            step["hypre_rebuilds"], step["hypre_solves"]]:
        raise RuntimeError("HYPRE rebuild/solve count differs from legacy baseline")

    outputs = sorted(
        path for path in result_dir.glob("*_000001.vts")
        if not path.name.startswith("._"))
    if len(outputs) != 1:
        raise RuntimeError("one final VTS output was expected")
    observed = arrays(outputs[0])
    for name, expected in baseline["arrays"].items():
        if name not in observed:
            raise RuntimeError(f"missing final array {name}")
        for label, actual in zip(("min", "max", "l2"), observed[name]):
            require_close(f"{name}.{label}", actual, expected[label])

    digest = hashlib.sha256(outputs[0].read_bytes()).hexdigest()
    print("PISO compiled-plan numerical regression passed")
    print(f"final VTS SHA-256: {digest}")
    if digest == baseline["legacy_vts_sha256"]:
        print("final VTS is byte-identical to the frozen legacy result")
    else:
        print("final arrays match the frozen legacy metrics within tolerance")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"PISO regression failed: {error}", file=sys.stderr)
        raise SystemExit(1)
