#!/usr/bin/env python3
"""Run the frozen WENO7/Rusanov/forward-Euler Sod t=0.2 baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def close(actual: float, expected: float, tolerance: float = 5.0e-12) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def require_close(name: str, actual: float, expected: float) -> None:
    if not close(actual, expected):
        raise RuntimeError(f"{name}: expected {expected:.17g}, got {actual:.17g}")


def read_arrays(path: Path) -> dict[str, tuple[int, list[float]]]:
    root = ET.parse(path).getroot()
    result: dict[str, tuple[int, list[float]]] = {}
    for array in root.findall(".//PointData/DataArray"):
        name = array.attrib.get("Name")
        if name:
            result[name] = (
                int(array.attrib.get("NumberOfComponents", "1")),
                [float(value) for value in (array.text or "").split()],
            )
    return result


def integral(
    arrays: dict[str, tuple[int, list[float]]], name: str, component: int = 0
) -> float:
    # The frozen structured mesh has 101 x 21 x 2 points at unit x/y spacing
    # and dz=0.1. Trapezoidal weights produce a stable physical-domain integral.
    nx, ny, nz = 101, 21, 2
    components, values = arrays[name]
    total = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                weight = (
                    (0.5 if i in (0, nx - 1) else 1.0)
                    * (0.5 if j in (0, ny - 1) else 1.0)
                    * 0.5
                    * 0.1
                )
                index = ((k * ny + j) * nx + i) * components + component
                total += weight * values[index]
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()
    solver = args.solver.resolve()
    case = args.case.resolve()
    baseline = json.loads((case / "baseline.json").read_text())

    with tempfile.TemporaryDirectory(prefix="sonic-sod-rusanov-") as directory:
        working_case = Path(directory) / "case"
        shutil.copytree(case, working_case, ignore=shutil.ignore_patterns("result"))
        run = subprocess.run(
            [str(solver), "run", str(working_case)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if run.returncode != 0:
            print(run.stdout)
            raise RuntimeError(f"sonicSolver exited with status {run.returncode}")
        if "invalid state" in run.stdout.lower():
            raise RuntimeError("EOS invalid state occurred")

        steps = re.findall(
            r"\[SF\] Step time: time=([^,]+), dt=([^\n]+)", run.stdout
        )
        if len(steps) != baseline["steps"]:
            raise RuntimeError(
                f"step count: expected {baseline['steps']}, got {len(steps)}"
            )
        dt_sequence = [float(step[1]) for step in steps]
        for index, (actual, expected) in enumerate(
            zip(dt_sequence, baseline["dt_sequence"])
        ):
            require_close(f"dt[{index}]", actual, expected)
        require_close("final time", float(steps[-1][0]), baseline["final_time"])

        outputs = list((working_case / "result").glob("*_t0.2.vts"))
        if len(outputs) != 1:
            raise RuntimeError("one t=0.2 VTS output was expected")
        output = outputs[0]
        arrays = read_arrays(output)
        rho = arrays["rho"][1]
        pressure = arrays["p"][1]
        if not all(math.isfinite(value) for value in rho + pressure):
            raise RuntimeError("non-finite density or pressure in final state")
        if min(rho) <= 0.0 or min(pressure) <= 0.0:
            raise RuntimeError("non-positive density or pressure in final state")

        velocity_components, velocity = arrays["U"]
        speed = [
            math.sqrt(sum(velocity[i * velocity_components + d] ** 2
                          for d in range(velocity_components)))
            for i in range(len(velocity) // velocity_components)
        ]
        for name, values in (("rho", rho), ("p", pressure), ("speed", speed)):
            require_close(f"{name}.min", min(values), baseline["ranges"][name][0])
            require_close(f"{name}.max", max(values), baseline["ranges"][name][1])

        require_close("mass", integral(arrays, "rho"),
                      baseline["integrals"]["mass"])
        for component, name in enumerate(
            ("MomentumX", "MomentumY", "MomentumZ")
        ):
            require_close(
                f"momentum[{component}]", integral(arrays, name),
                baseline["integrals"]["momentum"][component],
            )
        require_close("energy", integral(arrays, "TotalEnergyDensity"),
                      baseline["integrals"]["energy"])
        digest = hashlib.sha256(output.read_bytes()).hexdigest()
        if digest != baseline["vts_sha256"]:
            raise RuntimeError(f"final VTS SHA-256 changed: {digest}")

    print("WENO7/Rusanov Sod t=0.2 numerical regression passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Sod Rusanov regression failed: {error}")
        raise SystemExit(1)
