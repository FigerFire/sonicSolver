#!/usr/bin/env python3
"""End-to-end serial U/p/rhoConst PISO and collocated-face regression."""

import argparse
import math
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"test"))
from thermophysical_fixture import declare_historical_thermophysical_fixture


def values(root, name):
    array = root.find(f'.//PointData/DataArray[@Name="{name}"]')
    if array is None:
        raise AssertionError(f"missing output field {name}")
    data = [float(value) for value in array.text.split()]
    if not data or not all(math.isfinite(value) for value in data):
        raise AssertionError(f"{name} contains non-finite values")
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="sonic-pressure-") as directory:
        case = Path(directory) / "cavity"
        shutil.copytree(
            args.case, case,
            ignore=shutil.ignore_patterns("result", "._*", ".DS_Store",
                                          "mesh.sfm", "mesh.vtm", "mesh.vts"))
        declare_historical_thermophysical_fixture(case)
        generated = subprocess.run(
            [str(args.solver.resolve()), "run", str(case)],
            capture_output=True, text=True, timeout=120)
        if generated.returncode or not (case / "mesh/mesh.sfm").is_file():
            raise AssertionError("cavity mesh generation failed:\n"
                                 + generated.stdout + generated.stderr)
        runtime = case / "solvers/runtime.yaml"
        runtime.write_text(runtime.read_text().replace("createMesh: true",
                                                       "createMesh: false"))
        explained = subprocess.run(
            [str(args.solver.resolve()), "explain", str(case)],
            capture_output=True, text=True, timeout=120)
        description = explained.stdout + explained.stderr
        for token in ("C_INCOMPRESSIBILITY", "rhoConst=1",
                      "pressure face coupling : RhieChow",
                      "provider : flow.pressure-operators",
                      "provider : flow.rhie-chow", "RUNTIME STATUS\n  runnable"):
            if explained.returncode or token not in description:
                raise AssertionError(f"constant-density explain lacks {token}:\n"
                                     + description[-5000:])
        run = subprocess.run(
            [str(args.solver.resolve()), "run", str(case)],
            capture_output=True, text=True, timeout=120)
        output = run.stdout + run.stderr
        if run.returncode:
            raise AssertionError(f"constant-density PISO exited {run.returncode}:\n"
                                 + output[-5000:])
        defects = [tuple(map(float, pair)) for pair in re.findall(
            r"maxDiv\(before/after\)=([\deE+.-]+)/([\deE+.-]+)", output)]
        if len(defects) != 4:
            raise AssertionError(f"expected two PISO corrections in each of two steps: {defects}")
        for first in (0, 2):
            before, after = defects[first]
            if not before > 0.0 or not after < before * 1.0e-4:
                raise AssertionError(f"first pressure correction did not close continuity: {defects[first]}")
        files = sorted((case / "result").glob("ConstantDensityPiso_*.vts"))
        if len(files) != 3:
            raise AssertionError(f"expected initial + two committed VTK states: {files}")
        root = ET.parse(files[-1]).getroot()
        p = values(root, "p")
        rho = values(root, "rho")
        u = values(root, "U")
        nx = ny = 81
        if len(p) != 2 * nx * ny or len(rho) != len(p) or len(u) != 3 * len(p):
            raise AssertionError("unexpected structured field shape")
        if any(abs(value - 1.0) > 1.0e-12 for value in rho):
            raise AssertionError("derived rhoConst changed during PISO")
        if abs(p[nx + 1] - 101325.0) > 1.0e-8:
            raise AssertionError("pressure reference cell drifted")
        if abs(u[3 * (80 * nx + 40)] - 1.0) > 1.0e-12:
            raise AssertionError("moving lid fixedValue U was not applied")
        if abs(u[3 * 40]) > 1.0e-12:
            raise AssertionError("bottom wall no-slip U was not applied")
        interior = [p[j * nx + i] - 101325.0
                    for j in range(1, ny - 1) for i in range(1, nx - 1)]
        rms = math.sqrt(sum(value * value for value in interior) / len(interior))
        parity = abs(sum((-1 if (i+j) % 2 else 1)
                         * (p[j * nx + i] - 101325.0)
                         for j in range(1, ny - 1) for i in range(1, nx - 1))
                     / len(interior))
        if rms <= 0.0 or parity / rms >= 0.05:
            raise AssertionError(f"collocated pressure checkerboard signal: {parity/rms}")
        print(f"constant-density PISO: first corrections {defects[0]}, {defects[2]}; "
              f"checkerboard/rms={parity/rms:.3e}; reference=101325")


if __name__ == "__main__":
    main()
