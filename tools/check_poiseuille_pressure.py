#!/usr/bin/env python3
"""Two-grid analytical channel regression for the serial pressure provider."""

import argparse
import math
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def run_case(solver: Path, source: Path, directory: Path, cells: int):
    case = directory / f"poiseuille-{cells}"
    shutil.copytree(source, case, ignore=shutil.ignore_patterns(
        "result", "mesh.sfm", "mesh.vtm", "mesh.vts", "._*", ".DS_Store"))
    mesh = case / "mesh/blockMeshDict"
    mesh.write_text(mesh.read_text().replace("(16 16 1)",
                                             f"({cells} {cells} 1)"))
    runtime = case / "solvers/runtime.yaml"
    output_step = 2048 * (cells // 16) ** 2
    runtime.write_text(runtime.read_text()
                       .replace("endStep: 2200", f"endStep: {output_step + 1}")
                       .replace("writeInterval: 1000",
                                f"writeInterval: {output_step}"))
    generated = subprocess.run([str(solver), "run", str(case)],
                               capture_output=True, text=True, timeout=120)
    if generated.returncode or not (case / "mesh/mesh.sfm").is_file():
        raise AssertionError("Poiseuille mesh generation failed:\n"
                             + generated.stdout + generated.stderr)
    runtime.write_text(runtime.read_text().replace("createMesh: true",
                                                   "createMesh: false"))
    run = subprocess.run([str(solver), "run", str(case)],
                         capture_output=True, text=True, timeout=180)
    log = run.stdout + run.stderr
    if run.returncode:
        raise AssertionError(f"Poiseuille {cells} failed:\n{log[-5000:]}")
    defects = [float(value) for value in re.findall(
        r"maxDiv\(before/after\)=[\deE+.-]+/([\deE+.-]+)", log)]
    if not defects or defects[-1] > 1.0e-6:
        raise AssertionError(f"Poiseuille {cells} continuity defect: {defects[-1:]}")
    file = case / "result" / f"Poiseuille_{output_step:06d}.vts"
    root = ET.parse(file).getroot()
    def data(name):
        node = root.find(f'.//PointData/DataArray[@Name="{name}"]')
        if node is None:
            raise AssertionError(f"Poiseuille output lacks {name}")
        values = [float(value) for value in node.text.split()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError(f"Poiseuille {name} has non-finite data")
        return values
    velocity, pressure = data("U"), data("p")
    n = cells + 1
    errors = []
    for j in range(1, cells):
        y = -0.5 + j / cells
        actual = velocity[3 * (j * n + cells // 2)]
        exact = (0.25 - y * y) / 2.0  # -dp/dx=0.1, mu=0.1
        errors.append(actual - exact)
    l2 = math.sqrt(sum(error * error for error in errors) / len(errors))
    linf = max(map(abs, errors))
    gradient = pressure[(cells // 2) * n] - pressure[(cells // 2) * n + cells]
    if abs(gradient - 0.1) > 1.0e-7:
        raise AssertionError(f"Poiseuille {cells} pressure drop: {gradient}")
    print(f"Poiseuille {cells}x{cells}: L2={l2:.9g} Linf={linf:.9g} "
          f"maxDiv={defects[-1]:.9g} pressureDrop={gradient:.9g}")
    return l2, linf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="sonic-poiseuille-") as directory:
        errors16 = run_case(args.solver.resolve(), args.case, Path(directory), 16)
        errors32 = run_case(args.solver.resolve(), args.case, Path(directory), 32)
        if not (errors16[0] < 1.2e-3 and errors32[0] < errors16[0] * 0.7
                and errors32[1] < errors16[1] * 0.7):
            raise AssertionError(f"Poiseuille two-grid error did not decrease: "
                                 f"{errors16}, {errors32}")


if __name__ == "__main__":
    main()
