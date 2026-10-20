#!/usr/bin/env python3
"""Verify that CLI-generated numerics can be resolved by the production reader."""

import argparse
import pathlib
import shutil
import subprocess
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", required=True)
    parser.add_argument("--case", required=True)
    args = parser.parse_args()

    temp_root = pathlib.Path(__file__).resolve().parents[1] / "test/t"
    temp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sonic-cli-recipe-", dir=temp_root) as directory:
        root = pathlib.Path(directory)
        generated = root / "generated"
        subprocess.run(
            [args.solver, "init", str(generated), "--recipe",
             "compressibleSingleFluid"], check=True)
        generated_numerics = generated / "solvers" / "numerics.yaml"
        text = generated_numerics.read_text()
        if "default: forwardEuler" not in text or "convection: teno5Steger" not in text:
            raise RuntimeError("CLI template did not emit canonical recipe names")

        runnable = root / "runnable"
        shutil.copytree(args.case, runnable)
        shutil.copy2(generated_numerics, runnable / "solvers" / "numerics.yaml")
        for file in ("state.yaml", "equations.yaml", "equationRegistry.yaml", "execution.yaml"):
            path = generated / "solvers" / file
            if "object: solver" not in path.read_text():
                raise RuntimeError("CLI registry object must match solvers directory")
            if "  type: " + path.stem + "\n" not in path.read_text():
                raise RuntimeError("CLI file identity must match its registry key")
            shutil.copy2(path, runnable / "solvers" / file)
        for old in ("state", "equations", "algorithms"):
            if (generated / old).exists():
                raise RuntimeError("CLI created a separate solver registry directory")
        shutil.copy2(generated / "models/thermoDynamics.yaml", runnable / "models/thermoDynamics.yaml")
        registry=runnable / "models/models.yaml"
        if "type: thermoDynamics" not in registry.read_text():
            registry.write_text(registry.read_text()+"thermoDynamics:\n  type: thermoDynamics\n  file: models/thermoDynamics.yaml\n")
        completed = subprocess.run(
            [args.solver, "check", str(runnable)], check=False,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            raise RuntimeError(
                "production case reader rejected CLI-generated TimeRecipe:\n"
                + completed.stdout)
        algorithm = runnable / "solvers" / "algorithm.yaml"
        algorithm.write_text("SonicFile:\n  object: solver\n  type: algorithm\n"
                             "type: densityBase\nalgorithm: PIMPLE\n"
                             "PIMPLE:\n  outerCorrectors: 1\n")
        (runnable / "solvers/execution.yaml").write_text(
            "SonicFile:\n  object: solver\n  type: execution\nPIMPLE: {}\n")
        # The native algorithm has complete numerical controls, independent of old type labels.
        pressure_controls=pathlib.Path(args.case).parents[1] / "pressure/constantDensityPiso/solvers/algorithm.yaml"
        if pressure_controls.exists():
            shutil.copy2(pressure_controls, algorithm)
        invalid = subprocess.run([args.solver, "check", str(runnable)],
                                 text=True, capture_output=True)
        if invalid.returncode == 0 or "constraint" not in (invalid.stdout + invalid.stderr).lower():
            raise RuntimeError("Explicit inapplicable coupling was silently ignored")
    print("CLI init recipe resolution passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
