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

    with tempfile.TemporaryDirectory(prefix="sonic-cli-recipe-") as directory:
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
        completed = subprocess.run(
            [args.solver, "check", str(runnable)], check=False,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            raise RuntimeError(
                "production case reader rejected CLI-generated TimeRecipe:\n"
                + completed.stdout)
    print("CLI init recipe resolution passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
