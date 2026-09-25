#!/usr/bin/env python3
"""Reorganize native cases into fields/solvers/models/mesh/result folders."""
import json
import shutil
import subprocess
from pathlib import Path

def load(path):
    text = "\n".join(x for x in path.read_text().splitlines()
                       if not x.lstrip().startswith("#"))
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        result = subprocess.run(["ruby", "-rjson", "-ryaml", "-e",
                                 "puts JSON.generate(YAML.load_file(ARGV[0]))", str(path)],
                                check=True, capture_output=True, text=True)
        return json.loads(result.stdout)

def dump(path, value, comment):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(comment + "\n" + json.dumps(value, indent=2, ensure_ascii=False) + "\n")

def reorganize(case):
    fields = case / "fields"
    if (fields / "registry.yaml").exists():
        (fields / "registry.yaml").rename(fields / "fields.yaml")
    solver_registry = case / "solvers/registry.yaml"
    runtime = load(case / "runtime.yaml") if (case / "runtime.yaml").exists() else {}
    numerics = load(case / "numerics.yaml") if (case / "numerics.yaml").exists() else {}
    solver = load(solver_registry).get("solver", {}) if solver_registry.exists() else {}
    dump(case / "solvers/solvers.yaml", {"object": "solvers", "type": "registry",
         "case": {"formatVersion": 1, "name": case.name}, "runtime": runtime,
         "numerics": numerics, "solver": solver}, "# Solver registry and execution inputs.")
    if solver_registry.exists(): solver_registry.unlink()

    models = {"object": "models", "type": "registry"}
    model_dir = case / "model"
    if model_dir.exists():
        for source in sorted(model_dir.glob("*.yaml")):
            category = source.stem
            if category not in {"equations", "constraints", "closures", "geometry"}: continue
            value = load(source)
            models[category] = value.get(category, value)
    geometry_dir = case / "geometry"
    mesh_geometry = case / "mesh/geometry"
    if geometry_dir.exists():
        mesh_geometry.mkdir(parents=True, exist_ok=True)
        for item in geometry_dir.iterdir():
            target = mesh_geometry / item.name
            if not target.exists(): shutil.move(str(item), str(target))
        shutil.rmtree(geometry_dir)
    def rewrite_paths(value):
        if isinstance(value, dict): return {k: rewrite_paths(v) for k, v in value.items()}
        if isinstance(value, list): return [rewrite_paths(v) for v in value]
        if isinstance(value, str) and value.startswith("geometry/"):
            return "mesh/geometry/" + value[len("geometry/"):]
        return value
    models = rewrite_paths(models)
    dump(case / "models/models.yaml", models, "# Physics model registry.")
    for module in (case / "solvers").glob("*.yaml"):
        if module.name == "solvers.yaml": continue
        dump(module, rewrite_paths(load(module)), f"# Solver module: {module.stem}.")
    if model_dir.exists(): shutil.rmtree(model_dir)
    for name in ("runtime.yaml", "numerics.yaml", "solver.yaml"):
        p = case / name
        if p.exists(): p.unlink()
    for p in case.glob("mesh.*"):
        if p.is_file(): p.unlink()

def main():
    cases = sorted(p.parent for p in Path("test").rglob("case.yaml") if "ioRegistry" not in p.parts)
    for case in cases: reorganize(case)
    print(f"reorganized {len(cases)} cases")

if __name__ == "__main__": main()
