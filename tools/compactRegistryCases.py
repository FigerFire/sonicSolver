#!/usr/bin/env python3
"""Convert registry cases to the compact fields/boundaries syntax."""
import json
import re
import shutil
from pathlib import Path
import sys

def write(path, value, comment):
    path.write_text(comment + "\n" + dump_yaml(value) + "\n")

def scalar(value):
    if value is None: return "null"
    if value is True: return "true"
    if value is False: return "false"
    if isinstance(value, str):
        # Keep condition names readable while quoting values that YAML could
        # reinterpret as a boolean/null/number.
        if (re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", value)
                and value not in {"true", "false", "null", "yes", "no", "on", "off"}):
            return value
        return json.dumps(value, ensure_ascii=False)
    return str(value)

def dump_yaml(value, level=0):
    indent = "  " * level
    if isinstance(value, dict):
        lines = []
        for key, child in value.items():
            name = str(key)
            if (isinstance(child, dict) and len(child) == 1
                    and isinstance(next(iter(child.values())), (str, int, float, bool, list))):
                child_key, child_value = next(iter(child.items()))
                lines.append(f"{indent}{name}: {{{child_key}: {scalar(child_value)}}}")
                continue
            if isinstance(child, (dict, list)) and child:
                if isinstance(child, list) and all(not isinstance(x, (dict, list)) for x in child):
                    lines.append(f"{indent}{name}: [{', '.join(scalar(x) for x in child)}]")
                else:
                    lines.append(f"{indent}{name}:")
                    lines.append(dump_yaml(child, level + 1))
            else:
                lines.append(f"{indent}{name}: {scalar(child)}")
        return "\n".join(lines)
    if isinstance(value, list):
        return "\n".join(f"{indent}- {scalar(x)}" for x in value)
    return indent + scalar(value)

def load(path):
    return json.loads("\n".join(line for line in path.read_text().splitlines()
                                  if not line.lstrip().startswith("#")))

def compact_boundary(value):
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise ValueError(f"invalid boundary value: {value!r}")
    if "type" in value:
        typ = value["type"]
        if typ in {"zeroGradient", "empty", "symmetry", "slip", "noSlip"}:
            return typ
        if "value" in value:
            return {typ: value["value"]}
        return typ
    if len(value) == 1:
        key, child = next(iter(value.items()))
        return {key: child}
    raise ValueError(f"boundary must contain one type: {value!r}")

def compact_case(case, legacy_root=None):
    if legacy_root:
        from migrateRegistryCases import parse
    fields_path = case / "model/fields.yaml"
    boundaries_path = case / "model/boundaries.yaml"
    old_fields = load(fields_path) if not legacy_root else {}
    old_boundaries = load(boundaries_path) if not legacy_root else {}
    old_fields = old_fields.get("fields", old_fields)
    old_boundaries = old_boundaries.get("boundaries", old_boundaries)
    fields = {}
    for name, value in old_fields.items():
        if name == "default":
            fields[name] = value
            continue
        if isinstance(value, dict):
            if "default" in value:
                fields[name] = value["default"]
            elif "value" in value:
                fields[name] = value["value"]
            elif "initial" in value:
                fields[name] = value["initial"]
            else:
                fields[name] = None
        else:
            fields[name] = value
    registry = {}
    internal = {}
    if legacy_root:
        legacy_zero = legacy_root / case.relative_to(Path("test")) / "0"
        if legacy_zero.is_dir():
            from migrateRegistryCases import parse
            for source in sorted(legacy_zero.iterdir()):
                if source.name.startswith(".") or not source.is_file():
                    continue
                data = parse(source)
                initial = data.get("internalField")
                sets = data.get("internalSets", {})
                boundaries = data.get("boundaryField", {})
                old_boundaries[source.name] = boundaries
                typ = "vector" if isinstance(initial, list) and len(initial) == 3 else "scalar"
                dimensions = data.get("dimensions", [0, 0, 0, 0, 0, 0, 0])
                if not isinstance(dimensions, list) or len(dimensions) != 7:
                    dimensions = [0, 0, 0, 0, 0, 0, 0]
                registry[source.name] = {"type": typ, "domain": "fluid", "location": "point",
                                         "dimensions": dimensions, "storage": "primary"}
                if sets:
                    compact = {"default": initial}
                    for set_name, condition in sets.items():
                        compact[set_name] = condition.get("value", condition)
                    internal[source.name] = compact
                else:
                    internal[source.name] = {"uniform": initial}
        solver_source = legacy_root / case.relative_to(Path("test")) / "system/solverProperties"
        solver = parse(solver_source) if solver_source.exists() else {"type": "densityBase"}
        phase_fields = set(registry)
        for name, descriptor in registry.items():
            base = name.split('.', 1)[0]
            descriptor["output"] = (name in {"rho", "U", "p"}
                                     or (name == "T" and not any(x.startswith("T.") for x in phase_fields))
                                     or ("." in name and base in {"T", "U", "rho", "alpha"}))
    else:
        for name, value in fields.items():
            if name == "default": continue
            typ = "vector" if isinstance(value, list) and len(value) == 3 else "scalar"
            registry[name] = {"type": typ, "domain": "fluid", "location": "point",
                              "dimensions": [0, 0, 0, 0, 0, 0, 0], "storage": "primary"}
            internal[name] = {"uniform": value}
        solver = {"type": "densityBase"}
    boundaries = {}
    for name, patch_map in old_boundaries.items():
        if name == "default":
            boundaries[name] = compact_boundary(patch_map)
            continue
        boundaries[name] = {patch: compact_boundary(condition)
                            for patch, condition in patch_map.items()}
    fields_dir = case / "fields"
    solvers_dir = case / "solvers"
    fields_dir.mkdir(exist_ok=True)
    solvers_dir.mkdir(exist_ok=True)
    write(fields_dir / "registry.yaml", {"object": "fields", "type": "registry", **registry},
          "# Field registry: names, types and storage semantics.")
    write(fields_dir / "internalField.yaml", {"object": "fields", "type": "internalField", **internal},
          "# Internal field values; uniform is the default initial condition.")
    write(fields_dir / "boundaries.yaml", {"object": "fields", "type": "boundary", **boundaries},
          "# Boundary conditions: patch: {fixedValue: value} or patch: zeroGradient.")
    write(solvers_dir / "registry.yaml", {"object": "solver", "type": "registry", "solver": solver},
          "# Solver module registry.")
    constraints_path = case / "model/constraints.yaml"
    if constraints_path.exists():
        constraints_doc = load(constraints_path)
        constraints_map = constraints_doc.get("constraints", constraints_doc)
        retained = {}
        for name, item in constraints_map.items():
            if item.get("type") in {"ImmersedBoundary", "ILW"}:
                module_name = "IBM" if item["type"] == "ImmersedBoundary" else "ILW"
                write(solvers_dir / f"{module_name}.yaml",
                      {"object": "solver", "type": module_name, **item.get("parameters", {})},
                      f"# Solver module: {module_name}.")
            else:
                retained[name] = item
        write(constraints_path, {"constraints": retained}, "# Generic model constraints.")
    for old in (fields_path, boundaries_path):
        if old.exists(): old.unlink()
    old_solver = case / "solver.yaml"
    if old_solver.exists(): old_solver.unlink()

def main():
    root = Path("test")
    legacy_root = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    cases = sorted(p.parent for p in root.rglob("case.yaml")
                   if "ioRegistry" not in p.parts)
    for case in cases:
        compact_case(case, legacy_root)
    print(f"compacted {len(cases)} cases")

if __name__ == "__main__":
    main()
