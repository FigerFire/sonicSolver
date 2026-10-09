#!/usr/bin/env python3
"""Read ASCII VTK cylinder outputs without averaging conflicting replicas.

The physical integrals use Cartesian nodal trapezoid weights. Curvilinear
outputs need their metric/volume metadata and are explicitly unsupported here.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

PHYSICAL = ("rho", "p", "U", "MomentumX", "MomentumY", "MomentumZ", "TotalEnergyDensity")


def datasets(result, time):
    entries = {}
    for path in sorted(result.glob("*.pvd")):
        if path.name.startswith("._"): continue
        for item in ET.parse(path).findall(".//DataSet"):
            entries.setdefault(float(item.attrib["timestep"]), set()).add(result / item.attrib["file"])
    if time not in entries:
        raise ValueError(f"Missing requested time {time} in {result}; available={sorted(entries)}")

    def expand(path):
        if path.suffix == ".vts":
            return {path}
        return set().union(*(expand(path.parent / item.attrib["file"])
                             for item in ET.parse(path).findall(".//DataSet")))

    return set().union(*(expand(path) for path in entries[time])), sorted(entries)


def load(result, time):
    paths, times = datasets(result, time)
    points, conflicts, duplicates = {}, {}, 0
    for path in sorted(paths):
        tree = ET.parse(path)
        coordinates = tree.find(".//Points/DataArray")
        if coordinates.attrib.get("format") != "ascii":
            raise ValueError("Only ASCII VTK is supported")
        xyz = list(map(float, coordinates.text.split()))
        n = len(xyz) // 3
        arrays = {}
        for array in tree.findall(".//PointData/DataArray"):
            name = array.attrib["Name"]
            values = list(map(float, array.text.split()))
            width = int(array.attrib.get("NumberOfComponents", 1))
            if len(values) != n * width or not all(map(math.isfinite, values)):
                raise ValueError(f"Invalid storage/nonfinite {name} in {path}")
            arrays[name] = [tuple(values[i * width:(i + 1) * width]) for i in range(n)]
        for i in range(n):
            coordinate = tuple(xyz[3 * i:3 * i + 3])
            row = {name: values[i] for name, values in arrays.items()}
            if coordinate in points:
                duplicates += 1
                if row.keys() != points[coordinate].keys():
                    raise ValueError("Replica output schemas differ")
                for name, values in row.items():
                    if values != points[coordinate][name]:
                        conflicts.setdefault(name, set()).add(coordinate)
            else:
                points[coordinate] = row
    return points, conflicts, {"times": times, "files": len(paths), "uniquePoints": len(points),
                              "duplicates": duplicates,
                              "conflictingReplicas": {key: len(value) for key, value in conflicts.items()}}


def volumes(points):
    axes = [sorted({p[d] for p in points}) for d in range(3)]
    if len(points) != math.prod(map(len, axes)) or any(len(axis) < 2 for axis in axes):
        raise ValueError("Physical volume quadrature requires a complete Cartesian grid of nonzero thickness")
    if any(p not in points for p in itertools.product(*axes)):
        raise ValueError("Curvilinear or incomplete Cartesian output is unsupported")
    weights = []
    for axis in axes:
        weights.append({value: .5 * ((axis[i + 1] - value if i + 1 < len(axis) else 0.)
                                    + (value - axis[i - 1] if i else 0.))
                        for i, value in enumerate(axis)})
    return {p: math.prod(weights[d][p[d]] for d in range(3)) for p in points}


def integrals(points, conflicts, weights):
    result = {"domainVolume": math.fsum(weights.values()),
              "domain": "canonical Eulerian grid, including virtual body fluid"}
    for name in ("rho", "MomentumX", "MomentumY", "MomentumZ", "TotalEnergyDensity",
                 "IBMMultiplierX", "IBMMultiplierY", "IBMMultiplierZ"):
        if name not in next(iter(points.values())):
            continue
        # A conflicting force-density output cannot define a canonical integral.
        result[name] = None if name in conflicts else math.fsum(weights[p] * row[name][0]
                                                                for p, row in points.items())
    for name in ("rho", "p"):
        values = [row[name][0] for row in points.values()]
        result[name + "Min"] = min(values)
        result[name + "Max"] = max(values)
    for name in ("rho", "U", "TotalEnergyDensity"):
        result[name + "WeightedL2"] = math.sqrt(math.fsum(
            weights[p] * math.fsum(v * v for v in row[name]) for p, row in points.items()))
    if "vtkGhostType" in next(iter(points.values())):
        visible = {p for p, row in points.items() if not int(row["vtkGhostType"][0]) & 2}
        result["visibleFluidVolume"] = math.fsum(weights[p] for p in visible)
        result["visibleFluidMass"] = math.fsum(weights[p] * points[p]["rho"][0] for p in visible)
    return result


def compare(current, reference):
    if current.keys() != reference.keys():
        raise ValueError("Physical coordinate sets differ")
    return {name: max(abs(a - b) for p in current for a, b in zip(current[p][name], reference[p][name]))
            for name in PHYSICAL}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--time", type=float, default=5.)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    points, conflicts, report = load(args.result, args.time)
    report["physicalVolumeIntegrals"] = integrals(points, conflicts, volumes(points))
    if args.reference:
        reference, reference_conflicts, ref_report = load(args.reference, args.time)
        if any(name in reference_conflicts for name in PHYSICAL):
            raise ValueError("Reference physical replicas conflict")
        report["referenceReplicaConflicts"] = ref_report["conflictingReplicas"]
        report["physicalMaxAbsDifference"] = compare(points, reference)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if conflicts:
        raise SystemExit("Conflicting replicas are reported and never averaged")


if __name__ == "__main__":
    main()
