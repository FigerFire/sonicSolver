#!/usr/bin/env python3
"""逐步比较两次 SonicSolver VTK 输出，不调用或修改求解器。

用法：
  python3 tools/compare_high_order_steps.py OLD_RESULT CURRENT_RESULT
  python3 tools/compare_high_order_steps.py OLD_RESULT CURRENT_RESULT \\
      --abs-tol 1e-12 --rel-tol 1e-10

默认只报告差异，不声明 FIRST_DIVERGENCE_STEP。必须由调查者显式提供绝对和
相对阈值，才会基于两者判定 first divergence；输出同时保留原始差值和相对
machine-roundoff scale，避免把任意阈值误写成数值结论。
"""

from __future__ import annotations

import argparse
import math
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

EPSILON = sys.float_info.epsilon
TARGETS = {
    "rho": ("rho",),
    "rhoU": ("rhoU", "MomentumX"),
    "rhoV": ("rhoV", "MomentumY"),
    "rhoW": ("rhoW", "MomentumZ"),
    "rhoE": ("rhoE", "TotalEnergyDensity"),
    "p": ("p", "Pressure"),
}


@dataclass(frozen=True)
class Frame:
    ordinal: int
    time: float | None
    pieces: tuple[Path, ...]


@dataclass
class Difference:
    count: int = 0
    l1: float = 0.0
    l2_squared: float = 0.0
    linf: float = 0.0
    max_rel: float = 0.0
    scale: float = 1.0

    def add(self, old: float, current: float) -> None:
        difference = abs(current - old)
        magnitude = max(abs(old), abs(current))
        self.count += 1
        self.l1 += difference
        self.l2_squared += difference * difference
        self.linf = max(self.linf, difference)
        self.scale = max(self.scale, magnitude)
        if difference == 0.0:
            return
        self.max_rel = math.inf if magnitude == 0.0 else max(
            self.max_rel, difference / magnitude)

    @property
    def l2(self) -> float:
        return math.sqrt(self.l2_squared)

    @property
    def roundoff_scale(self) -> float:
        return EPSILON * self.scale


def data_arrays(path: Path) -> dict[str, list[float]]:
    root = ET.parse(path).getroot()
    arrays: dict[str, list[float]] = {}
    for array in root.findall(".//PointData/DataArray"):
        name = array.get("Name")
        if not name:
            continue
        if array.get("format", "ascii") != "ascii":
            raise ValueError(f"{path}: DataArray '{name}' is not ASCII VTK data")
        arrays[name] = [float(value) for value in (array.text or "").split()]
    return arrays


def vtm_pieces(vtm: Path) -> tuple[Path, ...]:
    root = ET.parse(vtm).getroot()
    files = [item.get("file") for item in root.findall(".//DataSet")]
    pieces = tuple(vtm.parent / filename for filename in files if filename)
    if not pieces:
        raise ValueError(f"{vtm}: VTM does not reference any VTS piece")
    missing = [str(piece) for piece in pieces if not piece.is_file()]
    if missing:
        raise ValueError(f"{vtm}: missing VTS pieces: {', '.join(missing)}")
    return pieces


def series_from_pvd(directory: Path, pvd: Path) -> list[Frame]:
    root = ET.parse(pvd).getroot()
    frames: list[Frame] = []
    for ordinal, dataset in enumerate(root.findall(".//DataSet")):
        filename = dataset.get("file")
        if not filename:
            raise ValueError(f"{pvd}: DataSet {ordinal} has no file")
        output = directory / filename
        if not output.is_file():
            raise ValueError(f"{pvd}: output '{filename}' does not exist")
        time_text = dataset.get("timestep")
        time = float(time_text) if time_text is not None else None
        pieces = vtm_pieces(output) if output.suffix == ".vtm" else (output,)
        frames.append(Frame(ordinal, time, pieces))
    if not frames:
        raise ValueError(f"{pvd}: no DataSet entries")
    return frames


def load_series(directory: Path) -> list[Frame]:
    pvd_files = sorted(directory.glob("*.pvd"))
    if len(pvd_files) == 1:
        return series_from_pvd(directory, pvd_files[0])
    if len(pvd_files) > 1:
        raise ValueError(
            f"{directory}: multiple PVD files; supply a directory containing one output series")
    pieces = tuple(sorted(directory.glob("*.vts")))
    if not pieces:
        raise ValueError(f"{directory}: no PVD or VTS output found")
    return [Frame(index, None, (piece,)) for index, piece in enumerate(pieces)]


def target_array(arrays: dict[str, list[float]], semantic: str) -> tuple[str, list[float]]:
    for candidate in TARGETS[semantic]:
        if candidate in arrays:
            return candidate, arrays[candidate]
    choices = ", ".join(TARGETS[semantic])
    raise ValueError(f"missing required array for {semantic}; expected one of {choices}")


def compare_frame(old: Frame, current: Frame) -> dict[str, Difference]:
    if len(old.pieces) != len(current.pieces):
        raise ValueError(
            f"output {old.ordinal}: old has {len(old.pieces)} pieces, current has {len(current.pieces)}")
    differences = {name: Difference() for name in TARGETS}
    for piece_index, (old_piece, current_piece) in enumerate(zip(old.pieces, current.pieces)):
        old_arrays = data_arrays(old_piece)
        current_arrays = data_arrays(current_piece)
        for semantic in TARGETS:
            old_name, old_values = target_array(old_arrays, semantic)
            current_name, current_values = target_array(current_arrays, semantic)
            if len(old_values) != len(current_values):
                raise ValueError(
                    f"output {old.ordinal}, piece {piece_index}, {semantic}: "
                    f"old array {old_name} has {len(old_values)} values, current array "
                    f"{current_name} has {len(current_values)}")
            for old_value, current_value in zip(old_values, current_values):
                differences[semantic].add(old_value, current_value)
    return differences


def differs_beyond_tolerance(
    values: Iterable[Difference], abs_tol: float | None, rel_tol: float | None
) -> bool | None:
    if abs_tol is None or rel_tol is None:
        return None
    return any(
        item.linf > abs_tol + rel_tol * item.scale
        for item in values
    )


def format_time(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.12e}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old_result", type=Path)
    parser.add_argument("current_result", type=Path)
    parser.add_argument("--abs-tol", type=float, default=None,
                        help="first-divergence absolute tolerance; requires --rel-tol")
    parser.add_argument("--rel-tol", type=float, default=None,
                        help="first-divergence relative tolerance; requires --abs-tol")
    args = parser.parse_args()
    if (args.abs_tol is None) != (args.rel_tol is None):
        parser.error("--abs-tol and --rel-tol must be supplied together")
    if args.abs_tol is not None and (args.abs_tol < 0.0 or args.rel_tol < 0.0):
        parser.error("comparison tolerances must be non-negative")

    try:
        old_frames = load_series(args.old_result)
        current_frames = load_series(args.current_result)
        compared_frames = min(len(old_frames), len(current_frames))
        if compared_frames == 0:
            raise ValueError("output series have no shared frame ordinal")
        print("step old_time current_time time_abs rho_Linf rhoE_Linf p_Linf max_rel roundoff_scale classification")
        first_difference: int | None = None
        first_divergence: int | None = None
        first_time_difference: int | None = None
        for old, current in zip(old_frames[:compared_frames], current_frames[:compared_frames]):
            time_difference = (
                abs(old.time - current.time)
                if old.time is not None and current.time is not None else math.nan)
            if not math.isnan(time_difference) and time_difference != 0.0 and first_time_difference is None:
                first_time_difference = old.ordinal
            values = compare_frame(old, current)
            any_difference = any(item.linf != 0.0 for item in values.values())
            if any_difference and first_difference is None:
                first_difference = old.ordinal
            exceeds = differs_beyond_tolerance(values.values(), args.abs_tol, args.rel_tol)
            if exceeds and first_divergence is None:
                first_divergence = old.ordinal
            classification = (
                "difference" if exceeds is None and any_difference else
                "identical" if exceeds is None else
                "divergent" if exceeds else "within-tolerance")
            if not math.isnan(time_difference) and time_difference != 0.0:
                classification += "+time-different"
            max_rel = max(item.max_rel for item in values.values())
            roundoff = max(item.roundoff_scale for item in values.values())
            print(
                f"{old.ordinal} {format_time(old.time)} {format_time(current.time)} "
                f"{time_difference:.17e} "
                f"{values['rho'].linf:.17e} {values['rhoE'].linf:.17e} "
                f"{values['p'].linf:.17e} {max_rel:.17e} {roundoff:.17e} {classification}")
            for semantic, item in values.items():
                print(
                    f"  {semantic}: count={item.count} maxAbs={item.linf:.17e} "
                    f"maxRel={item.max_rel:.17e} L1={item.l1:.17e} "
                    f"L2={item.l2:.17e} roundoffScale={item.roundoff_scale:.17e}")
        print(f"FIRST_DIFFERENT_STEP = {first_difference if first_difference is not None else 'none'}")
        print(f"FIRST_TIME_DIFFERENT_STEP = {first_time_difference if first_time_difference is not None else 'none'}")
        if args.abs_tol is not None:
            print(f"FIRST_DIVERGENCE_STEP = {first_divergence if first_divergence is not None else 'none'}")
        else:
            print("FIRST_DIVERGENCE_STEP = not-classified (provide --abs-tol and --rel-tol)")
        if len(old_frames) != len(current_frames):
            print(
                "SERIES_LENGTH_MISMATCH = "
                f"old={len(old_frames)}, current={len(current_frames)}, "
                f"compared=0..{compared_frames - 1}")
    except (OSError, ValueError, ET.ParseError) as error:
        print(f"compare_high_order_steps.py: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
