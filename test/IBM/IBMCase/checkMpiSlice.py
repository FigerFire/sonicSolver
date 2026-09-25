from collections import defaultdict
import math
import os
import sys

try:
    from paraview import servermanager
    from paraview.simple import OpenDataFile, Slice
except Exception:
    servermanager = None
    OpenDataFile = None
    Slice = None


def iter_datasets(data):
    if data is None:
        return
    if data.IsA("vtkCompositeDataSet"):
        iterator = data.NewIterator()
        iterator.SkipEmptyNodesOn()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            child = iterator.GetCurrentDataObject()
            if child is not None:
                yield child
            iterator.GoToNextItem()
        return
    yield data


def tuple_value(array, index):
    if array.GetNumberOfComponents() == 1:
        return (array.GetTuple1(index),)
    return tuple(array.GetTuple(index))


def fields_to_compare():
    return [
        "Density",
        "Pressure",
        "Velocity",
        "BC_ID",
        "BC_ID_rho",
        "BC_ID_U",
        "BC_ID_p",
    ]


def field_scalar(name, value):
    if name == "Velocity":
        return math.sqrt(sum(component * component for component in value))
    return value[0]


def free_stream_thresholds():
    return {
        "Density": float(os.environ.get(
            "SLICE_MAX_DENSITY_JUMP", "1.0e-10")),
        "Pressure": float(os.environ.get(
            "SLICE_MAX_PRESSURE_JUMP", "1.0e-7")),
        "Velocity": float(os.environ.get(
            "SLICE_MAX_VELOCITY_JUMP", "1.0e-10")),
    }


def duplicate_thresholds():
    return {
        "Density": 1.0e-13,
        "Pressure": 1.0e-10,
        "Velocity": 1.0e-13,
        "BC_ID": 0.0,
        "BC_ID_rho": 0.0,
        "BC_ID_U": 0.0,
        "BC_ID_p": 0.0,
    }


def build_unique_points(groups):
    unique = []
    for key, samples in groups.items():
        point = key[:3]
        values = {}
        for name in ("Density", "Pressure", "Velocity"):
            field_values = [sample[name] for sample in samples if name in sample]
            if not field_values:
                continue
            # 求解器契约是 GlobalDof owner-to-copy；诊断也直接采用 canonical
            # 副本，不用平均掩盖副本不一致。严格一致性由
            # compare_grouped_points() 单独检查。
            values[name] = field_values[0]
        if values:
            unique.append((point, values))
    return unique


def nearest_point(points, target, predicate):
    tx, ty = target
    best = None
    best_distance = math.inf
    for point, values in points:
        x, y, _ = point
        if not predicate(x, y):
            continue
        distance = (x - tx) ** 2 + (y - ty) ** 2
        if distance < best_distance:
            best_distance = distance
            best = (point, values, math.sqrt(distance))
    return best


def compare_interface_neighbour_jumps(groups):
    points = build_unique_points(groups)
    if not points:
        print("interfaceNeighbourJump=skipped")
        return 1

    thresholds = free_stream_thresholds()
    seams = [
        ("centerTop", "h", 0.04, 1.0),
        ("centerBot", "h", -0.04, -1.0),
        ("centerLeft", "v", -0.04, -1.0),
        ("centerRight", "v", 0.04, 1.0),
    ]
    spacing = 0.002
    failed = False

    for seam_name, kind, coordinate, sign in seams:
        if kind == "h":
            seam_points = [
                (point, values)
                for point, values in points
                if abs(point[1] - coordinate) < 5.0e-10
                and -0.039 < point[0] < 0.039
            ]
        else:
            seam_points = [
                (point, values)
                for point, values in points
                if abs(point[0] - coordinate) < 5.0e-10
                and -0.039 < point[1] < 0.039
            ]

        print(f"interfaceSamples[{seam_name}]={len(seam_points)}")
        for field in ("Density", "Pressure", "Velocity"):
            max_jump = 0.0
            max_seam_deviation = 0.0
            max_deviation_sample = None
            for point, values in seam_points:
                x, y, _ = point
                if kind == "h":
                    inner = nearest_point(
                        points,
                        (x, y - sign * spacing),
                        lambda xx, yy: abs(xx - x) < 4.0e-4
                        and (yy - coordinate) * sign < -1.0e-9,
                    )
                    outer = nearest_point(
                        points,
                        (x, y + sign * spacing),
                        lambda xx, yy: abs(xx - x) < 4.0e-4
                        and (yy - coordinate) * sign > 1.0e-9,
                    )
                else:
                    inner = nearest_point(
                        points,
                        (x - sign * spacing, y),
                        lambda xx, yy: abs(yy - y) < 4.0e-4
                        and (xx - coordinate) * sign < -1.0e-9,
                    )
                    outer = nearest_point(
                        points,
                        (x + sign * spacing, y),
                        lambda xx, yy: abs(yy - y) < 4.0e-4
                        and (xx - coordinate) * sign > 1.0e-9,
                    )
                if inner is None or outer is None:
                    continue
                seam_value = field_scalar(field, values[field])
                inner_value = field_scalar(field, inner[1][field])
                outer_value = field_scalar(field, outer[1][field])
                jump = abs(outer_value - inner_value)
                seam_deviation = abs(seam_value - 0.5 * (inner_value + outer_value))
                max_jump = max(max_jump, jump)
                if seam_deviation > max_seam_deviation:
                    max_seam_deviation = seam_deviation
                    max_deviation_sample = (
                        point, seam_value, inner_value, outer_value,
                        inner[0], outer[0],
                    )

            print(
                f"interfaceNeighbourJump[{seam_name}][{field}]="
                f"{max_jump:.12e}"
            )
            print(
                f"interfaceSeamDeviation[{seam_name}][{field}]="
                f"{max_seam_deviation:.12e}"
            )
            if max_deviation_sample is not None:
                print(
                    f"interfaceSeamDeviationSample[{seam_name}][{field}]="
                    f"{max_deviation_sample}"
                )
            if max_jump > thresholds[field] \
                    or max_seam_deviation > thresholds[field]:
                failed = True

    return 1 if failed else 0


def compare_z_partition_jumps(groups):
    points = build_unique_points(groups)
    if not points:
        print("mpiPartitionJump=skipped")
        return 1

    by_xy = {
        (round(point[0], 10), round(point[1], 10), round(point[2], 10)): values
        for point, values in points
    }
    z0 = 0.095
    dz = 0.001
    samples = [
        (point, values)
        for point, values in points
        if abs(point[2] - z0) < 5.0e-10
    ]
    print(f"mpiPartitionSamples[z={z0}]={len(samples)}")

    thresholds = free_stream_thresholds()
    failed = False
    for field in ("Density", "Pressure", "Velocity"):
        max_jump = 0.0
        max_seam_deviation = 0.0
        for point, values in samples:
            x, y, _ = point
            left = by_xy.get((round(x, 10), round(y, 10), round(z0 - dz, 10)))
            right = by_xy.get((round(x, 10), round(y, 10), round(z0 + dz, 10)))
            if left is None or right is None:
                continue
            seam_value = field_scalar(field, values[field])
            left_value = field_scalar(field, left[field])
            right_value = field_scalar(field, right[field])
            max_jump = max(max_jump, abs(right_value - left_value))
            max_seam_deviation = max(
                max_seam_deviation,
                abs(seam_value - 0.5 * (left_value + right_value)),
            )

        print(f"mpiPartitionJump[z][{field}]={max_jump:.12e}")
        print(f"mpiPartitionSeamDeviation[z][{field}]={max_seam_deviation:.12e}")
        if max_jump > thresholds[field] \
                or max_seam_deviation > thresholds[field]:
            failed = True

    return 1 if failed else 0


def compare_grouped_points(groups, point_count, source_name):
    fields = fields_to_compare()
    tolerances = duplicate_thresholds()
    shared_groups = 0
    inconsistent_groups = 0
    max_difference = {name: 0.0 for name in fields}
    for samples in groups.values():
        if len(samples) < 2:
            continue
        shared_groups += 1
        for name in fields:
            values = [sample[name] for sample in samples if name in sample]
            if len(values) < 2:
                continue
            reference = values[0]
            difference = max(
                max(abs(a - b) for a, b in zip(reference, value))
                for value in values[1:]
            )
            max_difference[name] = max(max_difference[name], difference)
            if difference > tolerances[name]:
                inconsistent_groups += 1
                break

    print(f"slice={source_name}")
    print(f"slicePoints={point_count}")
    print(f"uniqueCoordinates={len(groups)}")
    print(f"sharedCoordinateGroups={shared_groups}")
    print(f"inconsistentGroups={inconsistent_groups}")
    for name in fields:
        if math.isfinite(max_difference[name]):
            print(f"maxDiff[{name}]={max_difference[name]:.12e}")

    expected = {
        "Density": 5.0,
        "Pressure": 10000.0,
        "Velocity": 0.0,
    }
    free_stream_error = {name: 0.0 for name in expected}
    for samples in groups.values():
        for sample in samples:
            # STL 内部诊断点在输出中保留零状态；自由流回归只评价具有
            # 正密度、正压力的 Eulerian 流体/IBM ghost 点。
            if "Density" not in sample or "Pressure" not in sample \
                    or sample["Density"][0] <= 0.0 \
                    or sample["Pressure"][0] <= 0.0:
                continue
            for name, reference in expected.items():
                if name not in sample:
                    continue
                free_stream_error[name] = max(
                    free_stream_error[name],
                    abs(field_scalar(name, sample[name]) - reference),
                )
    free_stream_status = 0
    free_stream_limits = free_stream_thresholds()
    for name, error in free_stream_error.items():
        print(f"freeStreamError[{name}]={error:.12e}")
        if error > free_stream_limits[name]:
            free_stream_status = 1

    duplicate_status = 0 if shared_groups > 0 and inconsistent_groups == 0 else 1
    interface_status = compare_interface_neighbour_jumps(groups)
    return 0 if duplicate_status == 0 \
        and interface_status == 0 \
        and free_stream_status == 0 else 1



def check_with_paraview(vtm):
    source = OpenDataFile(vtm)
    cut = Slice(Input=source)
    cut.SliceType = "Plane"
    cut.SliceType.Origin = [0.0, 0.0, 0.095]
    cut.SliceType.Normal = [0.0, 0.0, 1.0]
    cut.UpdatePipeline()
    result = servermanager.Fetch(cut)

    fields = fields_to_compare()
    groups = defaultdict(list)
    point_count = 0
    for dataset in iter_datasets(result):
        point_data = dataset.GetPointData()
        arrays = {name: point_data.GetArray(name) for name in fields}
        for point_id in range(dataset.GetNumberOfPoints()):
            point = dataset.GetPoint(point_id)
            key = tuple(round(value, 10) for value in point)
            values = {
                name: tuple_value(array, point_id)
                for name, array in arrays.items()
                if array is not None
            }
            groups[key].append(values)
            point_count += 1

    return compare_grouped_points(groups, point_count, vtm)


def check_with_vtk(vtm):
    import vtk

    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName(vtm)
    reader.Update()
    data = reader.GetOutput()

    fields = fields_to_compare()
    groups = defaultdict(list)
    volume_groups = defaultdict(list)
    point_count = 0
    target_z = 0.095
    tolerance = 1.0e-9

    for dataset in iter_datasets(data):
        point_data = dataset.GetPointData()
        arrays = {name: point_data.GetArray(name) for name in fields}
        for point_id in range(dataset.GetNumberOfPoints()):
            point = dataset.GetPoint(point_id)
            key = tuple(round(value, 10) for value in point)
            values = {
                name: tuple_value(array, point_id)
                for name, array in arrays.items()
                if array is not None
            }
            # z 分核检查只需要 z0-dz、z0、z0+dz 三层。不要把近百万
            # 全场点复制进 Python 字典，否则本地等价检查会无谓耗尽内存。
            if abs(abs(point[2] - target_z) - 0.001) <= tolerance \
                    or abs(point[2] - target_z) <= tolerance:
                volume_groups[key].append(values)
            if abs(point[2] - target_z) > tolerance:
                continue
            groups[key].append(values)
            point_count += 1

    slice_status = compare_grouped_points(groups, point_count, vtm)
    mpi_status = compare_z_partition_jumps(volume_groups)
    return 0 if slice_status == 0 and mpi_status == 0 else 1


def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    vtm = (
        sys.argv[1]
        if len(sys.argv) > 1
        else os.path.join(case_dir, "result", "IBMCylinder_000005.vtm")
    )

    # IBMCase 的 z=0.095 本来就是网格点层。优先直接读取该层，避免 ParaView
    # Slice 对流体点和零值固体诊断点做可视化插值，制造 rho/2、p/2 这类并非
    # 求解状态的中间值。只有运行环境没有 VTK reader 时才采用 Slice 后备路径。
    try:
        return check_with_vtk(vtm)
    except ImportError:
        if OpenDataFile is not None and Slice is not None:
            return check_with_paraview(vtm)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
