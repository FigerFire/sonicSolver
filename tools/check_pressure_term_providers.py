#!/usr/bin/env python3
"""Exercise compiled Momentum sources without changing pressure-plan topology."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile

from check_distributed_pressure import compare, fields, invoke, prepare


def configure(case: Path, source: str) -> None:
    numerics = case / 'solvers/numerics.yaml'
    numerics.write_text(numerics.read_text() + f'\nsources:\n  default: {source}\n')
    registry = case / 'models/models.yaml'
    registry.write_text(registry.read_text() +
                        f'{source}:\n  type: {source}\n  file: models/{source}.yaml\n')
    if source == 'gravity':
        contents = ('SonicFile:\n  object: models\n  type: gravity\n'
                    'value: [0.01, 0.0, 0.0]\n')
    elif source == 'MRF':
        contents = ('SonicFile:\n  object: models\n  type: MRF\nMRF:\n'
                    '  zone: all\n  center: [0, 0, 0]\n'
                    '  axis: [0, 0, 1]\n  omega: 0.1\n'
                    '  velocity: [0, 0, 0]\n')
    else:
        raise AssertionError(source)
    (case / f'models/{source}.yaml').write_text(contents)


def explain(solver: Path, case: Path) -> str:
    result = subprocess.run([str(solver), 'explain', str(case)],
                            capture_output=True, text=True, timeout=60)
    output = result.stdout + result.stderr
    if result.returncode:
        raise AssertionError(f'explain failed for {case.name}:\n{output[-3000:]}')
    return output


def plan(report: str) -> str:
    found = re.search(r'COMPILED SOLVE PLAN\s*\n(.*?)\nRUNTIME STATUS',
                      report, re.DOTALL)
    if not found:
        raise AssertionError('explain omitted the compiled plan')
    return found.group(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--mpi', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    solver = args.solver.resolve()
    directory = args.output or Path(tempfile.mkdtemp(prefix='sonic-pressure-terms-'))
    directory.mkdir(parents=True, exist_ok=True)

    for preset in ('PISO', 'SIMPLE', 'PIMPLE'):
        control = prepare(solver, args.case, directory, f'{preset}-control',
                          (1, 1, 1), preset, 2)
        reference = explain(solver, control)
        before = invoke(solver, control)
        for source, provider in (('gravity', 'source.gravity.primitive'),
                                 ('MRF', 'source.mrf.primitive')):
            case = prepare(solver, args.case, directory,
                           f'{preset}-{source}', (1, 1, 1), preset, 2)
            configure(case, source)
            report = explain(solver, case)
            if plan(report) != plan(reference):
                raise AssertionError(f'{preset}/{source} changed pressure-plan topology')
            if provider not in report or 'RUNTIME STATUS\n  runnable' not in report:
                raise AssertionError(f'{preset}/{source} provider is not runnable')
            log = invoke(solver, case)
            clocks = re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)', log)
            if len(clocks) != 2 or clocks != re.findall(
                    r'Step time: time=([^,\n]+), dt=([^\n]+)', before):
                raise AssertionError(f'{preset}/{source} changed physical clock')
            left, right = fields(control, 2), fields(case, 2)
            if left.keys() != right.keys():
                raise AssertionError(f'{preset}/{source} changed mesh')
            difference = max(abs(left[key][c] - right[key][c])
                             for key in left for c in range(3))
            if difference <= 1e-10 or not math.isfinite(difference):
                raise AssertionError(f'{preset}/{source} source was not effective')
            print(f'{preset}/{source}: source active; U delta={difference:.9g}; '
                  'clock and plan unchanged', flush=True)

    turbulent = prepare(solver, args.case, directory, 'PISO-turbulence',
                        (1, 1, 1), 'PISO', 2)
    settings = turbulent / 'models/turbulence.yaml'
    settings.write_text(settings.read_text().replace('enabled: false',
                                                    'enabled: true'))
    report = explain(solver, turbulent)
    if 'RUNTIME STATUS\n  unsupported' not in report or not all(
            item in report for item in ('closure.turbulence',
                                        'flow.turbulence', 'turbulence.advance')):
        raise AssertionError('pressure turbulence was reported runnable without '
                             'its transported-equation runtime')
    attempted = subprocess.run([str(solver), 'run', str(turbulent)],
                               capture_output=True, text=True, timeout=60)
    if attempted.returncode == 0:
        raise AssertionError('unsupported pressure turbulence executed silently')
    print('PISO/turbulence: explicitly Unsupported; no silent stress-only run',
          flush=True)

    if args.mpi:
        serial = directory / 'PISO-gravity'
        parallel = prepare(solver, args.case, directory, 'gravity-211',
                           (2, 1, 1), 'PISO', 2)
        configure(parallel, 'gravity')
        serial_log = invoke(solver, serial, trace=True)
        parallel_log = invoke(solver, parallel, 2, trace=True)
        du, dp = compare(serial, parallel, serial_log, parallel_log, 2)
        print(f'gravity MPI: matrix/RHS equivalent; U delta={du:.9g}; '
              f'p shape delta={dp:.9g}', flush=True)
    print(f'results: {directory}')


if __name__ == '__main__':
    main()
