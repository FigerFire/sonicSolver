#!/usr/bin/env python3
"""Compare one physical pressure system under serial and MPI execution.

Thresholds are fixed before the run. Pressure is compared modulo the existing
per-inactive-layer nullspace; configured gauge rows and reference-layer pressure are also checked. Matrix rows are keyed by physical point, never backend row.
"""
import argparse
import math
import os
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"test"))
from thermophysical_fixture import declare_historical_thermophysical_fixture
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

U_ATOL = 1e-8
P_ATOL = 1e-6
MATRIX_ATOL = 1e-12


def invoke(solver, case, ranks=1, trace=False):
    command = [str(solver), 'run', str(case)]
    if ranks > 1:
        command = ['mpirun', '-np', str(ranks), *command]
    env = os.environ.copy()
    env['SF_PLAN_TRACE'] = '1'
    if trace:
        env['SF_PRESSURE_TRACE_DIR'] = str(case / 'result/assembly')
    try:
        run = subprocess.run(command, capture_output=True, text=True,
                             timeout=180, env=env)
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(f'{case.name}: timeout (possible collective mismatch)') from exc
    log = run.stdout + run.stderr
    (case / 'run.log').write_text(log)
    if run.returncode:
        if any(token in log for token in ('No network interfaces', 'Unable to start',
                                         'Operation not permitted', 'PRTE ERROR')):
            raise RuntimeError('MPI runtime launch unavailable:\n' + log[-3000:])
        raise AssertionError(f'{case.name} numerical/runtime failure:\n{log[-5000:]}')
    return log


def prepare(solver, source, parent, name, split, preset='PISO', steps=2,
            reference=0, cells=16):
    case = parent / name
    shutil.copytree(source, case, ignore=shutil.ignore_patterns(
        'result', '._*', '.DS_Store', '*.sfm', '*.vts', '*.vtm'))
    declare_historical_thermophysical_fixture(case)
    mesh = case / 'mesh/blockMeshDict'
    text = mesh.read_text()
    text = re.sub(r'\((?:80 80|16 16) 1\)', f'({cells} {cells} 1)', text)
    mesh.write_text(text)
    algorithm = case / 'solvers/algorithm.yaml'
    text = algorithm.read_text().replace('algorithm: PISO', f'algorithm: {preset}').replace('PISO:', preset + ':')
    text = text.replace('referenceCell: 0', f'referenceCell: {reference}')
    if preset != 'PISO':
        text = text.replace('outerCorrectors: 1', 'outerCorrectors: 3')
        text = text.replace('momentumRelaxation: 1.0', 'momentumRelaxation: 0.8')
        text = text.replace('pressureRelaxation: 1.0', 'pressureRelaxation: 0.8')
        text = text.replace('  pressureRelaxation: 0.8', '  pressureRelaxation: 0.8\n  relativeTolerance: 0.001\n  absoluteTolerance: 1e-8')
        if preset == 'SIMPLE':
            text = text.replace('pressureCorrectors: 2', 'pressureCorrectors: 1')
    algorithm.write_text(text)
    semantic=case/'algorithms/algorithms.yaml'
    if semantic.exists():
        semantic.write_text('SonicFile:\n  object: algorithms\n  type: registry\n'+preset+': {}\n')
    if name.startswith('equilibrium'):
        boundary = case / 'fields/boundaries.yaml'
        boundary.write_text(boundary.read_text().replace('[1, 0, 0]', '[0, 0, 0]'))
    runtime = case / 'solvers/runtime.yaml'
    text = runtime.read_text()
    text = re.sub(r'endStep:.*', f'endStep: {steps}', text)
    if steps<=8:
        text=re.sub(r'endTime:.*','endTime: 1.0',text)
    text = re.sub(r'writeInterval:.*', f'writeInterval: {steps}', text)
    # Same grid, initial condition, physical timestep controls for every partition.
    text = text.replace('createMesh: false', 'createMesh: true')
    runtime.write_text(text)
    invoke(solver, case)
    text = text.replace('createMesh: true', 'createMesh: false')
    ranks = math.prod(split)
    if ranks > 1:
        text = text.replace('enabled: false', 'enabled: true')
        text = re.sub(r'nProcs:.*', f'nProcs: {ranks}', text)
        text = re.sub(r'split:.*', 'split: [' + ', '.join(map(str, split)) + ']', text)
    runtime.write_text(text)
    return case


def fields(case, step):
    result = {}
    for path in (case / 'result').rglob('*.vts'):
        if path.name.startswith('._') or f'_{step:06d}' not in path.name:
            continue
        root = ET.parse(path).getroot()
        coords = list(map(float, root.find('.//Points/DataArray').text.split()))
        u = list(map(float, root.find('.//PointData/DataArray[@Name="U"]').text.split()))
        p = list(map(float, root.find('.//PointData/DataArray[@Name="p"]').text.split()))
        for index, pressure in enumerate(p):
            key = tuple(round(v, 12) for v in coords[3*index:3*index+3])
            value = (*u[3*index:3*index+3], pressure)
            if not all(math.isfinite(v) for v in value):
                raise AssertionError(f'nonfinite output at {key}')
            if key in result and max(abs(a-b) for a,b in zip(result[key],value)) > 1e-12:
                raise AssertionError(f'owner COPY mismatch at shared point {key}')
            result[key] = value
    if not result:
        raise AssertionError(f'{case}: no final VTS')
    return result


def matrices(case, index):
    rows = {}
    for path in (case / 'result/assembly').glob(f'matrix-{index}-patch-*.txt'):
        current = None
        for line in path.read_text().splitlines():
            tokens = line.split()
            key = tuple(round(float(v),12) for v in tokens[1:4])
            value = float(tokens[4])
            if tokens[0] == 'row':
                if key in rows:
                    raise AssertionError(f'duplicate physical pressure row {key}')
                rows[key] = (value,{})
                current = key
            else:
                rows[current][1][key] = value
    if not rows:
        raise AssertionError('missing pressure assembly trace')
    return rows


def compare(serial, parallel, before, after, steps):
    clock = lambda log: re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)', log)
    if clock(before) != clock(after) or len(clock(after)) != steps:
        raise AssertionError('global dt/physical clock differs')
    a,b = fields(serial,steps), fields(parallel,steps)
    if a.keys() != b.keys():
        raise AssertionError('physical output mesh differs')
    du=max(abs(a[k][c]-b[k][c]) for k in a for c in range(3))
    if du > U_ATOL:
        raise AssertionError(f'first final velocity difference {du} > {U_ATOL}')
    pressure_span=0.0
    for z in sorted({k[2] for k in a}):
        delta=[b[k][3]-a[k][3] for k in a if k[2] == z]
        span=max(delta)-min(delta)
        pressure_span=max(pressure_span,span)
        if span>P_ATOL:
            raise AssertionError(f'pressure shape differs on plane z={z}: {span}')
    for index in [0]:  # Initial matrix/RHS precedes differences in iterative linear solves.
        ma,mb=matrices(serial,index),matrices(parallel,index)
        if ma.keys()!=mb.keys():
            raise AssertionError('owner row identities differ')
        gauges=[key for key,(rhs,entries) in ma.items()
                if rhs == 0.0 and entries == {key: 1.0}]
        if parallel.name == 'cavity-121':
            owners=[]
            for path in (parallel / 'result/assembly').glob('matrix-0-patch-*.txt'):
                for line in path.read_text().splitlines():
                    tokens=line.split()
                    if tokens and tokens[0]=='row' and tuple(map(float,tokens[1:4])) in gauges:
                        owners.append(int(path.stem.rsplit('-',1)[1]))
            if len(owners)!=1 or owners[0]==0:
                raise AssertionError('reference-row regression did not exercise a nonzero owner')
        for key in gauges:
            if abs(a[key][3]-b[key][3]) > P_ATOL:
                raise AssertionError('configured reference pressure differs')
            if max(abs(a[k][3]-b[k][3]) for k in a if k[2] == key[2]) > P_ATOL:
                raise AssertionError('reference-layer absolute pressure differs')
        for key,(rhs,entries) in ma.items():
            rb,eb=mb[key]
            if entries.keys()!=eb.keys():
                raise AssertionError(f'first matrix sparsity difference at {key}')
            if abs(rhs-rb)>MATRIX_ATOL or any(abs(v-eb[col])>MATRIX_ATOL for col,v in entries.items()):
                raise AssertionError(f'first matrix/RHS difference at {key}')
    serial_iterations=before.count('op=pressure.iteration.begin')
    ranks=math.prod(tuple(map(int, parallel.name.rsplit('-',1)[1])))
    if after.count('op=pressure.iteration.begin') != ranks*serial_iterations:
        raise AssertionError('rank/serial outer-iteration control differs')
    if 'equilibrium' in serial.name and serial_iterations != steps:
        raise AssertionError('fixed-point regression never exercised early exit')
    if 'channel' in serial.name:
        errors=[a[k][0] - 0.5*(0.25-k[1]*k[1]) for k in a
                if abs(k[0])<1e-10 and -0.5<k[1]<0.5]
        l2=math.sqrt(sum(v*v for v in errors)/len(errors))
        if l2>1.2e-3:
            raise AssertionError(f'Poiseuille analytical profile failed: {l2}')
    return du,pressure_span


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--solver',type=Path,required=True)
    parser.add_argument('--cavity',type=Path,required=True)
    parser.add_argument('--poiseuille',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    directory=args.output or Path(tempfile.mkdtemp(prefix='sonic-distributed-pressure-'))
    directory.mkdir(parents=True,exist_ok=True)
    solver=args.solver.resolve()
    for label,source,preset,steps in [
            ('cavity',args.cavity,'PISO',2),
            ('channel',args.poiseuille,'PISO',2048),
            ('simple',args.poiseuille,'SIMPLE',8),
            ('pimple',args.poiseuille,'PIMPLE',8),
            ('equilibrium-simple',args.cavity,'SIMPLE',2),
            ('equilibrium-pimple',args.cavity,'PIMPLE',2)]:
        reference=prepare(solver,source,directory,label+'-serial',(1,1,1),preset,steps,reference=150)
        before=invoke(solver,reference,trace=True)
        for split in ([(2,1,1),(1,2,1),(2,2,1)] if preset=='PISO' else [(2,1,1)]):
            name=label+'-'+''.join(map(str,split))
            case=prepare(solver,source,directory,name,split,preset,steps,reference=150)
            after=invoke(solver,case,math.prod(split),trace=True)
            du,dp=compare(reference,case,before,after,steps)
            print(f'{name}: same dt; matrix/RHS equivalent; U maxDelta={du:.9g}; p shape delta={dp:.9g}',flush=True)
    print('results:',directory)

if __name__=='__main__':
    main()
