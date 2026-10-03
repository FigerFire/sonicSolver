#!/usr/bin/env python3
"""Frozen single-fluid SST arithmetic regression and compiled operation-count check."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def capture(case, log):
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    outputs = {}
    fields = {}
    boundary = {}
    for path in sorted((case / 'result').glob('*.vts'), key=lambda p: float(p.stem.rsplit('_t',1)[1])):
        if path.name.startswith('._'):
            continue
        outputs[path.name] = digest(path)
        root = ET.parse(path)
        piece = root.find('.//Piece')
        extent = list(map(int, piece.get('Extent').split()))
        nx, ny, nz = (extent[i + 1] - extent[i] + 1 for i in (0, 2, 4))
        arrays = {}
        edges = {}
        for array in root.findall('.//PointData/DataArray'):
            name = array.get('Name')
            values = list(map(float, (array.text or '').split()))
            if not values:
                continue
            if not all(math.isfinite(value) for value in values):
                raise AssertionError('Non-finite output: ' + name)
            def statistics(v):
                return {'min': min(v), 'max': max(v),
                        'l2': math.sqrt(sum(x*x for x in v)), 'sum': sum(v)}
            arrays[name] = statistics(values)
            nc = int(array.get('NumberOfComponents', '1'))
            edge = [values[((k*ny+j)*nx+i)*nc+c]
                    for k in range(nz) for j in range(ny) for i in range(nx)
                    if i in (0,nx-1) or j in (0,ny-1) for c in range(nc)]
            edges[name] = statistics(edge)
        fields[path.name] = arrays
        boundary[path.name] = edges
    inputs = {str(path.relative_to(case)): digest(path) for path in sorted(case.rglob('*'))
              if path.is_file() and 'result' not in path.relative_to(case).parts
              and path.name != 'baseline.json' and not path.name.startswith('._')}
    clocks = [list(pair) for pair in re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)', log)]
    checkpoints = [re.sub(r'0x[0-9a-fA-F]+', '<address>', line)
                   for line in log.splitlines() if 'State closure:' in line or '[SF TRACE]' in line]
    plan = [line for line in log.splitlines() if '[SF PLAN]' in line]
    operations = [re.search(r'op=([^\s]+)', line).group(1) for line in plan]
    return {'inputs': inputs, 'vts': outputs, 'fields': fields, 'boundary': boundary,
            'time_dt': clocks, 'diagnostics': checkpoints, 'plan': plan, 'operations': operations}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    baseline = json.loads((args.case / 'baseline.json').read_text())
    with tempfile.TemporaryDirectory(prefix='sonic-native-sst-') as directory:
        case = Path(directory) / 'case'
        shutil.copytree(args.case, case, ignore=shutil.ignore_patterns('result', '._*'))
        run = subprocess.run([str(args.solver.resolve()), 'run', str(case)], text=True,
                             capture_output=True, env={**os.environ, 'SF_PLAN_TRACE':'1',
                                                       'SF_HIGH_ORDER_TRACE':'1'})
        log = run.stdout + run.stderr
        if run.returncode:
            raise RuntimeError(log)
        actual = capture(case, log)
        if args.output:
            args.output.mkdir(parents=True, exist_ok=True)
            shutil.copytree(case,args.output/'case',dirs_exist_ok=True)
            (args.output/'run.log').write_text(log)
            (args.output/'evidence.json').write_text(json.dumps(actual,indent=2)+'\n')
        for key in ('inputs','vts','fields','boundary','time_dt','diagnostics'):
            if actual[key] != baseline[key]:
                raise AssertionError('First differing frozen evidence category: '+key)
        without_advance = [line for line in actual['plan'] if 'op=turbulence.advance' not in line]
        if without_advance != baseline['plan']:
            raise AssertionError('Flow plan clocks/stage timings changed')
        steps = len(actual['time_dt'])
        expected = ['flow.step.prepare','flow.dt.compute','turbulence.advance','flow.step.begin']
        expected += ['explicit.stage.execute']*4 + ['flow.step.commit','time.commit']
        if actual['operations'] != expected*steps:
            raise AssertionError('RAS count or physical-step order changed')
        first, last = next(iter(actual['fields'].values())), list(actual['fields'].values())[-1]
        for name in ('TurbulenceK','omega','TurbulentViscosity'):
            if first[name] == last[name]:
                raise AssertionError('Transport/closure publication did not change '+name)
    print('SST Passed identical: 6 steps, 4 VTS; 6 pair advances, 24 flow stages; clocks, fields, boundaries and checkpoints unchanged')


if __name__ == '__main__':
    main()
