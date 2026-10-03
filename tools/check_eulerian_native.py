#!/usr/bin/env python3
"""Frozen existing Eulerian shared-pressure CFD evidence and operation-order regression."""
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
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    inputs = {str(p.relative_to(case)): digest(p) for p in sorted(case.rglob('*'))
              if p.is_file() and 'result' not in p.relative_to(case).parts
              and p.name != 'baseline.json' and not p.name.startswith('._')}
    fields, outputs, boundaries, reconstructed = {}, {}, {}, {}
    for path in sorted((case/'result').glob('*.vts')):
        if path.name.startswith('._'):
            continue
        outputs[path.name] = digest(path)
        root = ET.parse(path)
        extent = list(map(int, root.find('.//Piece').get('Extent').split()))
        nx, ny, nz = (extent[i+1]-extent[i]+1 for i in (0, 2, 4))
        arrays, edges, raw = {}, {}, {}
        for array in root.findall('.//PointData/DataArray'):
            values = list(map(float, (array.text or '').split()))
            if not values:
                continue
            assert all(math.isfinite(x) for x in values), array.get('Name')
            stats = lambda v: {'min': min(v), 'max': max(v),
                              'l2': math.sqrt(sum(x*x for x in v)), 'sum': sum(v)}
            name = array.get('Name')
            raw[name] = values
            arrays[name] = stats(values)
            nc = int(array.get('NumberOfComponents', '1'))
            edge = [values[((k*ny+j)*nx+i)*nc+c]
                    for k in range(nz) for j in range(ny) for i in range(nx)
                    if i in (0,nx-1) or j in (0,ny-1) for c in range(nc)]
            edges[name] = stats(edge)
        reconstructed[path.name] = {}
        for name, mass in raw.items():
            if not name.startswith('phaseMass.'):
                continue
            phase = name[len('phaseMass.'):]
            velocity, enthalpy = raw['U.'+phase], raw['phaseEnthalpy.'+phase]
            reconstructed[path.name]['momentum.'+phase] = stats([mass[i//3]*v for i,v in enumerate(velocity)])
            reconstructed[path.name]['h.'+phase] = stats([h/m for h,m in zip(enthalpy,mass)])
        fields[path.name], boundaries[path.name] = arrays, edges
    plan = [line for line in log.splitlines() if '[SF PLAN]' in line]
    diagnostics = [line for line in log.splitlines() if '[SF] Eulerian step:' in line]
    operations = [re.search(r'op=(\S+)', line).group(1) for line in plan]
    return {'inputs': inputs, 'vts': outputs, 'fields': fields, 'boundary': boundaries,
            'plan': plan, 'operations': operations, 'diagnostics': diagnostics, 'reconstructed_from_vts': reconstructed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    baseline = json.loads((args.case/'baseline.json').read_text())
    with tempfile.TemporaryDirectory(prefix='sonic-eulerian-native-') as directory:
        case = Path(directory)/'case'
        shutil.copytree(args.case, case, ignore=shutil.ignore_patterns('result', '._*'))
        run = subprocess.run([str(args.solver.resolve()), 'run', str(case)], text=True,
                             capture_output=True, env={**os.environ, 'SF_PLAN_TRACE':'1'})
        log = run.stdout+run.stderr
        if run.returncode:
            raise RuntimeError(log)
        actual = capture(case, log)
        if args.output:
            args.output.mkdir(parents=True, exist_ok=True)
            shutil.copytree(case,args.output/'case',dirs_exist_ok=True)
            (args.output/'run.log').write_text(log)
            (args.output/'evidence.json').write_text(json.dumps(actual,indent=2)+'\n')
        for key in ('inputs','vts','fields','boundary','plan','operations','diagnostics','reconstructed_from_vts'):
            if actual[key] != baseline[key]:
                raise AssertionError('First differing frozen evidence category: '+key)
        assert len(actual['diagnostics']) == 6 and len(actual['vts']) == 7
        assert len(actual['operations']) == 168
    print('Eulerian Passed identical: 6 steps; 7 VTS; 168 operations; phase/pressure fields, boundaries, clocks, closure and linear diagnostics')


if __name__ == '__main__':
    main()
