#!/usr/bin/env python3
"""Run the bounded stationary SST/Ghost contract through default and authored HOW."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True, type=Path)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='sf-viscous-wall-') as temporary:
        folder = Path(temporary)
        def clone(name, target):
            case = folder / target
            shutil.copytree(args.root / 'test/IBM' / name, case,
                            ignore=shutil.ignore_patterns('result', '._*', '*.tar.gz', '*.vts', '*.vtm'))
            runtime = case / 'solvers/runtime.yaml'
            runtime.write_text(runtime.read_text().replace('writeInterval: 0.01', 'writeInterval: 0.00001'))
            return case
        def run(case, negative=None):
            env = dict(os.environ, SF_PLAN_TRACE='1')
            command = [str(args.solver), 'run', '--steps', '1' if negative else '4']
            result = subprocess.run(command + [str(case)], env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            if negative:
                if result.returncode == 0 or negative.lower() not in result.stdout.lower():
                    raise RuntimeError('Expected rejection: ' + negative + '\n' + result.stdout)
                print('Rejected:', negative)
                return
            if result.returncode:
                raise RuntimeError(result.stdout)
            times = re.findall(r'Step time: time=([^,]+), dt=([^\n]+)', result.stdout)
            trace = re.findall(r'\[SF PLAN\] (.*)', result.stdout)
            if len(times) != 4 or sum('op=turbulence.advance' in item for item in trace) != 4:
                raise RuntimeError('RAS did not execute exactly once per physical step')
            if sum('op=explicit.stage.execute' in item for item in trace) != 16:
                raise RuntimeError('SST composition changed RK4 stage count')
            hashes = [hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in sorted((case/'result').glob('*.vts')) if not path.name.startswith('._')]
            if len(hashes) != 5:
                raise RuntimeError('Missing step outputs for numerical schedule comparison')
            return times, trace, hashes
        default = clone('cylinderFlowSSTGhost', 'default')
        authored = clone('cylinderFlowSSTGhostAuthored', 'authored')
        if run(default) != run(authored):
            raise RuntimeError('Default/authored SST numerical results or operation contexts differ')
        print('Default/authored: identical dt, 16 stages, 4 RAS updates and five complete fields')
        slip = clone('cylinderFlowSSTGhost', 'slip')
        ibm = slip/'models/IBM.yaml'
        ibm.write_text(ibm.read_text().replace('stationaryNoSlipAdiabatic', 'eulerSlip'))
        run(slip, 'turbulence.boundary.immersed')
        ilw = clone('cylinderFlowSSTGhost', 'ilw')
        path = ilw/'models/ILW.yaml'
        path.write_text(path.read_text().replace('enabled: false', 'enabled: true').replace('value: 0', 'value: 5'))
        run(ilw, 'ILW')
        parallel = clone('cylinderFlowSSTGhost', 'parallel')
        path = parallel/'solvers/runtime.yaml'
        path.write_text(path.read_text().replace('enabled: false', 'enabled: true'))
        run(parallel, 'serial')
        multi = clone('cylinderFlowSSTGhost', 'multi')
        path = multi/'mesh/mesh.yaml'
        path.write_text(path.read_text().replace('files: [mesh/mesh.sfm]', 'files: [mesh/mesh.sfm, mesh/mesh.sfm]'))
        run(multi, 'patch')
        interleaved = clone('cylinderFlowSSTGhostAuthored', 'interleaved')
        path = interleaved/'solvers/how.yaml'
        how = json.loads(path.read_text())
        children = how['root']['children']
        children[1], children[2] = children[2], children[1]
        path.write_text(json.dumps(how, indent=2))
        run(interleaved, 'HOW placement violation')

if __name__ == '__main__':
    main()
