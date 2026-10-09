#!/usr/bin/env python3
"""Exercise native interface RK/pseudo-time and authored HOW/WHICH on shared fixtures."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def output(path):
    return {str(p.relative_to(path / 'result')): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((path / 'result').rglob('*.vts')) if not p.name.startswith('._')}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--work-dir')
    args = parser.parse_args()
    solver, root = Path(args.solver).resolve(), Path(args.root).resolve()
    fixture = root / 'test/interface/nativeLevelSet'
    with tempfile.TemporaryDirectory(prefix='sonic-interface-', dir=args.work_dir) as temporary:
        workspace = Path(temporary)
        def make_case(name):
            case = workspace / name
            shutil.copytree(root / 'test/Sod/sodCase_weno7', case,
                ignore=shutil.ignore_patterns('result', '._*'))
            shutil.copyfile(fixture / 'mesh.sfm', case / 'mesh/mesh.sfm')
            for file, directory in [('multiPhase.yaml', 'models'), ('numerics.yaml', 'solvers'), ('runtime.yaml', 'solvers')]:
                shutil.copyfile(fixture / file, case / directory / file)
            with (case / 'models/models.yaml').open('a') as registry:
                registry.write('multiPhase:\n  type: multiPhase\n  file: models/multiPhase.yaml\n')
            return case
        def run(case, success=True, reason=None):
            result = subprocess.run([str(solver), 'run', str(case)], capture_output=True, text=True, timeout=90, env={**os.environ, 'SF_PLAN_TRACE': '1'})
            (case / 'run.log').write_text(result.stdout + result.stderr)
            if success != (result.returncode == 0):
                raise AssertionError(result.stdout[-2000:] + result.stderr)
            if reason and reason not in result.stdout + result.stderr:
                raise AssertionError('Failure did not report capability: ' + reason + '\n' + result.stdout + result.stderr)
            return re.findall(r'Step time: time=([^,]+), dt=([^\n]+)', result.stdout)
        reference = make_case('preset')
        clock = run(reference)
        assert len(clock) == 3, clock
        expected = json.loads((fixture / 'frozen-kernel.json').read_text())
        assert output(reference) == expected['vts_sha256'], 'Native preset differs from isolated old numerical kernel.'
        assert clock == [tuple(item) for item in expected['time_dt']], 'RK clock differs from old kernel.'
        authored = make_case('authored')
        children = [{'kind': 'EquationCall', 'equation': equation, 'target': target, 'order': order}
                    for equation, target, order in [('continuity', 'rho', 10), ('momentum', 'rhoU', 20), ('energy', 'rhoE', 60), ('levelSetAdvection', 'phi', 70)]]
        children += [{'kind': 'Sequence', 'id': 'interface.finalize', 'order': 990000, 'children': [
            {'kind': 'EquationCall', 'equation': 'levelSetReference', 'target': 'phi0', 'targetKind': 'Workspace'},
            {'kind': 'Loop', 'id': 'interface.pseudoTime', 'repetitions': 3, 'children': [
                {'kind': 'EquationCall', 'equation': 'levelSetReinitialization', 'target': 'phi'}]},
            {'kind': 'EquationCall', 'equation': 'interfaceGeometry', 'target': 'interfaceCurvature'}]},
            {'kind': 'Commit', 'id': 'physicalStep.commit', 'order': 1000000}]
        bindings = [{'equation': equation, 'method': 'ConservativeResidual'} for equation in ['continuity', 'momentum', 'energy']]
        params = {'order': 5, 'epsilon': 1e-6, 'power': 2, 'pseudoDt': .01, 'signFactor': 1}
        bindings += [{'equation': 'levelSetAdvection', 'method': 'LevelSetAdvection', 'parameters': {'order': 5, 'epsilon': 1e-6, 'power': 2, 'csf': 0, 'ghostFluid': 0, 'sigma': 0, 'width': 0}},
                     {'equation': 'levelSetReference', 'method': 'LevelSetReference', 'parameters': params},
                     {'equation': 'levelSetReinitialization', 'method': 'LevelSetReinitialization', 'inputs': ['levelSet.phi0']},
                     {'equation': 'interfaceGeometry', 'method': 'LevelSetGeometry'}]
        def write_how(case, nodes=children, choices=bindings):
            (case / 'solvers/how.yaml').write_text(json.dumps({'SonicFile': {'object': 'solver', 'type': 'executionProgram'}, 'root': {'kind': 'Sequence', 'children': nodes}}))
            (case / 'solvers/which.yaml').write_text(json.dumps({'SonicFile': {'object': 'solver', 'type': 'providerBindings'}, 'bindings': choices}))
            with (case / 'solvers/solvers.yaml').open('a') as registry:
                registry.write('how:\n  type: executionProgram\n  file: solvers/how.yaml\nwhich:\n  type: providerBindings\n  file: solvers/which.yaml\n')
        write_how(authored)
        assert run(authored) == clock
        assert output(authored) == output(reference), 'Authored and preset contracts produce different fields.'
        def operations(case):
            return re.findall(r'\[SF PLAN\].*?op=(\S+)', (case / 'run.log').read_text())
        assert operations(authored) == operations(reference), 'Authored and preset compiled operations differ.'
        assert operations(authored).count('time.commit') == 3, 'Pseudo-time gained physical-clock authority.'

        for name, change, reason in [
            ('unknown', lambda c, b: b[0].update(method='UnimplementedProvider'), 'UnimplementedProvider'),
            ('badEquation', lambda c, b: (c[0].update(equation='unknownEquation'), b[0].update(equation='unknownEquation')), 'unknownEquation'),
            ('missingWhich', lambda c, b: b.clear(), 'WHICH'),
            ('badTarget', lambda c, b: c[0].update(target='unknownState'), 'unknownState'),
            ('badTargetKind', lambda c, b: c[0].update(targetKind='InvalidKind'), 'InvalidKind'),
            ('reversedHOW', lambda c, b: c.reverse(), 'placement violation'),
            ('emptyLoop', lambda c, b: c[4]['children'][1].update(repetitions=0), 'Loop'),
            ('unusedParameter', lambda c, b: b[3]['parameters'].update(unknownOption=1), 'unknownOption'),
            ('hiddenPhi', lambda c, b: (c.pop(3), b.pop(3)), 'levelSetAdvection'),
        ]:
            case = make_case(name)
            nodes, choices = json.loads(json.dumps(children)), json.loads(json.dumps(bindings))
            change(nodes, choices);write_how(case, nodes, choices);run(case, False, reason)
        print('Native preset/explicit HOW+WHICH match frozen kernel; unsupported provider/equation/target and hidden scalar paths fail visibly.')


if __name__ == '__main__':
    main()
