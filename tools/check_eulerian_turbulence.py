#!/usr/bin/env python3
"""Native Eulerian RAS against frozen pre-migration kernel/wiring evidence."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from check_eulerian_native import capture


def compare(actual, baseline):
    for key in ('inputs', 'plan', 'operations', 'diagnostics', 'vts', 'fields',
                'boundary', 'reconstructed_from_vts'):
        if actual[key] != baseline[key]:
            # Report the first step/array when possible; never adjust a tolerance.
            if isinstance(actual[key], list):
                for i, (a, b) in enumerate(zip(actual[key], baseline[key])):
                    if a != b:
                        raise AssertionError(f'First differing {key}[{i}]: {a!r} != {b!r}')
            if isinstance(actual[key], dict):
                for name in sorted(set(actual[key]) | set(baseline[key])):
                    if actual[key].get(name) != baseline[key].get(name):
                        raise AssertionError(f'First differing {key}: {name}')
            raise AssertionError('First differing frozen evidence category: '+key)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    baseline = json.loads((args.case/'baseline.json').read_text())
    assert baseline['provenance']['is_historical_production_baseline'] is False
    with tempfile.TemporaryDirectory(prefix='sonic-eulerian-turbulence-') as directory:
        case = Path(directory)/'case'
        shutil.copytree(args.case, case, ignore=shutil.ignore_patterns('result', '._*'))
        run = subprocess.run([str(args.solver.resolve()), 'run', str(case)], text=True,
                             capture_output=True, env={**os.environ, 'SF_PLAN_TRACE': '1'})
        log = run.stdout+run.stderr
        if run.returncode:
            raise RuntimeError(log)
        actual = capture(case, log)
        if args.output:
            args.output.mkdir(parents=True, exist_ok=True)
            shutil.copytree(case, args.output/'case', dirs_exist_ok=True)
            (args.output/'run.log').write_text(log)
            (args.output/'evidence.json').write_text(json.dumps(actual, indent=2)+'\n')
        compare(actual, baseline)
        steps = len(actual['diagnostics'])
        assert steps == 6 and len(actual['vts']) == 7
        assert actual['operations'].count('ee.turbulence.prepare') == steps
        assert actual['operations'].count('ee.turbulence.solve') == steps
        first, last = next(iter(actual['fields'].values())), list(actual['fields'].values())[-1]
        transported = [name for name in first if name.startswith(('k.', 'omega.', 'epsilon.'))]
        assert len(transported) == 4
        for name in transported:
            assert first[name] != last[name], 'Transport did not advance '+name
    print(f'{args.case.name}: Passed identical; 6 steps, 7 VTS, 180 operations; phase/RAS fields, boundary, dt/time and linear diagnostics')


if __name__ == '__main__':
    main()
