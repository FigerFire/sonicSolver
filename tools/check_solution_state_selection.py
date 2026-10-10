#!/usr/bin/env python3
"""Exercise explicit WHAT/STATE/HOW through the production native reader."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"test"))
from thermophysical_fixture import declare_historical_thermophysical_fixture
import shutil
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--solver',required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--output')
    args=parser.parse_args()
    records=[]
    with tempfile.TemporaryDirectory(prefix='sonic-solution-state-') as directory:
        root=Path(directory)
        def case(name,states=None,algorithm=None,eos=None,equations=None):
            target=root/name
            shutil.copytree(args.case,target,ignore=shutil.ignore_patterns('result','._*','.DS_Store'))
            declare_historical_thermophysical_fixture(target)
            if states is not None:
                (target/'state/state.yaml').write_text('SonicFile:\n  object: state\n  type: registry\nuse: ['+', '.join(states)+']\n')
            if algorithm is not None:
                (target/'algorithms/algorithms.yaml').write_text('SonicFile:\n  object: algorithms\n  type: registry\n'+algorithm+(': {pressureCorrectors: 1}\n' if algorithm=='SIMPLE' else ': {}\n'))
            if eos is not None:
                eos_properties=('    rho: 1.0\n' if eos=='rhoConst' else '    gamma: 1.4\n    R: 287.05\n')
                (target/'models/thermoDynamics.yaml').write_text('SonicFile:\n  object: models\n  type: thermoDynamics\nthermoDynamics:\n  equationOfState: '+eos+'\n  thermo: hConst\n  transport: const\nproperties:\n  equationOfState:\n'+eos_properties+'  transport:\n    mu: 0.01\n    Pr: 0.72\n')
            if equations is not None:
                (target/'equations/equations.yaml').write_text('SonicFile:\n  object: equations\n  type: registry\nuse: ['+', '.join(equations)+']\n')
            return target
        def check(name,target,expected=None):
            result=subprocess.run([args.solver,'explain',str(target)],text=True,capture_output=True)
            text=result.stdout+result.stderr
            records.append({'name':name,'exit_code':result.returncode,'output':text})
            if expected is None:
                if result.returncode or 'use : U p' not in text or 'require :' not in text:
                    raise AssertionError(name+': '+text)
            elif result.returncode==0 or expected not in text:
                raise AssertionError(name+' did not reject the requested contract: '+text)
        # A: no solver-family identity participates in native parsing.
        positive=case('constant-simple',algorithm='SIMPLE')
        old=positive/'solvers/algorithm.yaml'
        old.write_text(old.read_text().replace('type: pressureBase\n',''))
        check('A rhoConst + U,p + SIMPLE',positive)
        check('B rhoConst + conservative + Explicit',
              case('constant-conservative',['rho','rhoU','rhoE'],'Explicit',
                   eos='rhoConst', equations=['Continuity','Momentum','Energy']),
              'closure-consistent density evolution')
        check('C perfectGas + U,p + PISO',case('variable-primitive',eos='perfectGas'),
              'variable-density primitive momentum/continuity')
        check('D ambiguous momentum',case('ambiguous',['U','p','rhoU']),
              'Ambiguous Momentum target')
        missing=case('missing-state');(missing/'state/state.yaml').unlink()
        check('E missing STATE',missing,'requires explicit solution STATE')
        missing=case('missing-how');(missing/'algorithms/algorithms.yaml').unlink()
        check('missing HOW',missing,'requires explicit HOW')
        check('duplicate STATE',case('duplicate',['U','p','U']),'Duplicate solution STATE')
        check('unowned solution',case('unowned',['U','p','T']),
              "No compiled HOW output owns selected solution STATE 'T'")
        caloric=case('missing-caloric',['rho','rhoU','rhoE'],'Explicit','perfectGas',
                     ['Continuity','Momentum','Energy'])
        model=caloric/'models/thermoDynamics.yaml'
        model.write_text(model.read_text().replace('thermo: hConst','thermo: janaf'))
        check('unimplemented caloric closure',caloric,"Unsupported caloric provider 'janaf'")
        transport=case('missing-transport')
        model=transport/'models/thermoDynamics.yaml'
        model.write_text(model.read_text().replace('transport: const','transport: sutherland'))
        check('unimplemented transport closure',transport,"Unsupported transport provider 'sutherland'")
        # Canonical physical IDs remain unchanged; unknown catalogue metadata is not guessed.
        check('unknown STATE',case('unknown',['U','p','custom']),
              'Solution STATE has no metadata: custom')
    if args.output:
        output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(records,indent=2)+'\n')
    print('Explicit solution STATE native contracts passed (11 cases)')


if __name__=='__main__':
    main()
