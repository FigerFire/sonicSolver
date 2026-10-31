"""Targeted current tests; all work/logs below lowering, frozen evidence untouched."""
from pathlib import Path
import json,subprocess,os,ast
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')}
inventorySource=ast.parse((root/'test/linearA8/run_tests.py').read_text())
selected=next(ast.literal_eval(n.value) for n in inventorySource.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='selected' for t in n.targets))|{'howStateContracts','nativeHowStateContracts','loweringContracts','nativeLowering','scopeResources','nativeScopes'}
inventory=json.loads(subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],env=env,text=True))
q=lambda v:json.dumps(str(v));lines=[];found=[]
for item in inventory['tests']:
 if item['name'] not in selected:continue
 name=item['name'];found.append(name);command=list(item['command'])
 if name=='nativeLowering':command.append(str(base/'nativeLowering'))
 if name=='nativeHowStateContracts':command.append(str(base/'nativeA9'))
 if name=='nativeEquationTargets':command.append(str(base/'targetsAfter'))
 if name=='nativeLinearBlock':command.append(str(base/'after/blocks'))
 if name=='nativeDiffusionGhostBoundary':command.append(str(base/'after/boundaries'))
 if name=='nonlinearScalarExecution':command[-1]=str(base/'after/nonlinear')
 cwd=base if name in {'explicitEquationCapabilities','coefficientStateExecution'} else base/'evidence/unit'/name;cwd.mkdir(parents=True,exist_ok=True)
 lines+=['add_test('+q(name)+' '+' '.join(q(c) for c in command)+')','set_tests_properties('+q(name)+' PROPERTIES WORKING_DIRECTORY '+q(cwd)+' TIMEOUT 180)']
assert set(found)==selected,selected-set(found)
folder=base/'ctest';folder.mkdir(exist_ok=True);(folder/'CTestTestfile.cmake').write_text('\n'.join(lines));(base/'evidence/test-selection.json').write_text(json.dumps(found,indent=2))
with (base/'evidence/ctest-after.log').open('w') as log:r=subprocess.run(['ctest','--test-dir',str(folder),'--output-on-failure','--parallel','1'],env=env,cwd=base,stdout=log,stderr=subprocess.STDOUT)
print((base/'evidence/ctest-after.log').read_text());raise SystemExit(r.returncode)
