from pathlib import Path
import json,subprocess,os
root=Path(__file__).resolve().parents[2];base=root/'test/nativeArchitectureClosure';out=base/'ctest-run';out.mkdir(exist_ok=True)
(base/'evidence').mkdir(exist_ok=True)
raw=subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],cwd=root,text=True)
(base/'evidence/ctest-inventory.json').write_text(raw)
inventory=json.loads(raw);selected=[];excluded=[]
def quote(s):return '"'+str(s).replace('\\','\\\\').replace('"','\\"')+'"'
lines=[]
for item in inventory['tests']:
 command=item['command'];name=item['name'];launcher=Path(command[0]).name
 if 'python' in launcher or 'sonicSolver'==launcher:
  excluded.append({'test':name,'reason':'Existing CFD/static-script test excluded from this native CTest subset; no historical long-time/temp-case script run. Relevant short CFD copies run separately.'});continue
 directory=base/'evidence'/'unit'/name;directory.mkdir(parents=True,exist_ok=True)
 if name in ('stageSSTTransport','sharedFluxTransport'):directory=base/'evidence'
 lines.append('add_test('+quote(name)+' '+' '.join(quote(x) for x in command)+')')
 properties=['WORKING_DIRECTORY',quote(directory),'TIMEOUT','180']
 # Preserve expected launch/environment behavior but relocate any temporary evidence.
 for prop in item['properties']:
  if prop['name'] in ('WILL_FAIL','PASS_REGULAR_EXPRESSION','FAIL_REGULAR_EXPRESSION'):
   value=prop['value'];properties += [prop['name'],quote(';'.join(value) if isinstance(value,list) else value)]
 lines.append('set_tests_properties('+quote(name)+' PROPERTIES '+' '.join(properties)+')');selected.append({'name':name,'command':command,'directory':str(directory),'MPI':'mpirun' in launcher or 'mpiexec' in launcher})
(out/'CTestTestfile.cmake').write_text('\n'.join(lines)+'\n');(base/'evidence/short-test-manifest.json').write_text(json.dumps({'selected':selected,'excluded':excluded},indent=2))
command=['ctest','--test-dir',str(out),'--output-on-failure','--parallel','1','-E','^(distributed|haloFailureContract|peskinPartialRows|sparseCanonicalCopy)']
with (base/'evidence/short-unit-ctest.log').open('w') as log:r=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')})
print('native subset CTest exit',r.returncode);raise SystemExit(r.returncode)
