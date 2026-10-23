from pathlib import Path
import json,shutil,subprocess,os,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')};binary=root/'build/sonicSolver';records=[]
for name,action,reason in [('noWhich','explain','nonempty'),('missingField','explain','unknown field dye'),('missingInitial','explain','initialization'),('missingBoundary','run','Missing scalar boundary backing'),('unsupportedImplicit','explain','Unsupported'),('unknownMethod','explain','registered implementation')]:
 p=base/'negative'/name
 if not p.exists():shutil.copytree(base/'composition/flowPair',p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 rel={'noWhich':'solvers/which.yaml','missingField':'fields/fields.yaml','missingInitial':'fields/internalField.yaml','missingBoundary':'fields/boundaries.yaml','unsupportedImplicit':'solvers/which.yaml','unknownMethod':'solvers/which.yaml'}[name]
 f=p/rel;d=json.loads((base/'composition/flowPair'/rel).read_text())
 if name=='noWhich':d['bindings']=[]
 if name=='missingField':d.pop('dye')
 if name=='missingInitial':d.pop('dye')
 if name=='missingBoundary':d['dye'].pop('Left')
 if name=='unsupportedImplicit':d['bindings'][-1]['method']='LinearEquation'
 if name=='unknownMethod':d['bindings'][-1]['method']='UnimplementedA5Method'
 f.write_text(json.dumps(d,indent=2));cmd=[str(binary),action,str(p)]
 r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(p/(action+'.log')).write_text(r.stdout)
 assert r.returncode and reason in r.stdout,(name,r.stdout[-1600:])
 assert not list((p/'result').glob('*.vts'))
 records.append(dict(case=name,exit=r.returncode,command=cmd,expected=reason,reason=r.stdout.split('Reason:')[-1].strip()[-600:]));print(name,'rejected',flush=True)
(base/'evidence/negative.json').write_text(json.dumps(dict(binary=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=records),indent=2))
