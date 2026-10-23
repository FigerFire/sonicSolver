from pathlib import Path
import json,shutil,subprocess,os
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
for name,source,flow in [('referenceScalar','scalarPair',False),('referenceFlow','flowPair',True)]:
 p=base/'composition'/name
 if not p.exists():shutil.copytree(base/'composition'/source,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 if flow:
  def load(f):return json.loads(f.read_text())
  def save(f,d):f.write_text(json.dumps(d,indent=2))
  f=p/'solvers/state.yaml';d=load(f);d['use']=['rho','rhoU','rhoE'];save(f,d)
  f=p/'solvers/equationRegistry.yaml';d=load(f);d['use']=['Continuity','Momentum','Energy'];d.pop('add');save(f,d)
  f=p/'solvers/how.yaml';d=load(f);d['root']['children']=[x for x in d['root']['children'] if x.get('equation') not in ['spreadDye','spreadSalt']];save(f,d)
  f=p/'solvers/which.yaml';d=load(f);d['bindings']=[x for x in d['bindings'] if x['equation'] in ['continuity','momentum','energy']];save(f,d)
  for rel in ['fields/fields.yaml','fields/internalField.yaml','fields/boundaries.yaml']:
   f=p/rel;d=load(f);d.pop('dye',None);d.pop('salt',None);save(f,d)
 for action in ['explain','run']:
  cmd=[str(base/'baseline/sonicSolver'),action,str(p)]
  r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_bytes(r.stdout)
  print(name,action,r.returncode,flush=True)
  if r.returncode:raise RuntimeError(r.stdout[-2000:].decode())
