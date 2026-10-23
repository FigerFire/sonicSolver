from pathlib import Path
import shutil,os,subprocess,sys,hashlib,json,re
base=Path(__file__).resolve().parent;root=base.parents[1];os.chdir(base)
label=sys.argv[1];binary=base/'baseline/sonicSolver' if label=='before' else root/'build/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
sources={n:root/'test/coreSinglePath/before'/n for n in ['PISO','SIMPLE','PIMPLE','conservativePressure','sod','stageSST','eulerianSST']}
sources.update({n:root/'test/nativeEquationClosure/cases'/n for n in ['pairEuler','pairRK4','renamedRK4','tripleRK4','linearTermsRK4']})
sources['levelSet']=root/'test/coreSinglePath/interface/after/final/preset'
records=[]
for name,source in sources.items():
 p=base/label/name
 if not p.exists():shutil.copytree(source if label=='before' else base/'before'/name,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store','numerical.txt'))
 if label=='before' and name in ['pairEuler','pairRK4','renamedRK4','tripleRK4','linearTermsRK4']:
  field=p/'fields/fields.yaml';d=json.loads(field.read_text())
  for v in d.values():
   if isinstance(v,dict) and v.get('new'):v['dimensions']=[0]*7
  field.write_text(json.dumps(d,indent=2))
  eq=p/'solvers/equationRegistry.yaml';d=json.loads(eq.read_text())
  def annotate(x):
   if not isinstance(x,dict):return
   if x.get('op')=='diffusion':x['args'][0]['dimensions']=[0,2,-1,0,0,0,0]
   if x.get('op')=='multiply' and 'constant' in x['args'][0]:x['args'][0]['dimensions']=[0,0,-1,0,0,0,0]
   if x.get('op')=='add' and 'constant' in x['args'][1]:x['args'][1]['dimensions']=[0,0,-1,0,0,0,0]
   for v in x.values():
    if isinstance(v,list):
     for y in v:annotate(y)
    else:annotate(v)
  annotate(d);eq.write_text(json.dumps(d,indent=2))
 e={**env,'SF_PRESSURE_TRACE_DIR':str(p/'result/assembly')} if name in ['PISO','SIMPLE','PIMPLE'] else env
 for action in ['explain','run']:
  r=subprocess.run([str(binary),action,str(p)],cwd=p,env=e,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  (p/(action+'.log')).write_bytes(r.stdout)
  if r.returncode:raise RuntimeError(name+' '+action+' '+r.stdout[-2200:].decode(errors='replace'))
 records.append({'case':name,'source':str(source),'binary':hashlib.sha256(binary.read_bytes()).hexdigest(),'path':str(p)})
 print(label,name,'completed',flush=True)
(base/'evidence'/('cases-'+label+'.json')).write_text(json.dumps(records,indent=2))
