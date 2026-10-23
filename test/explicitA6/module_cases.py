from pathlib import Path
import shutil,json,subprocess,os,sys,re
base=Path(__file__).resolve().parent;root=base.parents[1];label=sys.argv[1]
binary=base/'baseline/sonicSolver' if label=='before' else base/'candidate/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
def write(p,object,identity,**body):p.write_text(json.dumps({'SonicFile':{'object':object,'type':identity},**body},indent=2))
records=[]
for name in ['absent','laminar','gravity','kepsilon','kepsilonGravity','disabled','override','sstOnce']:
 p=base/label/name
 if not p.exists():
  source=root/'test/resourceA5/after/stageSST' if label=='before' else base/'before'/name
  shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 if label=='before':
  time='classicalRK4' if name in ['override','sstOnce'] else 'forwardEuler'
  write(p/'solvers/numerics.yaml','solver','numerics',time={'default':time},terms={'convection':'weno5Rusanov','diffusion':'central2Explicit'})
  write(p/'solvers/runtime.yaml','solver','runtime',startTime=0,endTime=.001,maxDeltaT=.0001,CFL=.5,writeControl='timeStep',writeInterval=1,writeInitial=True,parallel={'enabled':False},output={'jobName':'C5','outputDir':'result'})
  params={'simulationType':'RAS','RAS':{'model':'kOmegaSST' if name=='sstOnce' else 'kEpsilon'},'transportMethod':'legacyOnceStep'}
  if name=='disabled':params['enabled']=False
  if name in ['absent','laminar','gravity']:params={'simulationType':'laminar'}
  write(p/'models/turbulence.yaml','models','turbulence',**params)
  write(p/'models/gravity.yaml','models','gravity',enabled=name in ['gravity','kepsilonGravity'],value=[.01,0,0])
  # Keep the baseline-native epsilon registration, but provide the actual active transport inputs.
  f=p/'fields/internalField.yaml';s=f.read_text();s=re.sub(r'(epsilon:\s*\n\s*(?:uniform|default):)[^\n]*',r'\1 0.001',s);f.write_text(s)
  f=p/'fields/fields.yaml';s=f.read_text()
  for variable in ['k','epsilon','omega']:
   active=name not in ['absent','laminar','gravity','disabled'] and (variable!='omega' if name!='sstOnce' else variable!='epsilon')
   s=re.sub(r'(^'+variable+r':\n(?:(?:  [^\n]*\n))*)',lambda m:re.sub(r'output: (true|false)','output: '+str(active).lower(),m.group(0)),s,flags=re.M)
  f.write_text(s)
  if name=='absent':
   f=p/'models/models.yaml';s=f.read_text();s=re.sub(r'^turbulence:\n(?:  [^\n]*\n)*','',s,flags=re.M);f.write_text(s)
  if name=='override':
   # An explicit supported WHICH overrides the module defaults, without a second HOW authority.
   f=p/'solvers/solvers.yaml';f.write_text((root/'test/resourceA5/after/stageSST/solvers/solvers.yaml').read_text()+'\nwhich:\n  type: providerBindings\n  file: solvers/which.yaml\n')
   write(p/'solvers/which.yaml','solver','which',bindings=[{'equation':v,'method':'TurbulenceTransport','inputs':['k','epsilon']} for v in ['k','epsilon']])
 binarySha=__import__('hashlib').sha256(binary.read_bytes()).hexdigest()
 for mode in ['explain','run']:
  cmd=[str(binary),mode,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(mode+'.log')).write_bytes(r.stdout)
  print(label,name,mode,r.returncode,flush=True)
  if r.returncode:raise RuntimeError(r.stdout[-1800:].decode())
 records.append({'case':name,'commands':[['explain',str(p)],['run',str(p)]],'binary':str(binary),'sha256':binarySha})
(base/'evidence'/('cases-'+label+'.json')).write_text(json.dumps(records,indent=2))
