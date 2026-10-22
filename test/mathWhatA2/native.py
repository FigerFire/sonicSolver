from pathlib import Path
import json,shutil,os,subprocess,hashlib,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1];binary=root/'build/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
case=base/'nativeReaction'
if not case.exists():shutil.copytree(base/'after/pairRK4',case,ignore=shutil.ignore_patterns('result','*.log','._*'))
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2))
p=case/'solvers/runtime.yaml';v=read(root/'test/mathWhatA1/after/pairRK4/solvers/runtime.yaml');v.update(endTime=.5,maxDeltaT=.05);write(p,v)
p=case/'fields/internalField.yaml';v=read(p);v['dye']={'default':1};v['salt']={'default':1.3};write(p,v)
p=case/'fields/boundaries.yaml';v=read(p)
for target in ['dye','salt']:v[target]={side:{'type':'zeroGradient'} for side in ['Left','Right','Bottom','Top','Front','Back']}
write(p,v)
p=case/'solvers/equationRegistry.yaml';v=read(root/'test/mathWhatA1/after/pairRK4/solvers/equationRegistry.yaml')
for equation,target,k in [('spreadDye','dye',.7),('spreadSalt','salt',1.1)]:
 v['add'][equation]['rhs']={'op':'subtract','args':[v['add'][equation]['rhs'],{'op':'multiply','args':[{'constant':k,'dimensions':[0,0,-1,0,0,0,0]},{'op':'pow','args':[{'symbol':target},{'constant':2}]}]}]}
write(p,v)
records=[]
def execute(p,action,success=True,reason=''):
 command=[str(binary),action,str(p)];r=subprocess.run(command,cwd=p,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_text(r.stdout)
 assert (r.returncode==0)==success,(command,r.stdout[-1500:])
 if reason:assert reason in r.stdout,r.stdout[-1500:]
 records.append(dict(case=str(p),command=command,exit=r.returncode,expected=reason))
for action in ['check','explain','run']:execute(case,action)
pvd=next((case/'result').glob('*.pvd'));frames=E.parse(pvd).findall('.//DataSet');assert len(frames)==11
field=E.parse(case/'result'/frames[-1].get('file'));stats={}
for target,c0,k in [('dye',1.,.7),('salt',1.3,1.1)]:
 values=list(map(float,next(a for a in field.findall('.//PointData/DataArray') if a.get('Name')==target).text.split()));exact=c0/(1+k*c0*.5);error=max(abs(x-exact) for x in values)
 assert error<2e-7,(target,error)
 stats[target]=dict(initial=c0,k=k,exact=exact,min=min(values),max=max(values),Linf=error)
for name,reason in [('implicit','Unsupported: LinearEquation'),('wrongUnits','Incompatible mathematical signatures')]:
 p=base/('native-'+name)
 if not p.exists():shutil.copytree(case,p,ignore=shutil.ignore_patterns('result','*.log','._*'))
 if name=='implicit':
  f=p/'solvers/which.yaml';v=read(f)
  for binding in v['bindings']:binding['method']='LinearEquation'
 else:
  f=p/'solvers/equationRegistry.yaml';v=read(f);v['add']['spreadDye']['rhs']['args'][1]['args'][0]['dimensions']=[0]*7
 write(f,v);execute(p,'explain',False,reason);assert not list((p/'result').glob('*.vts'))
(base/'evidence/native-reaction.json').write_text(json.dumps(dict(binary=hashlib.sha256(binary.read_bytes()).hexdigest(),runs=records,frames=len(frames),finalTime=frames[-1].get('timestep'),stats=stats),indent=2));print(json.dumps(stats,indent=2));print('native reaction and two negative bindings passed')
