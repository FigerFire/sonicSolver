"""Exercise the ordinary CLI; C++ injection is an explicitly separate numerical reference."""
from pathlib import Path
from collections import defaultdict,Counter
import subprocess,shutil,json,os,hashlib,re,math,xml.etree.ElementTree as ET
from create_cases import CASES
base=Path(__file__).resolve().parent;root=base.parents[1]
solver=root/'build/sonicSolver';reference=base/'bin/test_nativeEquationReference'
runs=base/'evidence/runs';runs.mkdir(parents=True,exist_ok=True)
run=runs/f'{len(list(runs.iterdir()))+1:03d}';run.mkdir()
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
records=[];negatives=[]
def invoke(args,cwd,log,success=True):
 r=subprocess.run([str(a) for a in args],cwd=cwd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 log.write_text(r.stdout)
 assert (r.returncode==0)==success,(args,r.returncode,r.stdout[-2500:])
 return r.stdout

def copy(source,p):
 shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','reference','*.log','._*'))
 # These test fixtures declare dimensionless tracers and SI diffusivity/rate constants.
 # Historical source inputs remain frozen; no production reader fallback is involved.
 path=p/'fields/fields.yaml';doc=json.loads(path.read_text())
 for value in doc.values():
  if isinstance(value,dict) and value.get('new'):value['dimensions']=[0]*7
 path.write_text(json.dumps(doc,indent=2))
 path=p/'solvers/equationRegistry.yaml';doc=json.loads(path.read_text())
 def declarations(node):
  if not isinstance(node,dict):return
  if node.get('op')=='diffusion':node['args'][0]['dimensions']=[0,2,-1,0,0,0,0]
  if node.get('op')=='multiply' and 'constant' in node['args'][0]:node['args'][0]['dimensions']=[0,0,-1,0,0,0,0]
  if node.get('op')=='add' and 'constant' in node['args'][1]:node['args'][1]['dimensions']=[0,0,-1,0,0,0,0]
  for value in node.values():
   if isinstance(value,list):
    for item in value:declarations(item)
   else:declarations(value)
 declarations(doc);path.write_text(json.dumps(doc,indent=2))
 return p

def jsonfile(p):return json.loads(p.read_text())
def edit(p,action):
 v=jsonfile(p);action(v);p.write_text(json.dumps(v,indent=2)+'\n')
def content(p):
 x=ET.parse(p);return {a.attrib['Name']:list(map(float,a.text.split())) for a in x.findall('.//PointData/DataArray')}
specs=[dict(name=n,recipe=r,ids=s,equations=e,diffusivities=[.03,.08,.11][:len(s)],path=str(base/'cases'/n)) for n,r,s,e in CASES]
for spec in specs:
 p=copy(Path(spec['path']),run/spec['name'])
 for action in ['check','explain','run']:invoke([solver,action,p],p,p/(action+'.log'))
 explain=(p/'explain.log').read_text();log=(p/'run.log').read_text()
 assert 'Q=0, EOS=none, scalar arrays='+str(len(spec['ids'])) in log
 assert 'flow.conservative' not in explain and 'thermodynamics.single-fluid' not in explain
 assert '/add/' in explain and '/solvers/equationRegistry.yaml' in explain
 ref=p/'reference';ref.mkdir()
 args=[reference,p,ref]
 for slot,(s,e,d) in enumerate(zip(spec['ids'],spec['equations'],spec['diffusivities'])):
  coefficient=str(d)+(':0.2:0' if slot==0 else ':0:0.15') if spec['name']=='linearTermsRK4' else str(d)
  args += [s,e,coefficient]
 reflog=invoke(args,ref,ref/'run.log')
 frames=sorted(f for f in (p/'result').glob('*.vts') if not f.name.startswith('._'));other=sorted(f for f in (ref/'result').glob('*.vts') if not f.name.startswith('._'))
 assert len(frames)==len(other)==11
 for a,b in zip(frames,other):
  aa,bb=content(a),content(b)
  for name in spec['ids']:assert aa[name]==bb[name],('First output difference',a.name,name)
  assert all(math.isfinite(v) for values in aa.values() for v in values)
 for name in spec['ids']:assert content(frames[0])[name]!=content(frames[-1])[name],('Scalar never updated',name)
 times=[float(x.attrib['timestep']) for x in ET.parse(next(f for f in (p/'result').glob('*.pvd') if not f.name.startswith('._'))).findall('.//DataSet')]
 assert len(times)==11 and times[0]==0 and times[-1]==.005
 plan=lambda text:[x for x in text.splitlines() if x.startswith('[SF PLAN]')]
 assert plan(log)==plan(reflog),'Native/C++ execution order or clock differs'
 stages=lambda text:[x for x in text.splitlines() if x.startswith('[SF SCALAR stage]')]
 assert stages(log)==stages(reflog),'Native/C++ actual stage time differs'
 count=1 if spec['recipe']=='forwardEuler' else 4
 abscissae=[0] if count==1 else [0,.5,.5,1]
 stageRows=re.findall(r'\[SF SCALAR stage\] target=(\S+) step=(\d+) stage=(\d+) time=(\S+) dt=(\S+)',log)
 assert len(stageRows)==10*count*len(spec['ids'])
 for name,step,stage,t,dt in stageRows:
  assert name in spec['ids'] and math.isclose(float(t),times[int(step)]+abscissae[int(stage)]*float(dt),rel_tol=0,abs_tol=2e-17)
 ops=defaultdict(list)
 for step,operation in re.findall(r'\[SF PLAN\] step=(\d+).*?op=(\S+)',log):ops[int(step)].append(operation)
 assert len(ops)==10
 for operations in ops.values():
  assert Counter(operations)['time.commit']==1
  for name in spec['ids']:
   assert sum(x.startswith('scalar.stage.advance@') and x.endswith('advance.'+name) for x in operations)==count
   assert sum(x.startswith('scalar.publish@') and x.endswith('advance.'+name) for x in operations)==1
 pointers=re.findall(r'\[SF SCALAR storage\] target=\S+ old=(\S+) stage=(\S+) rhs=(\S+)',log)
 assert len(pointers)==len(spec['ids']) and len({p for row in pointers for p in row})==3*len(spec['ids'])
 final=content(frames[-1]);stats={name:{'min':min(final[name]),'max':max(final[name]),'L1':math.fsum(abs(x) for x in final[name]),'L2':math.sqrt(math.fsum(x*x for x in final[name]))} for name in spec['ids']}
 records.append({**spec,'case':str(p),'commands':[['sonicSolver',a,str(p)] for a in ['check','explain','run']],'referenceCommand':list(map(str,args)),'steps':10,'finalTime':times[-1],'frames':11,'allScalarValuesExact':True,'planAndStageTimesExact':True,'stats':stats})
 print(spec['name'],'CLI and C++ reference exact',flush=True)

final=lambda name,symbol:content(sorted(f for f in (run/name/'result').glob('*.vts') if not f.name.startswith('._'))[-1])[symbol]
for old,new in [('dye','colourant'),('salt','dopant')]:
 assert final('pairRK4',old)==final('renamedRK4',new)==final('tripleRK4',new),'Name/instance count affected an independent scalar'
assert max(abs(a-b-.15*.005) for a,b in zip(final('linearTermsRK4','salt'),final('pairRK4','salt')))<1e-14,'Constant source was not integrated'
assert max(abs(a-b) for a,b in zip(final('linearTermsRK4','dye'),final('pairRK4','dye')))>1e-6,'Linear reaction never affected the solution'

source=Path(specs[0]['path'])
def negative(name,relative,action,expected,mode='check',extra=None):
 p=copy(source,run/('negative-'+name));edit(p/relative,action)
 if extra:edit(p/extra[0],extra[1])
 text=invoke([solver,mode,p],p,p/(mode+'.log'),False)
 assert expected in text,(name,expected,text[-2500:])
 assert not list((p/'result').glob('*.vts')),'Rejected binding wrote a physical output frame'
 negatives.append({'name':name,'case':str(p),'mode':mode,'changedFile':relative,'additionalChangedFile':extra[0] if extra else None,'expected':expected,'actual':text[-1800:]})
q='solvers/equationRegistry.yaml';h='solvers/how.yaml';s='solvers/state.yaml';w='solvers/which.yaml'
negative('unknown-symbol',q,lambda x:x['add']['spreadDye']['rhs']['args'][1].update(symbol='notRegistered'),'Unknown declared mathematical symbol')
negative('duplicate-equation',q,lambda x:x['use'].append('spreadDye'),'duplicate/empty WHAT')
negative('duplicate-state',s,lambda x:x['use'].append('dye'),'Duplicate solution STATE')
negative('missing-target',h,lambda x:x['root']['children'][0].pop('target'),'requires equation/target')
negative('wrong-target',h,lambda x:x['root']['children'][0].update(target='notRegistered'),'Unknown')
negative('no-provider',w,lambda x:x['bindings'].pop(0),'no provider for equation occurrence')
negative('unknown-provider',w,lambda x:x['bindings'][0].update(method='NotImplemented'),'has no registered implementation')
negative('unknown-operator',q,lambda x:x['add']['spreadDye']['rhs'].update(op='curl'),'Math operator curl')
negative('bad-arity',q,lambda x:x['add']['spreadDye']['rhs']['args'].pop(),'Wrong operand count')
negative('unsupported-local-operator',q,lambda x:x['add']['spreadDye'].update(rhs={'op':'add','args':[x['add']['spreadDye']['rhs'],{'op':'multiply','args':[{'constant':1,'dimensions':[0,0,-1,0,0,0,0]},{'op':'sin','args':[{'symbol':'dye'}]}]}]}),'has no compatible numerical provider')
negative('unsupported-recipe','solvers/numerics.yaml',lambda x:x['time'].update(default='SSPRK3'),'recipe')
negative('missing-initial','fields/internalField.yaml',lambda x:x.pop('dye'),'missing initialization backing','run')
negative('missing-bc','fields/boundaries.yaml',lambda x:x.pop('dye'),'Missing scalar boundary backing','run')
negative('wrong-bc','fields/boundaries.yaml',lambda x:x['dye']['Left'].update(type='symmetry'),'fixedValue/zeroGradient only','run')
negative('duplicate-commit',h,lambda x:x['root']['children'].append({'kind':'Commit','id':'another'}),'Commit')
negative('missing-commit',h,lambda x:x['root']['children'].pop(),'Commit')
negative('operation-collision',h,lambda x:x['root']['children'][1].update(occurrence='advance.dye'),'Duplicate')
negative('target-workspace-collision',h,lambda x:x['root']['children'][1].update(target='dye'),'positive additive ddt(target)')
negative('duplicate-physical-writer',h,lambda x:x['root']['children'][1].update(target='dye'),'Duplicate temporal physical target',extra=(q,lambda x:x['add'].update(spreadSalt={'lhs':{'op':'ddt','args':[{'symbol':'dye'}]},'rhs':{'op':'diffusion','args':[{'constant':.08,'dimensions':[0,2,-1,0,0,0,0]},{'symbol':'dye'}]}})))
negative('unregistered-initial','fields/fields.yaml',lambda x:x.pop('internalField'),'missing initialization backing','run')
negative('wrong-child-identity',s,lambda x:x['SonicFile'].update(type='stateRegistry'),'expected object solver type state')
negative('incomplete-fluid',q,lambda x:(x.update(use=['Momentum']),x.pop('add')),'requires Continuity and Momentum')
negative('fluid-no-eos',q,lambda x:(x.update(use=['Continuity','Momentum','Energy']),x.pop('add')),'equationOfState')
negative('unsupported-replace',q,lambda x:x.update(replace={}),'lowering custom modifications is not implemented')
negative('duplicate-which',w,lambda x:x['bindings'].append(x['bindings'][0]),'Conflicting highest-priority')
negative('unselected-target',s,lambda x:x['use'].pop(),'unknown declared STATE symbol')
negative('missing-state-registration','solvers/solvers.yaml',lambda x:x.pop('state'),'explicit solution STATE')
negative('scalar-mpi','solvers/runtime.yaml',lambda x:x.update(parallel={'enabled':True,'nProcs':2}),'serial')
negative('negative-diffusion',q,lambda x:x['add']['spreadDye']['rhs']['args'][0].update(constant=-.03),'diffusivity must be finite and positive')
negative('missing-time-recipe','solvers/numerics.yaml',lambda x:x.pop('time'),'explicit time recipe')
negative('repeated-physical-loop',h,lambda x:x['root'].update(children=[{'kind':'Loop','repetitions':2,'children':x['root']['children'][:-1]},x['root']['children'][-1]]),'authored HOW topology')
negative('numerics-cannot-replace-D',w,lambda x:x['bindings'][0].update(parameters={'D':.5}),'D belongs to WHAT')

# A pre-existing native parser cannot run these newly supported inputs. This is not a numerical before baseline.
old=base/'baseline/sonicSolver'
if old.exists():invoke([old,'check',run/specs[0]['name']],run,run/'baseline-native-rejection.log',False)
manifest={'solver':str(solver),'sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),'referenceSha256':hashlib.sha256(reference.read_bytes()).hexdigest(),'positive':records,'negative':negatives,'mpiNumericalRun':False}
(run/'manifest.json').write_text(json.dumps(manifest,indent=2));(base/'evidence/latest.json').write_text(json.dumps({'run':str(run),'manifest':str(run/'manifest.json')},indent=2))
print('Passed',len(records),'native/reference cases and',len(negatives),'negative inputs')
