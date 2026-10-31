"""Explicit Step/Stage scopes with independent discrete references and negative contracts."""
from pathlib import Path
import ast,json,subprocess,os,math,hashlib,shutil,sys,xml.etree.ElementTree as E
root=Path(__file__).resolve().parents[2];base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'native';base.mkdir(parents=True,exist_ok=True);(base/'evidence').mkdir(exist_ok=True)
source=ast.parse((root/'test/equationBinding/cases.py').read_text());header=[]
for node in source.body:
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='records' for t in node.targets):break
 header.append(node)
exec(compile(ast.Module(body=header,type_ignores=[]),'input-generator','exec'))
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'native';binary=root/'build/sonicSolver';env={**os.environ,'TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','DEVELOPER_DIR':'/Library/Developer/CommandLineTools'}
def read(p):return json.loads(p.read_text())
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def edit(p,name,fn):f=p/'solvers'/name;v=read(f);fn(v);save(f,v)
def execute(p,mode='run',expected=None):
 cmd=[str(binary),mode,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(mode+'.log')).write_bytes(r.stdout)
 assert (r.returncode==0) if expected is None else (r.returncode!=0 and expected in r.stdout.decode()),r.stdout[-5000:].decode()
 return dict(command=cmd,exit=r.returncode)
def values(p):
 frames=E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=frames[-1];tree=E.parse(p/'result'/last.get('file'))
 return {v.get('Name'):list(map(float,v.text.split())) for v in tree.findall('.//PointData/DataArray')},float(last.get('timestep'))
def setup(name,held=False,nested=False):
 p,n=make(name,stage=True);save(p/'case.yaml',{**read(p/'case.yaml'),'formatVersion':2})
 calls=[{'id':id,'equation':eq,'target':target} for id,eq,target in [('coefficient','D','D'),('diffuse','G','a'),('decay','B','b')]]
 if held:
  calls[0]['frequency']='Step';tree=[calls[0],{'kind':'Sequence','id':'transport','frequency':'Stage','children':calls[1:]}]
 else:
  if nested:calls[0]={'kind':'Loop','id':'outer','repetitions':2,'children':[{'kind':'Loop','id':'inner','repetitions':2,'children':[calls[0]]}]}
  tree=[{'kind':'Sequence','id':'synchronous','frequency':'Stage','children':calls}]
 edit(p,'how.yaml',lambda v:v.update(root=tree+[{'kind':'Commit','id':'accepted'}]))
 edit(p,'state.yaml',lambda v:v.update(calls=[{'equation':'D','reads':[{'symbol':'D0','version':'Current'},{'symbol':'b','version':'Current' if held else 'Stage'}]},{'equation':'G','reads':[{'symbol':'D','version':'Current' if held else 'Stage'}]}]))
 return p,n
R=lambda z:1+z+z*z/2+z**3/6+z**4/24
records=[]
for name,held,nested in [('held',True,False),('stage',False,False),('nested',False,True)]:
 p,n=setup(name,held,nested);execute(p,'explain');record=execute(p);data,t=values(p);assert abs(t-.2)<1e-14
 lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2+1600;q=R(-.7*.025)
 if held:
  amp=math.prod(R(-lam*.001*(1+q**k)*.025) for k in range(8));d=.001*(1+q**7)
 else:
  # Independent RK4 of two amplitudes; D is evaluated on the synchronous b stage.
  a,b=1.,1.;dt=.025
  f=lambda a,b:(-lam*.001*(1+b)*a,-.7*b)
  for _ in range(8):
   k1=f(a,b);k2=f(a+dt/2*k1[0],b+dt/2*k1[1]);k3=f(a+dt/2*k2[0],b+dt/2*k2[1]);k4=f(a+dt*k3[0],b+dt*k3[1]);a+=dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);b+=dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
  amp=a;d=.001*(1+b)
 error=max(abs(data['a'][4*n+i]-amp*math.sin(math.pi*i/(n-1))) for i in range(n));relation=max(abs(v-d) for v in data['D']);assert error<3e-14 and relation<3e-17,(name,error,relation)
 trace=(p/'run.log').read_text();explain=(p/'explain.log').read_text();assert 'frequency=Stage' in explain
 if held:assert 'frequency=Step' in explain
 if nested:assert 'Iteration within Stage' in explain
 operations=[line for line in trace.splitlines() if line.startswith('[SF PLAN]')]
 solves=[line for line in operations if 'op=affine.solve@' in line]
 assert len(solves)==(8 if held else 128 if nested else 32),(name,len(solves))
 assert all(('stage=' in line)!=held for line in solves)
 assert sum('op=time.commit' in line for line in operations)==8
 record.update(coefficientEvaluations=len(solves),physicalCommits=8)
 record.update(case=name,error=error,coefficientError=relation,finalTime=t);records.append(record);print(name,error,relation,flush=True)
assert values(base/'cases/stage')[0]==values(base/'cases/nested')[0]
assert values(base/'cases/held')[0]['a']!=values(base/'cases/stage')[0]['a'],'Different declared time semantics collapsed'
negative=[]
def reject(name,change,expected,source='stage'):
 p=base/'negative'/name
 if p.exists():shutil.rmtree(p,onerror=lambda fn,path,err:None if isinstance(err[1],FileNotFoundError) else (_ for _ in ()).throw(err[1]))
 p.parent.mkdir(exist_ok=True);shutil.copytree(base/'cases'/source,p,ignore=shutil.ignore_patterns('result','*.log','._*'));change(p);record=execute(p,'explain',expected);record.update(case=name,diagnostic=expected);negative.append(record);print('rejected',name,flush=True)
reject('ambiguous',lambda p:edit(p,'how.yaml',lambda v:v['root'][0].pop('frequency')),'Ambiguous HOW frequency')
reject('ddtStep',lambda p:edit(p,'how.yaml',lambda v:v.update(root=[{**v['root'][0]['children'][0],'frequency':'Step'},{**v['root'][0]['children'][1],'frequency':'Step'},{**v['root'][0]['children'][2],'frequency':'Stage'},v['root'][-1]])),'temporal HOW scope/frequency')
reject('escape',lambda p:edit(p,'how.yaml',lambda v:v['root'][0]['children'][0].update(frequency='Step')),'cannot escape')
reject('noncontiguous',lambda p:edit(p,'how.yaml',lambda v:v.update(root=[{**v['root'][0]['children'][1],'frequency':'Stage'},{**v['root'][0]['children'][0],'frequency':'Step'},{**v['root'][0]['children'][2],'frequency':'Stage'},{**v['root'][0]['children'][0],'id':'anotherCoefficient','frequency':'Stage'},v['root'][-1]])),'Non-contiguous Stage')
reject('legacyFrequency',lambda p:save(p/'case.yaml',{**read(p/'case.yaml'),'formatVersion':1}),'frequency requires formatVersion: 2')
reject('unknownFrequency',lambda p:edit(p,'how.yaml',lambda v:v['root'][0].update(frequency='automatic')),'frequency')
reject('commitFrequency',lambda p:edit(p,'how.yaml',lambda v:v['root'][-1].update(frequency='Step')),'frequency')
reject('stageWithoutTime',lambda p:edit(p,'how.yaml',lambda v:v['root'][0].update(children=v['root'][0]['children'][:1])),'Stage frequency requires an implemented temporal group')
reject('repeatedTime',lambda p:edit(p,'how.yaml',lambda v:v['root'][0].update(kind='Loop',repetitions=2)),'temporal RHS cannot repeat')
reject('stepReadsStage',lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][1].update(version='Stage')),'LocalAffine',source='held')
# Round-trip preserves explicit scope declarations through production IO.
roundtrip=base/'roundtrip';r=subprocess.run([str(root/'build/test_lowering'),str(base/'cases/held'),str(roundtrip)],cwd=base,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(base/'evidence/roundtrip.log').write_bytes(r.stdout);assert r.returncode==0,r.stdout.decode()
execute(roundtrip);assert values(roundtrip)==values(base/'cases/held')
(base/'evidence/verification.json').write_text(json.dumps(dict(binary=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=records,negative=negative),indent=2)+'\n')
print('Explicit frequency, independent stage/step references, nested iteration and negative contracts passed')
