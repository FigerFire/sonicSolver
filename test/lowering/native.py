"""Lowering production contracts: minimal HOW, explicit references and real OldTime reads."""
from pathlib import Path
import ast,json,subprocess,os,math,hashlib,shutil,sys,xml.etree.ElementTree as E
root=Path(__file__).resolve().parents[2];base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'native';base.mkdir(exist_ok=True);(base/'evidence').mkdir(exist_ok=True)
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
 assert r.returncode==0 if expected is None else r.returncode!=0 and expected in r.stdout.decode(),r.stdout[-6000:].decode()
 return dict(command=cmd,exit=r.returncode)
def values(p):
 frames=E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=frames[-1];tree=E.parse(p/'result'/last.get('file'))
 return {v.get('Name'):list(map(float,v.text.split())) for v in tree.findall('.//PointData/DataArray')},float(last.get('timestep'))
def setup(name,old=False,rename=False):
 p,n=make(name,stage=True);v=read(p/'case.yaml');v['formatVersion']=2;save(p/'case.yaml',v)
 ids=['z-coefficient','a-residual','m-partner'] if rename else ['coefficient','diffuse','partner']
 eqs=['relation','transport','decay'] if rename else ['D','G','B']
 edit(p,'how.yaml',lambda v:v.update(root=[{'id':id,'equation':eq,'target':target,'frequency':'Stage'} for id,eq,target in zip(ids,eqs,['D','a','b'])]+[{'kind':'Commit','id':'accepted'}]))
 if rename:
  edit(p,'equationRegistry.yaml',lambda v:v.update(use=eqs,add={new:v['add'][old] for old,new in zip(['D','G','B'],eqs)}))
 edit(p,'equationRegistry.yaml',lambda v:v['add'][eqs[1]]['rhs'].update(occurrence='mainDiffusion'))
 edit(p,'state.yaml',lambda v:v.update(calls=[{'equation':eqs[0],'reads':[{'symbol':'D0','version':'Current'},{'symbol':"b'" if old else 'b','version':'OldTime' if old else 'Stage'}]},{'call':ids[1],'reads':[{'symbol':'D','version':'Stage'}]}]))
 # A Step override changes the numerical site without repeating an equation method.
 edit(p,'which.yaml',lambda v:v.update(bindings=[{'equation':eqs[0],'method':'LocalAffine'},{'equation':eqs[1],'method':'ExplicitCentral2Arithmetic','operators':[{'occurrence':'diffusion','provider':'scalar.central2.arithmetic-faces'}]},{'step':ids[1],'operators':[{'occurrence':'mainDiffusion','provider':'scalar.central2.arithmetic-faces'}],'time':'classicalRK4'},{'equation':eqs[2],'method':'ScalarDiffusionCentral2'}]))
 return p,n
records=[]
def polynomial(z):return 1+z+z*z/2+z**3/6+z**4/24
for name,old,renamed in [('stage',False,False),('renamed',False,True),('overlappingIds',False,False),('oldTime',True,False),('loadedClock',True,False)]:
 p,n=setup(name,old,renamed)
 if name=='overlappingIds':
  edit(p,'how.yaml',lambda v:[item.update(id=eq) for item,eq in zip(v['root'],['D','G','B'])])
  edit(p,'state.yaml',lambda v:v['calls'][1].update(call='G'))
  edit(p,'which.yaml',lambda v:v['bindings'][2].update(step='G'))
 if name=='loadedClock':edit(p,'runtime.yaml',lambda v:v.update(startTime=3.,endTime=3.2))
 execute(p,'explain');record=execute(p);data,t=values(p);assert abs(t-(3.2 if name=='loadedClock' else .2))<1e-14
 if old:
  q=polynomial(-.7*.025);lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2+1600
  amp=math.prod(polynomial(-lam*.001*(1+q**k)*.025) for k in range(8));err=max(abs(data['a'][(1*3+1)*n+i]-amp*math.sin(math.pi*i/(n-1))) for i in range(n));assert err<2e-14,err
  relation=max(abs(d-.001*(1+q**7)) for d in data['D']);assert relation<2e-17,relation
 else:err=0;relation=max(abs(d-.001*(1+b)) for d,b in zip(data['D'],data['b']));assert relation<2e-17
 assert 'recipe.stages' in (p/'explain.log').read_text() and 'mainDiffusion' in (p/'explain.log').read_text()
 record.update(case=name,error=err,relation=relation,time=t,oldTime=old);records.append(record);print(name,err,relation,flush=True)
left,_=values(base/'cases/stage');right,_=values(base/'cases/renamed');assert left==right==values(base/'cases/overlappingIds')[0],'Renaming or overlapping namespace IDs changed fields'
assert values(base/'cases/oldTime')[0]['a']!=left['a'],'OldTime was silently replaced by Stage'
negatives=[]
def negative(name,edit_case,expected,source='stage'):
 p=base/'negative'/name;p.parent.mkdir(exist_ok=True)
 if p.exists():shutil.rmtree(p,onerror=lambda fn,path,err:None if isinstance(err[1],FileNotFoundError) else (_ for _ in ()).throw(err[1]))
 shutil.copytree(base/'cases'/source,p,ignore=shutil.ignore_patterns('result','*.log','._*'));edit_case(p);r=execute(p,'explain',expected);r.update(case=name,diagnostic=expected);negatives.append(r);print('rejected',name,flush=True)
negative('legacyPrime',lambda p:(save(p/'case.yaml',{**read(p/'case.yaml'),'formatVersion':1}),edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][1].update(symbol="b'",version='OldTime'))),'Legacy prime-as-Correction')
negative('primeWrongVersion',lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][1].update(symbol="b'",version='Stage')),'OldTime')
negative('primeTarget',lambda p:edit(p,'how.yaml',lambda v:v['root'][0].update(target="D'")),'OldTime prime is read-only')
negative('authoredStages',lambda p:edit(p,'how.yaml',lambda v:v.update(root={'kind':'StageLoop','children':v['root'][:-1]})),'recipe')
negative('duplicateStep',lambda p:edit(p,'how.yaml',lambda v:v['root'][2].update(id='diffuse')),'Duplicate')
negative('missingStepReference',lambda p:edit(p,'state.yaml',lambda v:v['calls'][1].update(call='G')),'missing or ambiguous')
negative('timeConflict',lambda p:edit(p,'which.yaml',lambda v:v['bindings'][2].update(time='forwardEuler')),'Time recipe conflict')
negative('unknownOperator',lambda p:edit(p,'which.yaml',lambda v:v['bindings'][2]['operators'][0].update(provider='unimplemented')),'provider')
negative('ambiguousOperator',lambda p:(edit(p,'equationRegistry.yaml',lambda v:v['add']['G'].update(rhs=op('add',v['add']['G']['rhs'],op('diffusion',sym('D'),sym('a'))))),edit(p,'which.yaml',lambda v:v.update(bindings=[{'equation':'D','method':'LocalAffine'},{'equation':'G','method':'ExplicitCentral2Arithmetic','operators':[{'occurrence':'diffusion','provider':'scalar.central2.arithmetic-faces'}]},{'equation':'B','method':'ScalarDiffusionCentral2'}]))),'ambiguous operator')
negative('duplicateWhich',lambda p:edit(p,'which.yaml',lambda v:v['bindings'].append(v['bindings'][2])),'Conflicting')
negative('futureCoefficient',lambda p:edit(p,'how.yaml',lambda v:v.update(root=[v['root'][1],v['root'][0],*v['root'][2:]])),'preceding Stage algebraic publication')
negative('unsupportedHistory',lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][0].update(version='OldTime')),'old-time STATE')
roundtrip=base/'roundtrip'
r=subprocess.run([str(root/'build/test_lowering'),str(base/'cases/oldTime'),str(roundtrip)],env=env,cwd=base,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(base/'evidence/roundtrip.log').write_bytes(r.stdout);assert r.returncode==0,r.stdout.decode()
execute(roundtrip,'explain');execute(roundtrip);assert values(roundtrip)==values(base/'cases/oldTime'),'Read/write/read changed actual execution'
(base/'evidence/verification.json').write_text(json.dumps(dict(binary=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=records,negative=negatives),indent=2))
print('Minimal HOW/rename, partial WHICH, frozen OldTime, loaded clock and negative contracts passed')
