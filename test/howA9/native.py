"""A9 native HOW/STATE execution, exact small-system references; all runs in test/."""
from pathlib import Path
import ast,json,shutil,subprocess,os,math,re,hashlib,sys,xml.etree.ElementTree as E
ROOT=Path(__file__).resolve().parents[2]
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'native'
BASE.mkdir(exist_ok=True);(BASE/'evidence').mkdir(exist_ok=True)
source=ast.parse((ROOT/'test/equationBinding/cases.py').read_text())
header=[]
for node in source.body:
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='records' for t in node.targets):break
 header.append(node)
# Use only the independent native input generator, never its executable or result references.
exec(compile(ast.Module(body=header,type_ignores=[]),'native-fixtures','exec'))
base=BASE;root=ROOT;binary=ROOT/'build/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(BASE/'evidence'),'SF_PLAN_TRACE':'1'}
def read(p):return json.loads(p.read_text())
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def edit(p,name,fn):q=p/'solvers'/name;v=read(q);fn(v);save(q,v)
def execute(p,mode='run',expected=None):
 cmd=[str(binary),mode,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(mode+'.log')).write_bytes(r.stdout)
 if expected is None:assert r.returncode==0,r.stdout[-5000:].decode()
 else:assert r.returncode!=0 and expected in r.stdout.decode(),r.stdout[-5000:].decode()
 return {'command':cmd,'code':r.returncode,'sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
def values(p):
 frames=E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=frames[-1];tree=E.parse(p/'result'/last.get('file'))
 return {v.get('Name'):list(map(float,v.text.split())) for v in tree.findall('.//PointData/DataArray')},float(last.get('timestep'))
def calls(p,items,equation=False):edit(p,'state.yaml',lambda v:v.update(calls=[{'equation' if equation else 'call':key,'reads':[{'symbol':name,'version':version} for name,version in inputs]} for key,inputs in items]))
def call(eq,target,id):return {'kind':'EquationCall','equation':eq,'target':target,'occurrence':id,'id':id}
def commit(id='publish'):return {'kind':'Commit','id':id,'sources':['x*','y*'],'targets':['x','y']}
def res(id,occurrence,kind='DiscreteEquation',tol=1e-11):return {'id':id,'occurrence':occurrence,'kind':kind,'norm':'Maximum','normalization':'None','absoluteTolerance':tol,'relativeTolerance':0}
def iter_case(name,jacobi=False,nested=False,stateChange=False):
 p,_=make(name,loop=True);how=read(p/'solvers/how.yaml');loop=how['root']['children'][0]
 loop.update(id='solve',repetitions=100,minimumIterations=1,residuals=['rx','ry'],children=[call('X','x*','xSolve'),call('Y','y*','ySolve'),commit()])
 if nested:loop.update(minimumIterations=2);how['root']['children']=[{'kind':'Loop','id':'outer','repetitions':3,'children':[loop]}]
 how['residuals']=[res('rx','xSolve','StateChange' if stateChange else 'DiscreteEquation'),res('ry','ySolve','StateChange' if stateChange else 'DiscreteEquation')]
 save(p/'solvers/how.yaml',how);calls(p,[('xSolve',[('y','Current')]),('ySolve',[('x' if jacobi else 'x*','Current')])])
 edit(p,'runtime.yaml',lambda v:v.update(startTime=3.,endTime=4.))
 return p
records=[]
for name,kw in [('gaussSeidel',{}),('previousIteration',{'jacobi':True}),('stateChange',{'stateChange':True}),('nested',{'nested':True})]:
 p=iter_case(name,**kw);execute(p,'explain');record=execute(p);data,t=values(p);err=max(abs(x-12/7) for x in data['x'])+max(abs(y-20/7) for y in data['y']);assert err<4e-11 and t==3,(name,err,t)
 lines=[s for s in (p/'run.log').read_text().splitlines() if s.startswith('[SF RES]')];rx=[float(re.search(r'norm=([^ ]+)',s).group(1)) for s in lines if 'metric=rx ' in s]
 # Independently execute the defined local linear recurrence, not the production AST.
 x=y=0.;reference=[]
 for n in range(100):
  xx=1+.25*y;yy=2+.5*(x if kw.get('jacobi') else xx)
  delta=max(abs(xx-x),abs(yy-y));x,y=xx,yy
  residual=max(abs(x-.25*y-1),abs(y-.5*x-2))
  reference.append((x,y))
  if (delta if kw.get('stateChange') else residual)<=1e-11:break
 expectedIterations=len(reference)+(4 if kw.get('nested') else 0)
 assert len(rx)==expectedIterations,(name,len(rx),expectedIterations)
 record.update(case=name,error=err,time=t,iterations=len(rx),firstResidual=rx[0]);records.append(record);print(name,err,len(rx),flush=True)
# Occurrence-specific inputs for repeated use of exactly the same WHAT equation.
p=iter_case('repeatOccurrence');edit(p,'how.yaml',lambda v:v.update(residuals=[],root={'kind':'Sequence','id':'HOW','children':[call('X','x*','predictX'),call('Y','y*','ySolve'),call('X','x*','correctX'),commit()]}))
calls(p,[('predictX',[('y','Current')]),('ySolve',[('x*','Current')]),('correctX',[('y*','Current')])]);record=execute(p);data,t=values(p);assert data['x']==[1.625]*45 and data['y']==[2.5]*45 and t==3;record.update(case=p.name,time=t);records.append(record)
# The correction is now an independently registered mathematical symbol.
p=BASE/'cases/correctionInput'
if p.exists():shutil.rmtree(p)
shutil.copytree(BASE/'cases/repeatOccurrence',p,ignore=shutil.ignore_patterns('result','*.log','._*'))
q=p/'fields/fields.yaml';v=read(q);v['deltaX']={'new':True,'type':'scalar','output':False,'dimensions':[0]*7};save(q,v)
q=p/'fields/internalField.yaml';v=read(q);v['deltaX']={'default':0.};save(q,v)
edit(p,'state.yaml',lambda v:v['use'].append('deltaX'))
edit(p,'equationRegistry.yaml',lambda v:(v['use'].append('PredictX'),v['add'].update(PredictX={'lhs':op('subtract',sym('deltaX'),op('multiply',num(.25),sym('y'))),'rhs':num(1)}),v['add']['Y'].update(lhs=op('subtract',sym('y'),op('multiply',num(.5),sym('deltaX'))))))
edit(p,'which.yaml',lambda v:v['bindings'].append({'equation':'PredictX','method':'LocalAffine'}))
edit(p,'how.yaml',lambda v:v['root']['children'][0].update(equation='PredictX',target='deltaX*'))
calls(p,[('predictX',[('y','Current')]),('ySolve',[('deltaX*','Current')]),('correctX',[('y*','Current')])]);record=execute(p);data,t=values(p);assert data['x']==[1.625]*45 and data['y']==[2.5]*45;record.update(case=p.name,time=t);records.append(record)
# Continuing after an iteration limit is explicit, observable, and does not claim convergence.
p=iter_case('continueAtLimit');edit(p,'how.yaml',lambda v:v['root']['children'][0].update(repetitions=1,failOnMaximum=False));record=execute(p);data,t=values(p);assert data['x']==[1.]*45 and data['y']==[2.5]*45 and t==3
assert 'maximum-continue-explicit' in (p/'run.log').read_text();record.update(case=p.name,time=t);records.append(record)
# Holding a coefficient between explicit writes is distinct from recomputing it inside an iteration.
for name,refresh in [('heldCoefficient',False),('updatedCoefficient',True)]:
 p,_=make(name);edit(p,'state.yaml',lambda v:v.update(use=['D','b']))
 q=p/'fields/fields.yaml';v=read(q)
 for key in ['D','D0']:v[key]['dimensions']=[0]*7
 save(q,v);q=p/'fields/internalField.yaml';v=read(q);v['D0']={'default':.2};v['b']={'default':0.};save(q,v)
 edit(p,'equationRegistry.yaml',lambda v:v.update(use=['C','B'],add={'C':{'lhs':sym('D'),'rhs':op('multiply',sym('D0'),op('add',num(1),sym('b')))},'B':{'lhs':sym('b'),'rhs':op('multiply',num(.5),sym('D'))}}))
 edit(p,'which.yaml',lambda v:v.update(bindings=[{'equation':e,'method':'LocalAffine'} for e in ['C','B']]))
 c=call('C','D*','coefficient');b=call('B','b','innerB');inner={'kind':'Loop','id':'inner','repetitions':4,'children':([c,b] if refresh else [b])}
 edit(p,'how.yaml',lambda v:v.update(root={'kind':'Sequence','id':'HOW','children':[{'kind':'Loop','id':'outer','repetitions':3,'children':([inner] if refresh else [c,inner])}]}))
 calls(p,[('coefficient',[('D0','Current'),('b','Current')]),('innerB',[('D*','Current')])]);record=execute(p);data,t=values(p);b=0.
 for outer in range(3):
  d=.2*(1+b)
  for inner in range(4):
   if refresh:d=.2*(1+b)
   b=.5*d
 assert max(abs(v-b) for v in data['b'])<2e-15 and t==0
 record.update(case=name,b=b,time=t);records.append(record);print(name,b,flush=True)
assert records[-1]['b']-records[-2]['b']>1e-4
# One real recipe stage nesting: bounded algebraic reevaluation, no repeated temporal RHS.
fields=[]
for name,nested,reverse in [('synchronous',False,False),('synchronousReverse',False,True),('stageInnerLoop',True,False),('stageInnerConvergence',True,False)]:
 p,_=make(name,stage=True)
 if nested:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(children=[{'kind':'Loop','id':'inner','repetitions':2,'children':['D -> D']},'G -> a','B -> b']))
 if name=='stageInnerConvergence':edit(p,'how.yaml',lambda v:(v.update(residuals=[res('stageD','D','StateChange')]),v['root']['children'][0]['children'][0].update(repetitions=10,minimumIterations=2,residuals=['stageD'])))
 if reverse:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(children=['D -> D','B -> b','G -> a']))
 execute(p,'explain');record=execute(p);data,t=values(p);assert abs(t-.2)<1e-15;fields.append(data);record.update(case=name,time=t);records.append(record);print(name,t,flush=True)
assert fields[0]==fields[1]==fields[2]==fields[3],'Stage order, LocalAffine terminal relation or publication differs'
# A satisfied initial relation must converge in one pass at every stage and step.
p,_=make('stageInitialGuess',stage=True)
edit(p,'equationRegistry.yaml',lambda v:v['add']['D'].update(rhs=sym('D0')))
q=p/'fields/internalField.yaml';v=read(q);v['D']={'default':.001};save(q,v)
calls(p,[('D',[('D0','Current')]),('G',[('D','Stage')])],equation=True)
edit(p,'how.yaml',lambda v:v.update(residuals=[res('stageD','D','StateChange')],root={'kind':'Sequence','id':'HOW','children':[{'kind':'StageLoop','id':'stages','children':[{'kind':'Loop','id':'initialGuess','repetitions':1,'minimumIterations':1,'residuals':['stageD'],'children':['D -> D']},'G -> a','B -> b']},{'kind':'Commit','id':'commit'}]}))
execute(p,'explain');record=execute(p);data,t=values(p)
lines=[s for s in (p/'run.log').read_text().splitlines() if s.startswith('[SF RES]')]
assert len(lines)==32 and all(float(re.search(r'norm=([^ ]+)',s).group(1))==0 for s in lines) and abs(t-.2)<1e-15,(lines,t)
lam=4*math.sin(math.pi/16)**2*64+1600;z=-lam*.001*.025
amplification=sum(z**k/math.factorial(k) for k in range(5))**8
error=max(abs(data['a'][(1*3+1)*9+i]-math.sin(math.pi*i/8)*amplification) for i in range(9))
assert error<1e-14,error
record.update(case=p.name,time=t,stageChecks=len(lines),error=error);records.append(record)
(BASE/'evidence/positive.json').write_text(json.dumps(records,indent=2))
negative=[]
def reject(name,source,edit,expect,mode='run'):
 p=BASE/'negative'/name
 if p.exists():shutil.rmtree(p,onerror=lambda f,q,e:None if isinstance(e[1],FileNotFoundError) else (_ for _ in ()).throw(e[1]))
 shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','._*'));edit(p);record=execute(p,mode,expect);record.update(case=name,expected=expect);negative.append(record);print('rejected',name,flush=True)
gs=BASE/'cases/gaussSeidel'
reject('whichReads',gs,lambda p:edit(p,'which.yaml',lambda v:v['bindings'][0].update(reads=[])),'retired')
reject('missingWorkWriter',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0]['children'].reverse()),'publication')
reject('missingResidual',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(residuals=['absent'])),'Missing/duplicate Loop residual')
reject('unimplementedMetric',gs,lambda p:edit(p,'how.yaml',lambda v:v['residuals'][0].update(kind='LinearSolver')),'cannot compute requested residual')
reject('noConvergence',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(repetitions=1)),'maximum iterations')
reject('badInputVersion',gs,lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][0].update(version='Frozen')),'input version/snapshot')
reject('unknownInput',gs,lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][0].update(symbol='undefined')),'lacks registered field metadata')
reject('targetIsInput',gs,lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'].append({'symbol':'x\'','version':'Current'})),'Legacy prime-as-Correction')
reject('duplicateState',gs,lambda p:edit(p,'state.yaml',lambda v:v['calls'].append(v['calls'][0])),'Duplicate STATE')
reject('ambiguousEquation',BASE/'cases/repeatOccurrence',lambda p:edit(p,'state.yaml',lambda v:(v['calls'][0].pop('call'),v['calls'][0].update(equation='X'))),'missing or ambiguous')
reject('invalidCommit',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0]['children'][-1].update(sources=['x','y*'])),'Working/Correction source')
reject('unknownEquation',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0]['children'][0].update(equation='unknown')),'not registered')
reject('resOutsideLoop',gs,lambda p:edit(p,'how.yaml',lambda v:v['root'].update(children=[v['root']['children'][0]['children'][0],{'kind':'Loop','id':'bad','repetitions':2,'residuals':['rx','ry'],'children':v['root']['children'][0]['children'][1:]}])),'outside its Loop')
reject('badStageRepeat',BASE/'cases/synchronous',lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(children=['D -> D',{'kind':'Loop','id':'wrongRhs','repetitions':2,'children':['G -> a']},'B -> b'])),'temporal RHS cannot repeat')
reject('doubleTimeCommit',BASE/'cases/synchronous',lambda p:edit(p,'how.yaml',lambda v:v['root']['children'].append({'kind':'Commit','id':'duplicateClock'})),'one terminal physical Commit')
reject('unsupportedSubcycle',gs,lambda p:edit(p,'how.yaml',lambda v:v['root']['children'][0].update(kind='Subcycle')),'Unknown HOW node kind')
reject('oldTimeWithoutSnapshot',gs,lambda p:edit(p,'state.yaml',lambda v:v['calls'][0]['reads'][0].update(version='OldTime')),'input version/snapshot')
reject('commitOutsideWriter',gs,lambda p:edit(p,'how.yaml',lambda v:v['root'].update(children=[commit('tooEarly'),v['root']['children'][0]])),'Variable Commit requires accepted')
reject('missingInitial',gs,lambda p:save(p/'fields/internalField.yaml',{k:v for k,v in read(p/'fields/internalField.yaml').items() if k!='y'}),'initialization')
reject('wrongDimensions',gs,lambda p:save(p/'fields/fields.yaml',{**read(p/'fields/fields.yaml'),'y':{'new':True,'type':'scalar','output':True,'dimensions':[0,1,0,0,0,0,0]}}),'Incompatible mathematical signatures')
(BASE/'evidence/negative.json').write_text(json.dumps(negative,indent=2))
print('A9 native positives and diagnostics passed')
