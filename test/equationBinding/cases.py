"""Unified Equation/HOW native production tests; independent references, no external fixtures."""
from pathlib import Path
import json,math,os,subprocess,hashlib,xml.etree.ElementTree as E,shutil,sys
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent;root=Path(__file__).resolve().parents[2];binary=root/'build/sonicSolver';(base/'evidence').mkdir(parents=True,exist_ok=True)
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
def write(p,obj,kind,body):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps({'SonicFile':{'object':obj,'type':kind},**body},indent=2)+'\n')
def op(name,*args):return {'op':name,'args':list(args)}
def sym(name):return {'symbol':name}
def num(v,d=None):return {'constant':v,**({'dimensions':d} if d is not None else {})}
dim=[0,2,-1,0,0,0,0];rate=[0,0,-1,0,0,0,0]
def make(name,stage=False,residual=False,target='D',loop=False,dt=.025,renamed=False):
 p=base/'cases'/name;n=9 if stage else 5;fields=['a','b','D','D0'] if stage else ['x','y'] if loop else ['D','b','D0'];targets=['a','b','D'] if stage else ['x','y'] if loop else [target]
 eqid='sameMeaningDifferentId' if renamed else 'D';D=sym('D');rhs=op('multiply',sym('D0'),op('add',num(1),sym('b')))
 equations={eqid:{'lhs':op('subtract',D,rhs) if residual else D,'rhs':num(0) if residual else rhs}}
 if stage:equations.update(G={'lhs':op('ddt',sym('a')),'rhs':op('diffusion',D,sym('a'))},B={'lhs':op('ddt',sym('b')),'rhs':op('multiply',num(-.7,rate),sym('b'))})
 if loop:equations={'X':{'lhs':op('subtract',sym('x'),op('multiply',num(.25),sym('y'))),'rhs':num(1)},'Y':{'lhs':op('subtract',sym('y'),op('multiply',num(.5),sym('x'))),'rhs':num(2)}}
 bindings=[{'equation':eqid,'method':'LocalAffine'}]
 callInputs=[{'equation':eqid,'reads':[{'symbol':'D0','version':'Current'},{'symbol':'b' if target=='D' else 'D','version':'Stage' if stage else 'Current'}]}]
 how=[eqid+' -> '+target]
 if stage:
  bindings+=[{'equation':'G','method':'ExplicitCentral2Arithmetic'},{'equation':'B','method':'ScalarDiffusionCentral2'}]
  callInputs+=[{'equation':'G','reads':[{'symbol':'D','version':'Stage'}]}]
  how=[{'kind':'StageLoop','id':'stages','children':[eqid+' -> D','G -> a','B -> b']},{'kind':'Commit','id':'commit'}]
 if loop:callInputs=[];bindings=[{'equation':id,'method':'LocalAffine'} for id in ['X','Y']];how=[{'kind':'Loop','id':'fixedPoint','repetitions':30,'children':['X -> x','Y -> y']}]
 write(p/'case.yaml','case','registry',{'formatVersion':1,'name':name})
 roles=[('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]
 write(p/'solvers/solvers.yaml','solver','registry',{k:{'type':v,'file':'solvers/'+k+'.yaml'} for k,v in roles})
 write(p/'solvers/runtime.yaml','solver','runtime',{'startTime':0,'endTime':.2,'maxDeltaT':dt,'CFL':.5,'writeControl':'timeStep','writeInterval':1 if stage else 1000,'writeInitial':False,'parallel':{'enabled':False},'output':{'jobName':name,'outputDir':'result'}})
 write(p/'solvers/numerics.yaml','solver','numerics',{'time':{'default':'classicalRK4'}})
 write(p/'solvers/state.yaml','solver','state',{'use':targets,'calls':callInputs});write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',{'use':list(equations),'add':equations})
 write(p/'solvers/execution.yaml','solver','execution',{'Explicit':{}});write(p/'solvers/how.yaml','solver','how',{'root':{'kind':'Sequence','id':'HOW','children':how}});write(p/'solvers/which.yaml','solver','which',{'bindings':bindings})
 registry={s:{'new':True,'type':'scalar','output':True,'dimensions':dim if s in ['D','D0'] else [0]*7} for s in fields};registry['internalField']={'type':'internalField','file':'fields/internalField.yaml'}
 if stage:registry['boundaries']={'type':'boundary','file':'fields/boundaries.yaml'}
 write(p/'fields/fields.yaml','fields','registry',registry)
 initial={s:{'default':1 if s=='b' and stage else .2 if s=='b' else (.001 if stage else .1) if s=='D0' else .3 if s=='D' and target=='b' else 0} for s in fields}
 lines=['#Information',f'{n} 3 3 1','','#Point']
 for k in range(3):
  for j in range(3):
   for i in range(n):
    idx=(k*3+j)*n+i;lines.append(f'{i/(n-1):.17g} {j*.05:.17g} {k*.05:.17g}')
    if stage:initial['a']['point'+str(idx)]=math.sin(math.pi*i/(n-1))*math.sin(math.pi*j/2)*math.sin(math.pi*k/2)
 lines+=['','#all',' '.join(map(str,range(n*9)))]
 for idx in range(n*9):lines+=['','#point'+str(idx),str(idx)]
 for side,select in [('Left',lambda i,j,k:i==0),('Right',lambda i,j,k:i==n-1),('Bottom',lambda i,j,k:j==0),('Top',lambda i,j,k:j==2),('Front',lambda i,j,k:k==0),('Back',lambda i,j,k:k==2)]:lines+=['','#'+side,' '.join(str((k*3+j)*n+i) for k in range(3) for j in range(3) for i in range(n) if select(i,j,k))]
 write(p/'fields/internalField.yaml','fields','internalField',initial)
 if stage:write(p/'fields/boundaries.yaml','fields','boundaries',{'a':{side:{'type':'fixedValue','value':0} for side in ['Left','Right','Bottom','Top','Front','Back']}})
 write(p/'models/models.yaml','models','registry',{});write(p/'mesh/mesh.yaml','mesh','registry',{'files':['mesh/mesh.sfm']});(p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
 return p,n
records=[]
def run(p,action='run',expect=None):
 cmd=[str(binary),action,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_bytes(r.stdout)
 if expect is None:assert r.returncode==0,r.stdout[-4500:].decode()
 else:assert r.returncode and expect in r.stdout.decode(),r.stdout[-4500:].decode()
 return cmd

def output(p):
 tree=E.parse(next((p/'result').glob('*.pvd')));frames=tree.findall('.//DataSet');data=E.parse(p/'result'/frames[-1].get('file'));return {v.get('Name'):list(map(float,v.text.split())) for v in data.findall('.//PointData/DataArray')},float(frames[-1].get('timestep'))
for name,args,expected in [('direct',{},.12),('residual',{'residual':True},.12),('reverseTarget',{'target':'b'},2.),('renamed',{'renamed':True},.12),('loop',{'loop':True},12/7)]:
 p,n=make(name,**args);run(p,'explain');cmd=run(p);data,time=output(p);id='x' if name=='loop' else 'b' if name=='reverseTarget' else 'D';error=max(abs(v-expected) for v in data[id]);assert error<1e-14,(name,error);assert time==0
 if name=='loop':assert max(abs(v-20/7) for v in data['y'])<1e-14
 records.append({'case':name,'path':str(p),'command':cmd,'maxError':error,'time':time,'kind':'local affine/no physical-time advancement'});print(name,error,flush=True)
# The same nonlinear equation can still be affine for a different HOW target.
p,n=make('nonlinearOtherTarget',target='b')
f=p/'solvers/equationRegistry.yaml';v=json.loads(f.read_text());v['add']['D']={'lhs':op('multiply',op('divide',sym('D'),sym('D0')),op('divide',sym('D'),sym('D0'))),'rhs':op('add',num(1),sym('b'))};f.write_text(json.dumps(v,indent=2)+'\n')
run(p,'explain');cmd=run(p);data,time=output(p);error=max(abs(b-8) for b in data['b']);assert error<1e-14 and time==0
records.append({'case':p.name,'path':str(p),'command':cmd,'maxError':error,'time':time,'kind':'same nonlinear WHAT is affine for target b'});print(p.name,error,flush=True)
p,n=make('nonzeroClock')
f=p/'solvers/runtime.yaml';v=json.loads(f.read_text());v.update(startTime=3.,endTime=4.,writeControl="runTime",writeInterval=.1);f.write_text(json.dumps(v,indent=2)+'\n')
run(p,'explain');cmd=run(p);data,time=output(p);error=max(abs(d-.12) for d in data['D']);assert time==3 and error<1e-14
assert E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet')[-1].get('file').startswith('nonzeroClock_t')
records.append({'case':p.name,'path':str(p),'command':cmd,'maxError':error,'time':time,'kind':'algebraic solve preserves loaded physical clock'});print(p.name,error,time,flush=True)
errors=[]
for dt in [.05,.025,.0125]:
 p,n=make('stage'+str(dt),stage=True,dt=dt);run(p,'explain');cmd=run(p);data,time=output(p);lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2+1600
 error=max(abs(data['a'][(1*3+1)*n+i]-math.sin(math.pi*i/(n-1))*math.exp(-lam*.001*(time+(1-math.exp(-.7*time))/.7))) for i in range(n));relation=max(abs(d-.001*(1+b)) for d,b in zip(data['D'],data['b']));assert relation<1e-17 and abs(time-.2)<1e-14
 errors.append(error);records.append({'case':p.name,'path':str(p),'command':cmd,'maxError':error,'algebraicResidual':relation,'time':time,'dt':dt,'nodes':[n,3,3]});print(p.name,error,relation,flush=True)
assert errors[0]/errors[1]>13 and errors[1]/errors[2]>13,errors
# A StageLoop is a HOW construct, independent of the presence of an affine method.
p,n=make('stageWithoutAffine',stage=True,dt=.025)
def local_change(file,modify):
 f=p/file;v=json.loads(f.read_text());modify(v);f.write_text(json.dumps(v,indent=2)+'\n')
local_change('solvers/state.yaml',lambda v:v.update(use=['a','b']))
local_change('fields/fields.yaml',lambda v:v.pop('D'))
local_change('fields/internalField.yaml',lambda v:v.pop('D'))
local_change('solvers/how.yaml',lambda v:v['root']['children'][0].update(children=['G -> a','B -> b']))
local_change('solvers/which.yaml',lambda v:v.update(bindings=[{'equation':'G','method':'ExplicitCentral2Arithmetic'},{'equation':'B','method':'ScalarDiffusionCentral2'}]))
local_change('solvers/state.yaml',lambda v:v.update(calls=[{'equation':'G','reads':[{'symbol':'D0','version':'Current'}]}]))
local_change('solvers/equationRegistry.yaml',lambda v:(v.update(use=['G','B']),v['add'].pop('D'),v['add']['G'].update(rhs=op('add',op('diffusion',sym('D0'),sym('a')),op('multiply',op('multiply',num(.05,rate),sym('b')),sym('a'))))))
run(p,'explain');cmd=run(p);data,time=output(p);lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2+1600
error=max(abs(data['a'][(1*3+1)*n+i]-math.sin(math.pi*i/(n-1))*math.exp(-lam*.001*time+.05/.7*(1-math.exp(-.7*time)))) for i in range(n));assert error<1e-8,error
records.append({'case':p.name,'path':str(p),'command':cmd,'maxError':error,'time':time,'kind':'coupled explicit stage RHS without affine action'});print(p.name,error,flush=True)
(base/'evidence/numerical.json').write_text(json.dumps({'binary':hashlib.sha256(binary.read_bytes()).hexdigest(),'cases':records},indent=2))
negatives=[]
def change(p,file,edit):
 f=p/file;v=json.loads(f.read_text());edit(v);f.write_text(json.dumps(v,indent=2)+'\n')
def negative(name,source,edit,expected,action='run'):
 p=base/'negative'/name
 def gone(function,path,error):
  if not isinstance(error[1],FileNotFoundError):raise error[1]
 if p.exists():shutil.rmtree(p,onerror=gone)
 shutil.copytree(base/'cases'/source,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'));edit(p);cmd=run(p,action,expected);negatives.append({'case':name,'path':str(p),'command':cmd,'diagnostic':expected});print('rejected',name,flush=True)
negative('unknownEquation','direct',lambda p:change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['Unknown -> D'])),'Unknown')
negative('unknownTarget','direct',lambda p:change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['D -> Unknown'])),'Unknown')
# A field can be writable yet not occur in the equation.
negative('absentTarget','direct',lambda p:(change(p,'solvers/equationRegistry.yaml',lambda x:x['add']['D'].update(lhs=sym('D0'),rhs=num(.1,dim))),None),'target is absent')
negative('nonlinear','direct',lambda p:change(p,'solvers/equationRegistry.yaml',lambda x:x['add']['D'].update(lhs=op('multiply',op('divide',sym('D'),sym('D0')),op('divide',sym('D'),sym('D0'))),rhs=op('add',num(1),sym('b')))),'nonlinear')
negative('targetDenominator','direct',lambda p:change(p,'solvers/equationRegistry.yaml',lambda x:x['add']['D'].update(lhs=op('divide',sym('D0'),sym('D')),rhs=sym('b'))),'nonlinear')
negative('singular','direct',lambda p:(change(p,'solvers/equationRegistry.yaml',lambda x:x['add']['D'].update(lhs=op('subtract',sym('D'),sym('D')),rhs=num(0))),change(p,'solvers/state.yaml',lambda x:x.update(calls=[]))),'singular/non-finite')
negative('missingInitial','direct',lambda p:change(p,'fields/internalField.yaml',lambda x:x.pop('D0')),'initialization')
negative('wrongDimension','direct',lambda p:change(p,'fields/fields.yaml',lambda x:x['D0'].update(dimensions=[0]*7)),'Incompatible mathematical signatures')
negative('frozenInput','direct',lambda p:change(p,'solvers/state.yaml',lambda x:x['calls'][0]['reads'][0].update(version='Frozen')),'snapshot')
negative('wrongMethod','direct',lambda p:change(p,'solvers/which.yaml',lambda x:x['bindings'][0].update(method='ScalarDiffusionCentral2')),'explicit RHS')
negative('stageOrder','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root']['children'][0].update(children=['G -> a','D -> D','B -> b'])),'preceding Stage algebraic publication')
negative('missingStageProducer','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root']['children'][0].update(children=['G -> a','B -> b'])),'missing or ambiguous')
negative('nestedTemporalRhs','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root']['children'][0].update(children=[{'kind':'Loop','children':['G -> a'],'repetitions':2},'D -> D','B -> b'])),'temporal RHS cannot repeat')
negative('independentStageCount','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root']['children'][0].update(repetitions=7)),'recipe-owned Stage')
negative('implicitBlockRK4','stage0.025',lambda p:(change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['[G, B, D] => [a, b, D]'])),change(p,'solvers/which.yaml',lambda x:x.update(bindings=[{'equation':'*','method':'ScalarDiffusionCentral2'}])),change(p,'solvers/state.yaml',lambda x:x.update(calls=[{'call':'HOW/0','reads':[{'symbol':'D0','version':'Current'}]}]))),'whole block requires an explicit Step binding')
negative('blockMissingEquation','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['[Unknown, B, D] => [a, b, D]'])),'Unknown')
negative('blockDuplicateUnknown','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['[G, B, D] => [a, a, D]'])),'Duplicate block unknown')
negative('unknownProvider','direct',lambda p:change(p,'solvers/which.yaml',lambda x:x['bindings'][0].update(method='UnimplementedNewton')),'no registered implementation')
negative('localBoundary','direct',lambda p:(write(p/'fields/boundaries.yaml','fields','boundaries',{'D':{'Left':{'type':'fixedValue','value':0}}}),change(p,'fields/fields.yaml',lambda x:x.update(boundaries={'type':'boundary','file':'fields/boundaries.yaml'}))),'boundary/halo resource')
negative('stageTermination','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root']['children'][0].update(terminationSignal='arbitrarySignal')),'terminationSignal is an internal')
negative('unsupportedParallel','direct',lambda p:change(p,'solvers/runtime.yaml',lambda x:x.update(parallel={'enabled':True,'nProcs':2})),'serial',action='explain')
negative('emptyBlockMember','stage0.025',lambda p:change(p,'solvers/how.yaml',lambda x:x['root'].update(children=['[G, B, D,] => [a, b, D]'])),'Empty implicit block member')
(base/'evidence/negative.json').write_text(json.dumps(negatives,indent=2))
