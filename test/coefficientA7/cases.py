"""A7 production native IO and independent semidiscrete reference; outputs stay in test/."""
from pathlib import Path
import hashlib,json,math,os,subprocess,shutil,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent;root=base.parents[1];binary=root/'build/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
def write(p,obj,identity,body):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps({'SonicFile':{'object':obj,'type':identity},**body},indent=2)+'\n')
def op(name,*args):return {'op':name,'args':list(args)}
def sym(name):return {'symbol':name}
p=base/'native/given';n=33;h=1/(n-1);dt=1e-5;end=.001
write(p/'case.yaml','case','registry',{'formatVersion':1,'name':'given'})
roles=[('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]
write(p/'solvers/solvers.yaml','solver','registry',{k:{'type':r,'file':'solvers/'+k+'.yaml'} for k,r in roles})
write(p/'solvers/runtime.yaml','solver','runtime',{'startTime':0,'endTime':end,'CFL':.5,'maxDeltaT':dt,'writeControl':'timeStep','writeInterval':100,'writeInitial':True,'parallel':{'enabled':False},'output':{'jobName':'given','outputDir':'result'}})
write(p/'solvers/numerics.yaml','solver','numerics',{'time':{'default':'classicalRK4'}})
write(p/'solvers/state.yaml','solver','state',{'use':['a'],'calls':[{'equation':'balance','reads':[{'symbol':'D','version':'Current'}]}]})
write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',{'use':['balance'],'add':{'balance':{'lhs':op('ddt',sym('a')),'rhs':op('diffusion',sym('D'),sym('a'))}}})
write(p/'solvers/execution.yaml','solver','execution',{'Explicit':{}})
write(p/'solvers/how.yaml','solver','how',{'root':{'kind':'Sequence','id':'A7','children':[{'kind':'EquationCall','equation':'balance','target':'a','occurrence':'advance.a'},{'kind':'Commit','id':'commit'}]}})
write(p/'solvers/which.yaml','solver','which',{'bindings':[{'equation':'balance','method':'ExplicitCentral2Arithmetic'}]})
write(p/'fields/fields.yaml','fields','registry',{'a':{'new':True,'type':'scalar','output':True,'dimensions':[0]*7},'D':{'new':True,'type':'scalar','output':True,'dimensions':[0,2,-1,0,0,0,0]},'internalField':{'type':'internalField','file':'fields/internalField.yaml'},'boundaries':{'type':'boundary','file':'fields/boundaries.yaml'}})
ic={'a':{'default':0},'D':{'default':0}}
lines=['#Information',f'{n} 3 3 1','','#Point']
for k in range(3):
 for j in range(3):
  for i in range(n):
   idx=(k*3+j)*n+i;x=i*h;y=j*.05;z=k*.05;lines.append(f'{x:.17g} {y:.17g} {z:.17g}')
   ic['a']['point'+str(idx)]=math.sin(math.pi*x)*math.sin(math.pi*j/2)*math.sin(math.pi*k/2)
   ic['D']['point'+str(idx)]=.03+.01*x
lines+=['','#all',' '.join(map(str,range(n*9)))]
for idx in range(n*9):lines+=['','#point'+str(idx),str(idx)]
for side,test in [('Left',lambda i,j,k:i==0),('Right',lambda i,j,k:i==n-1),('Bottom',lambda i,j,k:j==0),('Top',lambda i,j,k:j==2),('Front',lambda i,j,k:k==0),('Back',lambda i,j,k:k==2)]:
 lines+=['','#'+side,' '.join(str((k*3+j)*n+i) for k in range(3) for j in range(3) for i in range(n) if test(i,j,k))]
write(p/'fields/internalField.yaml','fields','internalField',ic)
write(p/'fields/boundaries.yaml','fields','boundaries',{'a':{side:{'type':'fixedValue','value':0} for side in ['Left','Right','Bottom','Top','Front','Back']}})
write(p/'models/models.yaml','models','registry',{});write(p/'mesh/mesh.yaml','mesh','registry',{'files':['mesh/mesh.sfm']});(p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
def run(p,action='run'):
 cmd=[str(binary),action,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_bytes(r.stdout);return r,cmd
for action in ['explain','run']:
 r,cmd=run(p,action);assert r.returncode==0,r.stdout[-3000:].decode()
# Independent e^(t A) on the x interior; y/z each have one nonzero Dirichlet interior node.
D=[.03+.01*i*h for i in range(n)];initial=[math.sin(math.pi*i*h) for i in range(1,n-1)]
def apply(v):
 v=[0]+v+[0];return [( .5*(D[i]+D[i+1])*(v[i+1]-v[i])-.5*(D[i-1]+D[i])*(v[i]-v[i-1]))/h**2-1600*D[i]*v[i] for i in range(1,n-1)]
reference=initial[:];term=initial[:]
for k in range(1,100):
 term=[end*v/k for v in apply(term)];reference=[a+b for a,b in zip(reference,term)]
 if max(map(abs,term))<1e-19:break
else:raise RuntimeError('Reference series did not converge')
frames=ET.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=p/'result'/frames[-1].get('file');tree=ET.parse(last)
a=list(map(float,tree.find('.//PointData/DataArray[@Name="a"]').text.split()));actualD=list(map(float,tree.find('.//PointData/DataArray[@Name="D"]').text.split()))
error=max(abs(a[(1*3+1)*n+i]-reference[i-1]) for i in range(1,n-1));dError=max(abs(v-D[i%n]) for i,v in enumerate(actualD));assert error<1e-12,(error,dError);assert dError<1e-15
assert abs(float(frames[-1].get('timestep'))-end)<1e-12
record={'path':str(p),'binary':hashlib.sha256(binary.read_bytes()).hexdigest(),'command':cmd,'endTime':end,'dt':dt,'steps':100,'nodes':[n,3,3],'reference':'independent matrix exponential (Taylor)','maxError':error,'coefficientError':dError,'seriesTerms':k}
(base/'evidence/native.json').write_text(json.dumps(record,indent=2));print('native given',error,'D unchanged',dError,flush=True)
negatives=[]
def negative(name,edit,expected,action="run"):
 q=base/'negative'/name
 if q.exists():
  def gone(function,path,error):
   if not isinstance(error[1],FileNotFoundError):raise error[1]
  shutil.rmtree(q,onerror=gone)
 shutil.copytree(p,q,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'));edit(q);r,cmd=run(q,action)
 assert '[SF] Scalar step:' not in r.stdout.decode(),name
 assert r.returncode!=0 and expected in r.stdout.decode(),(name,r.returncode,r.stdout[-2500:].decode());negatives.append({'name':name,'path':str(q),'command':cmd,'exit':r.returncode,'diagnostic':expected});print('rejected',name,flush=True)
def change(q,file,edit):
 f=q/file;v=json.loads(f.read_text());edit(v);f.write_text(json.dumps(v,indent=2)+'\n')
for version in ['Stage','Frozen','Lagged','OldTime']:
 negative('version'+version,lambda q,v=version:change(q,'solvers/state.yaml',lambda x:x['calls'][0]['reads'][0].update(version=v)),'Stage publisher' if version=='Stage' else 'requires a real snapshot')
negative('missingRead',lambda q:change(q,'solvers/state.yaml',lambda x:x.update(calls=[])),"unknown declared STATE symbol 'D'")
negative('missingField',lambda q:change(q,'fields/fields.yaml',lambda x:x.pop('D')),'D')
negative('missingInitial',lambda q:change(q,'fields/internalField.yaml',lambda x:x.pop('D')),'D')
negative('oldMethod',lambda q:change(q,'solvers/which.yaml',lambda x:x['bindings'][0].update(method='ScalarDiffusionCentral2')),'Unsupported scalar diffusivity')
negative('unsupportedBoundary',lambda q:change(q,'fields/boundaries.yaml',lambda x:x['a']['Left'].update(type='normalFlux')),'fixedValue/zeroGradient only')
negative('duplicateRead',lambda q:change(q,'solvers/state.yaml',lambda x:x['calls'][0]['reads'].append({'symbol':'D','version':'Current'})),'Duplicate STATE input read')
negative('unknownVersion',lambda q:change(q,'solvers/state.yaml',lambda x:x['calls'][0]['reads'][0].update(version='Newest')),'version')
negative('unsupportedTime',lambda q:change(q,'solvers/numerics.yaml',lambda x:x['time'].update(default='SSPRK3')),'cannot compile occurrence')
negative('unsupportedBackend',lambda q:change(q,'solvers/runtime.yaml',lambda x:x.update(parallel={'enabled':True,'nProcs':2})),'serial',action='explain')
(base/'evidence/negative.json').write_text(json.dumps(negatives,indent=2))
