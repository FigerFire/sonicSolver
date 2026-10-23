"""Native IO -> four contracts -> actual binder/Plan; independent spectral/analytic references."""
from pathlib import Path
import math,json,os,subprocess,xml.etree.ElementTree as ET,hashlib,shutil,sys
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'TMPDIR':str(base/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','SF_SCALAR_TRACE':'1','SF_PLAN_TRACE':'1'}
binary=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else root/'build/sonicSolver'
def write(p,obj,identity,body):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps({'SonicFile':{'object':obj,'type':identity},**body},indent=2)+'\n')
def op(name,*args):return {'op':name,'args':list(args)}
def symbol(name):return {'symbol':name}
def coefficient(v,dim):return {'constant':v,'dimensions':dim}
rate=lambda v:coefficient(v,[0,0,-1,0,0,0,0])
D=lambda v,x:op('diffusion',coefficient(v,[0,2,-1,0,0,0,0]),symbol(x))
def make(name):
 p=base/'native'/name;stencil=name not in ['reaction','logistic'];vector=name in ['vector','vectorFixed'];ids=['Vexperiment'] if vector else ['a','b'] if name in ['reaction','cross'] else ['a'];n=33 if stencil else 5;g=1;end=.05 if stencil else 1.;dt=.0001 if stencil else .025
 rhs={'reaction':[op('multiply',rate(.7),symbol('b')),op('multiply',rate(.3),symbol('a'))],
      'logistic':[op('multiply',op('multiply',rate(2),symbol('a')),op('subtract',{'constant':1},symbol('a')))],
      'multi':[op('add',D(.03,'a'),D(.07,'a'))],
      'cross':[op('add',D(.07,'a'),D(.02,'b')),op('add',D(.02,'a'),D(.04,'b'))],
      'vectorFixed':[D(.03,'Vexperiment')],
      'vector':[op('subtract',D(.03,'Vexperiment'),op('multiply',rate(.3),symbol('Vexperiment')))]}[name]
 write(p/'case.yaml','case','registry',{'formatVersion':1,'name':name})
 roles=[('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]
 write(p/'solvers/solvers.yaml','solver','registry',{k:{'type':r,'file':'solvers/'+k+'.yaml'} for k,r in roles})
 write(p/'solvers/runtime.yaml','solver','runtime',{'startTime':0,'endTime':end,'maxDeltaT':dt,'CFL':.5,'writeControl':'timeStep','writeInterval':500 if stencil else 40,'writeInitial':True,'parallel':{'enabled':False},'output':{'jobName':name,'outputDir':'result'}})
 write(p/'solvers/numerics.yaml','solver','numerics',{'time':{'default':'classicalRK4'}})
 write(p/'solvers/state.yaml','solver','state',{'use':ids});eqs=['equation.'+s for s in ids]
 write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',{'use':eqs,'add':{e:{'lhs':op('ddt',symbol(s)),'rhs':r} for s,e,r in zip(ids,eqs,rhs)}})
 write(p/'solvers/execution.yaml','solver','execution',{'Explicit':{}})
 write(p/'solvers/how.yaml','solver','how',{'root':{'kind':'Sequence','id':'A6','children':[{'kind':'EquationCall','equation':e,'target':s,'occurrence':'advance.'+s} for s,e in zip(ids,eqs)]+[{'kind':'Commit','id':'commit'}]}})
 write(p/'solvers/which.yaml','solver','which',{'bindings':[{'equation':e,'method':'ScalarDiffusionCentral2'} for e in eqs]})
 fields={s:{'new':True,'type':'vector' if vector else 'scalar','output':True,'dimensions':[0,1,-1,0,0,0,0] if vector else [0]*7} for s in ids};fields['internalField']={'type':'internalField','file':'fields/internalField.yaml'}
 if stencil:fields['boundaries']={'type':'boundary','file':'fields/boundaries.yaml'}
 write(p/'fields/fields.yaml','fields','registry',fields)
 initial={}
 for q,s in enumerate(ids):
  initial[s]={'default':([1,-2,.5] if name=='vectorFixed' else [0,0,0] if vector else .2 if name=='logistic' else 2 if q else 1)}
  if name=='vectorFixed':initial[s]['all']=[1,-2,.5]
  if stencil and name!='vectorFixed':
   for i in range(n):initial[s]['node'+str(i)]=[amp*math.sin((c+1)*math.pi*i/(n-1)) for c,amp in enumerate([1,-2,.5])] if vector else (q+1)*math.sin(math.pi*i/(n-1))
 write(p/'fields/internalField.yaml','fields','internalField',initial)
 if stencil:write(p/'fields/boundaries.yaml','fields','boundaries',{s:{**{side:{'type':'fixedValue','value':[1,-2,.5] if name=='vectorFixed' else [0,0,0] if vector else 0} for side in ['Left','Right']},**{side:{'type':'zeroGradient'} for side in ['Bottom','Top','Front','Back']}} for s in ids})
 write(p/'models/models.yaml','models','registry',{});write(p/'mesh/mesh.yaml','mesh','registry',{'files':['mesh/mesh.sfm']})
 lines=['#Information',f'{n} 3 3 {g}','','#Point']+[f'{i/(n-1):.17g} {j*.05:.17g} {k*.05:.17g}' for k in range(3) for j in range(3) for i in range(n)]
 lines+=['','#all',' '.join(map(str,range(n*9)))]
 for i in range(n):lines+=['','#node'+str(i),' '.join(str((k*3+j)*n+i) for k in range(3) for j in range(3))]
 for side,test in [('Left',lambda i,j,k:i==0),('Right',lambda i,j,k:i==n-1),('Bottom',lambda i,j,k:j==0),('Top',lambda i,j,k:j==2),('Front',lambda i,j,k:k==0),('Back',lambda i,j,k:k==2)]:lines+=['','#'+side,' '.join(str((k*3+j)*n+i) for k in range(3) for j in range(3) for i in range(n) if test(i,j,k))]
 (p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
 return p,n,end,ids
records=[]
for name in ['reaction','logistic','multi','cross','vector','vectorFixed']:
 p,n,t,ids=make(name)
 for action in ['explain','run']:
  cmd=[str(binary),action,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_bytes(r.stdout)
  if r.returncode:raise RuntimeError(name+' '+action+' '+r.stdout[-2500:].decode())
 frames=ET.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=p/'result'/frames[-1].get('file');tree=ET.parse(last)
 errors=[]
 for q,s in enumerate(ids):
  field=tree.find('.//PointData/DataArray[@Name="'+s+'"]');count=int(field.get('NumberOfComponents','1'));values=list(map(float,field.text.split()));assert count==(3 if name in ['vector','vectorFixed'] else 1)
  for i in range(len(values)//count):
   x=(i%n)/(n-1);lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2
   for c in range(count):
    if name=='reaction':z=math.sqrt(.21);ref=math.cosh(z)+1.4/z*math.sinh(z) if q==0 else 2*math.cosh(z)+.3/z*math.sinh(z)
    elif name=='logistic':ref=1/(1+4*math.exp(-2))
    elif name=='vectorFixed':ref=[1,-2,.5][c]
    elif name=='multi':ref=math.sin(math.pi*x)*math.exp(-.1*lam*t)
    elif name=='cross':ch=math.cosh(lam*t*.025);sh=math.sinh(lam*t*.025)/.025;ref=math.sin(math.pi*x)*math.exp(-.055*lam*t)*(ch-sh*.055 if q==0 else 2*ch+sh*.01)
    else:
     lk=4*math.sin((c+1)*math.pi/(2*(n-1)))**2*(n-1)**2;ref=[1,-2,.5][c]*math.sin((c+1)*math.pi*x)*math.exp(-(.03*lk+.3)*t)
    errors.append(abs(values[i*count+c]-ref))
 error=max(errors);assert error<(3e-8 if name=='reaction' else 5e-8 if name=='logistic' else 1e-12),(name,error)
 assert abs(float(frames[-1].get('timestep'))-t)<1e-12
 records.append(dict(name=name,path=str(p),binary=hashlib.sha256(binary.read_bytes()).hexdigest(),maxError=error,components=count,finalTime=t,frames=len(frames),command=[str(binary),'run',str(p)]));print(name,error,flush=True)
(base/'evidence/native.json').write_text(json.dumps(records,indent=2))
