"""A8 production native cases; independent analytic references, all evidence under test/."""
from pathlib import Path
import json,subprocess,os,math,hashlib,xml.etree.ElementTree as E,shutil,re,sys
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent;root=Path(__file__).resolve().parents[2];(base/'evidence').mkdir(parents=True,exist_ok=True);binary=root/'build/sonicSolver';env={**os.environ,'TMPDIR':str(base/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','SF_PLAN_TRACE':'1','SF_BLOCK_TRACE':'1'}
unit=[0]*7;rate=[0,0,-1,0,0,0,0];dim=[0,2,-1,0,0,0,0]
def sym(n):return {'symbol':n}
def num(n,d=None):return {'constant':n,**({'dimensions':d} if d is not None else {})}
def op(n,*args):return {'op':n,'args':list(args)}
def write(p,obj,id,**body):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps({'SonicFile':{'object':obj,'type':id},**body},indent=2)+'\n')
def change(p,file,edit):f=p/file;v=json.loads(f.read_text());edit(v);f.write_text(json.dumps(v,indent=2)+'\n')
def run(p,action='run',expect=None):
 cmd=[str(binary),action,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/(action+'.log')).write_bytes(r.stdout)
 if expect is None:assert r.returncode==0,r.stdout[-5000:].decode()
 else:assert r.returncode!=0 and expect in r.stdout.decode(),r.stdout[-5000:].decode()
 return {'path':str(p),'command':cmd,'returnCode':r.returncode,'binary':hashlib.sha256(binary.read_bytes()).hexdigest()}
def output(p):
 frames=E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');last=frames[-1];data=E.parse(p/'result'/last.get('file'));return {v.get('Name'):list(map(float,v.text.split())) for v in data.findall('.//PointData/DataArray')},float(last.get('timestep'))
def make(name,n=17,sine=False,renamed=False,reorder=False):
 p=base/'cases'/name;A,B,G,H=('u','v','first','second') if renamed else ('a','b','G','B');unknowns=[B,A] if reorder else [A,B];eqs=[H,G] if reorder else [G,H]
 a,b=sym(A),sym(B);da,db=.1,.2;alpha,beta,gamma,delta=2.,-.4,.3,1.5
 equations={G:{'lhs':op('add',op('add',op('negate',op('diffusion',num(da,dim),a)),op('multiply',num(alpha,rate),a)),op('multiply',num(beta,rate),b)),'rhs':sym('f')},H:{'lhs':op('add',op('add',op('negate',op('diffusion',num(db,dim),b)),op('multiply',num(gamma,rate),a)),op('multiply',num(delta,rate),b)),'rhs':sym('g')}}
 write(p/'case.yaml','case','registry',formatVersion=1,name=name)
 roles=[('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]
 write(p/'solvers/solvers.yaml','solver','registry',**{k:{'type':v,'file':'solvers/'+k+'.yaml'} for k,v in roles})
 write(p/'solvers/runtime.yaml','solver','runtime',startTime=3,endTime=4,maxDeltaT=.01,CFL=.5,writeControl='runTime',writeInterval=.1,writeInitial=False,parallel={'enabled':False},output={'jobName':name,'outputDir':'result'})
 write(p/'solvers/numerics.yaml','solver','numerics',time={'default':'classicalRK4'})
 write(p/'solvers/state.yaml','solver','state',use=unknowns,calls=[{'call':'coupled','reads':[{'symbol':name,'version':'Current'} for name in ['f','g']]}]);write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',use=eqs,add=equations)
 write(p/'solvers/execution.yaml','solver','execution',Explicit={});write(p/'solvers/how.yaml','solver','how',root={'kind':'Sequence','id':'HOW','children':[{'kind':'BlockSolve','id':'coupled','equations':eqs,'unknowns':unknowns}]})
 write(p/'solvers/which.yaml','solver','which',bindings=[{'equation':'*','occurrence':'coupled','method':'LinearBlockCentral2Hypre','parameters':{'relativeTolerance':1e-12,'absoluteTolerance':1e-13,'maxIterations':300,'krylovDimension':50}}])
 fields={s:{'new':True,'type':'scalar','output':True,'dimensions':unit if s in [A,B] else rate} for s in [A,B,'f','g']};fields.update(internalField={'type':'internalField','file':'fields/internalField.yaml'},boundaries={'type':'boundary','file':'fields/boundaries.yaml'});write(p/'fields/fields.yaml','fields','registry',**fields)
 initial={A:{'default':.3},B:{'default':.4},'f':{'default':0},'g':{'default':0}};points=[];refs={A:[],B:[]};ny=n if sine else 3;nz=ny
 for k in range(nz):
  for j in range(ny):
   for i in range(n):
    x,y,z=i/(n-1),j/(ny-1),k/(nz-1);idx=(k*ny+j)*n+i
    product=x*(1-x)*y*(1-y)*z*(1-z);lap=-2*(y*(1-y)*z*(1-z)+x*(1-x)*z*(1-z)+x*(1-x)*y*(1-y))
    av=math.sin(math.pi*x)*math.sin(math.pi*y)*math.sin(math.pi*z) if sine else 1+product;bv=.5*av if sine else 2+.5*product;d2a=-3*math.pi**2*av if sine else lap;d2b=-3*math.pi**2*bv if sine else .5*lap
    refs[A].append(av);refs[B].append(bv);initial['f']['p'+str(idx)]=-da*d2a+alpha*av+beta*bv;initial['g']['p'+str(idx)]=-db*d2b+gamma*av+delta*bv;points.append(f'{x:.17g} {y:.17g} {z:.17g}')
 write(p/'fields/internalField.yaml','fields','internalField',**initial);write(p/'fields/boundaries.yaml','fields','boundaries',**{s:{face:{'type':'fixedValue','value':0 if sine else 1 if s==A else 2} for face in ['Left','Right','Bottom','Top','Front','Back']} for s in [A,B]})
 lines=['#Information',f'{n} {ny} {nz} 1','','#Point',*points,'','#all',' '.join(map(str,range(n*ny*nz)))]
 for name,select in [('Left',lambda i,j,k:i==0),('Right',lambda i,j,k:i==n-1),('Bottom',lambda i,j,k:j==0),('Top',lambda i,j,k:j==ny-1),('Front',lambda i,j,k:k==0),('Back',lambda i,j,k:k==nz-1)]:lines+=['','#'+name,' '.join(str((k*ny+j)*n+i) for k in range(nz) for j in range(ny) for i in range(n) if select(i,j,k))]
 lines+=sum([['','#p'+str(i),str(i)] for i in range(n*ny*nz)],[])
 write(p/'models/models.yaml','models','registry');write(p/'mesh/mesh.yaml','mesh','registry',files=['mesh/mesh.sfm']);(p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
 return p,refs
records=[]
for name,kw in [('quadratic',{}),('renamed',{'renamed':True}),('reordered',{'reorder':True}),('arrow',{}),('equationOrder',{}),('columnOrder',{}),('sine5',{'n':5,'sine':True}),('sine9',{'n':9,'sine':True}),('sine17',{'n':17,'sine':True})]:
 p,refs=make(name,**kw)
 if name=='equationOrder':change(p,'solvers/how.yaml',lambda v:v['root']['children'][0].update(equations=['B','G']))
 if name=='columnOrder':change(p,'solvers/how.yaml',lambda v:v['root']['children'][0].update(unknowns=['b','a']))
 if name=='arrow':
  change(p,'solvers/how.yaml',lambda v:v['root'].update(children=['[G, B] => [a, b]']))
  change(p,'solvers/which.yaml',lambda v:v['bindings'][0].update(occurrence='HOW/0'))
  change(p,'solvers/state.yaml',lambda v:v['calls'][0].update(call='HOW/0'))
 run(p,'explain');record=run(p);data,time=output(p);err=max(abs(v-r) for s in refs for v,r in zip(data[s],refs[s]));assert time==3,(name,time)
 if not kw.get('sine'):assert err<2e-9,(name,err)
 n=kw.get('n',17);ny=n if kw.get('sine') else 3;nz=ny;N=n*ny*nz;fields=list(refs);physicalNodes=list(zip(data[fields[0]],data[fields[1]]))
 source=json.loads((p/'fields/internalField.yaml').read_text());residual=0.
 for k in range(1,nz-1):
  for j in range(1,ny-1):
   for i in range(1,n-1):
    cell=(k*ny+j)*n+i;lap=[]
    for target in fields:
     u=data[target];lap.append(sum((u[cell-offset]-2*u[cell]+u[cell+offset])/h**2 for offset,h in [(1,1/(n-1)),(n,1/(ny-1)),(n*ny,1/(nz-1))]))
    aa,bb=physicalNodes[cell];residual=max(residual,abs(-.1*lap[0]+2*aa-.4*bb-source['f']['p'+str(cell)]),abs(-.2*lap[1]+.3*aa+1.5*bb-source['g']['p'+str(cell)]))
 assert residual<2e-8,(name,residual)
 rows={}
 for line in (p/'run.log').read_text().splitlines():
  if line.startswith('[SF BLOCK ROW] '):
   parts=line.split();rows[int(parts[3])]={int(c.split(':')[0]):float(c.split(':')[1]) for c in parts[5:]}
 center=(1*ny+1)*n+n//2;aSlot=1 if kw.get('reorder') or name=='columnOrder' else 0;bSlot=1-aSlot
 assert rows[aSlot*N+center][bSlot*N+center]==-.4 and rows[bSlot*N+center][aSlot*N+center]==.3,'Missing off-diagonal block'
 record.update(case=name,maxError=err,independentResidual=residual,crossBlocks=[-.4,.3],time=time,nodes=[n,ny,nz]);records.append(record);print(name,err,residual,flush=True)
errors=[r['maxError'] for r in records if r['case'].startswith('sine')];assert errors[0]/errors[1]>3.8 and errors[1]/errors[2]>3.8,errors
# Reuse the pre-existing nodal zeroGradient mirror; the implicit row eliminates that same ghost.
naturalErrors=[]
for name,n,mixed in [('zero9',9,False),('zero17',17,False),('zero33',33,False),('mixedBlock',33,True)]:
 p,_=make(name,n=n);refs={'a':[],'b':[]};source=json.loads((p/'fields/internalField.yaml').read_text())
 for k in range(3):
  for j in range(3):
   for i in range(n):
    idx=(k*3+j)*n+i;x=i/(n-1);cos=math.cos(math.pi*x);a=1+.2*cos;b=2+.1*cos;refs['a'].append(a);refs['b'].append(b);source['f']['p'+str(idx)]=.02*math.pi**2*cos+2*a-.4*b;source['g']['p'+str(idx)]=.02*math.pi**2*cos+.3*a+1.5*b
 (p/'fields/internalField.yaml').write_text(json.dumps(source,indent=2)+'\n');bc={u:{face:{'type':'zeroGradient'} for face in ['Left','Right','Top','Bottom','Front','Back']} for u in ['a','b']}
 if mixed:
  for u,value in [('a',1.2),('b',2.1)]:bc[u]['Left']={'type':'fixedValue','value':value}
 write(p/'fields/boundaries.yaml','fields','boundaries',**bc);run(p,'explain');record=run(p);data,time=output(p);error=max(abs(value-ref) for u in refs for value,ref in zip(data[u],refs[u]));residual=0.
 for k in range(3):
  for j in range(3):
   for i in range(n):
    if mixed and i==0:continue
    cell=(k*3+j)*n+i;lap=[]
    for u in ['a','b']:
     total=0.
     for axis,(coord,size,stride,h) in enumerate([(i,n,1,1/(n-1)),(j,3,n,.5),(k,3,3*n,.5)]):
      for sign in [-1,1]:
       neighbor=coord+sign;step=sign if 0<=neighbor<size else -sign;total+=(data[u][cell+step*stride]-data[u][cell])/h**2
     lap.append(total)
    a,b=data['a'][cell],data['b'][cell];residual=max(residual,abs(-.1*lap[0]+2*a-.4*b-source['f']['p'+str(cell)]),abs(-.2*lap[1]+.3*a+1.5*b-source['g']['p'+str(cell)]))
 assert time==3 and residual<2e-8,(name,time,residual)
 if not mixed:naturalErrors.append(error)
 record.update(case=name,maxError=error,independentResidual=residual,time=time,nodes=[n,3,3],boundary='existing nodal mirror zeroGradient; left fixedValue only' if mixed else 'existing nodal mirror zeroGradient on all faces');records.append(record);print(name,error,residual,flush=True)
assert naturalErrors[0]/naturalErrors[1]>3.8 and naturalErrors[1]/naturalErrors[2]>3.8,naturalErrors
(base/'evidence/block-numerical.json').write_text(json.dumps(records,indent=2))
negatives=[]
def negative(name,edit,expected):
 p=base/'negative'/name
 if p.exists():shutil.rmtree(p,onerror=lambda f,p,e:None if isinstance(e[1],FileNotFoundError) else (_ for _ in ()).throw(e[1]))
 shutil.copytree(base/'cases/quadratic',p,ignore=shutil.ignore_patterns('result','*.log','._*'));edit(p);record=run(p,expect=expected);record.update(case=name,expected=expected);negatives.append(record);print('rejected',name,flush=True)
negative('unknownEquation',lambda p:change(p,'solvers/how.yaml',lambda v:v['root']['children'][0].update(equations=['Missing','B'])),'not registered')
negative('duplicateTarget',lambda p:change(p,'solvers/how.yaml',lambda v:v['root']['children'][0].update(unknowns=['a','a'])),'Duplicate block unknown')
negative('missingDof',lambda p:change(p,'solvers/how.yaml',lambda v:v['root']['children'][0].update(unknowns=['a'])),'square scalar')
negative('nonlinear',lambda p:change(p,'solvers/equationRegistry.yaml',lambda v:v['add']['G'].update(rhs=op('add',sym('f'),op('multiply',num(1,rate),op('multiply',sym('a'),sym('b')))))),'nonlinear block')
negative('nonconverged',lambda p:change(p,'solvers/which.yaml',lambda v:v['bindings'][0]['parameters'].update(maxIterations=1,relativeTolerance=1e-15,absoluteTolerance=0)),'did not converge')
negative('transientBlock',lambda p:change(p,'solvers/equationRegistry.yaml',lambda v:v['add']['G'].update(lhs=op('add',v['add']['G']['lhs'],op('ddt',sym('a'))))),'Time strategy conflict')
negative('missingBoundary',lambda p:change(p,'fields/boundaries.yaml',lambda v:v['a'].pop('Right')),'Missing scalar boundary')
negative('wrongDimensions',lambda p:change(p,'fields/fields.yaml',lambda v:v['f'].update(dimensions=unit)),'signature')
negative('laggedInput',lambda p:change(p,'solvers/state.yaml',lambda v:v['calls'][0].update(reads=[{'symbol':'f','version':'Lagged'},{'symbol':'g','version':'Current'}])),'input version')
negative('explicitStageBlock',lambda p:change(p,'solvers/how.yaml',lambda v:v['root'].update(children=[{'kind':'StageLoop','children':v['root']['children']}])),'Time strategy conflict')
negative('unsupportedParallel',lambda p:change(p,'solvers/runtime.yaml',lambda v:v.update(parallel={'enabled':True,'nProcs':2})),'serial single patch')
(base/'evidence/block-negative.json').write_text(json.dumps(negatives,indent=2))
