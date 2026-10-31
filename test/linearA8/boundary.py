"""Existing nodal ghost reuse: native scalar/vector, conservation and analytic convergence."""
from pathlib import Path
import json,math,os,subprocess,hashlib,xml.etree.ElementTree as E,shutil,sys
root=Path(__file__).resolve().parents[2];base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'boundary';binary=root/'build/sonicSolver';(base/'evidence').mkdir(parents=True,exist_ok=True)
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
unit=[0]*7;rate=[0,0,-1,0,0,0,0];dim=[0,2,-1,0,0,0,0]
def sym(n):return {'symbol':n}
def num(n,d):return {'constant':n,'dimensions':d}
def op(n,*args):return {'op':n,'args':list(args)}
def write(p,obj,id,**body):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps({'SonicFile':{'object':obj,'type':id},**body},indent=2)+'\n')
def change(p,file,edit):f=p/file;v=json.loads(f.read_text());edit(v);f.write_text(json.dumps(v,indent=2)+'\n')
def make(name,n=17,vector=False,constant=False,legacy=False,mms=False,mixed=False,noSource=False):
 p=base/'cases'/name;nc=3 if vector else 1;end=.05 if mms or mixed else .001;dt=.0001 if mms or mixed else .00005
 roles=[('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]
 write(p/'case.yaml','case','registry',formatVersion=1,name=name);write(p/'solvers/solvers.yaml','solver','registry',**{k:{'type':v,'file':'solvers/'+k+'.yaml'} for k,v in roles})
 write(p/'solvers/runtime.yaml','solver','runtime',startTime=0,endTime=end,maxDeltaT=dt,CFL=.5,writeControl='timeStep',writeInterval=1 if not mms and not mixed else 500,writeInitial=True,parallel={'enabled':False},output={'jobName':name,'outputDir':'result'})
 write(p/'solvers/numerics.yaml','solver','numerics',time={'default':'classicalRK4'});write(p/'solvers/state.yaml','solver','state',use=['u'],**({} if legacy else {'calls':[{'equation':'E','reads':[{'symbol':name,'version':'Current'} for name in (['D'] if noSource else ['D','S'])]}]}))
 rhs=op('diffusion',num(.04,dim) if legacy else sym('D'),sym('u'))
 if not noSource:rhs=op('add',rhs,sym('S'))
 write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',use=['E'],add={'E':{'lhs':op('ddt',sym('u')),'rhs':rhs}})
 write(p/'solvers/execution.yaml','solver','execution',Explicit={});write(p/'solvers/how.yaml','solver','how',root={'kind':'Sequence','id':'HOW','children':['E -> u',{'kind':'Commit','id':'commit'}]})
 # The old constant provider does not accept authored version overrides. Constant legacy uses a literal source too.
 if legacy:
  write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',use=['E'],add={'E':{'lhs':op('ddt',sym('u')),'rhs':op('add',op('diffusion',num(.04,dim),sym('u')),num(.07,rate))}})
 write(p/'solvers/which.yaml','solver','which',bindings=[{'equation':'E','method':'ScalarDiffusionCentral2' if legacy else 'ExplicitCentral2Arithmetic'}])
 fields={'u':{'new':True,'type':'vector' if vector else 'scalar','dimensions':unit,'output':True}}
 if not legacy:fields.update(D={'new':True,'type':'scalar','dimensions':dim,'output':True},S={'new':True,'type':'vector' if vector else 'scalar','dimensions':rate,'output':True})
 if noSource:fields.pop('S',None)
 fields.update(internalField={'type':'internalField','file':'fields/internalField.yaml'},boundaries={'type':'boundary','file':'fields/boundaries.yaml'});write(p/'fields/fields.yaml','fields','registry',**fields)
 ic={key:{'default':[0]*nc if vector and key in ['u','S'] else 0} for key in fields if key not in ['internalField','boundaries']};points=[];initial=[];coefficients=[];source=[]
 for k in range(3):
  for j in range(3):
   for i in range(n):
    x,y,z=i/(n-1),j/2,k/2;idx=(k*3+j)*n+i;D=.04 if constant or legacy else .04+.01*math.cos(math.pi*x) if mms or mixed else .03+.01*x+.005*y+.003*z
    if mms or mixed:
     wave=math.pi/2 if mixed else math.pi;u=math.cos(wave*x);s=D*wave**2*u-.01*math.pi*wave*math.sin(math.pi*x)*math.sin(wave*x)
     if mixed and i==n-1:u=0
     uv=[u]*nc;sv=[s]*nc
    else:
     uv=[1+c+(.3+.1*c)*math.cos(math.pi*x)+(.2-.03*c)*math.cos(math.pi*y)+(.1+.02*c)*math.cos(math.pi*z) for c in range(nc)];sv=[0]*nc if noSource else [.07,-.02,.04][:nc]
    ic['u']['p'+str(idx)]=uv if vector else uv[0];initial.append(uv);coefficients.append(D);source.append(sv)
    if not legacy:ic['D']['p'+str(idx)]=D
    if not legacy and not noSource:ic['S']['p'+str(idx)]=sv if vector else sv[0]
    points.append(f'{x:.17g} {y:.17g} {z:.17g}')
 write(p/'fields/internalField.yaml','fields','internalField',**ic)
 bc={face:{'type':'zeroGradient'} for face in ['Left','Right','Bottom','Top','Front','Back']}
 if mixed:bc['Right']={'type':'fixedValue','value':[0]*nc if vector else 0}
 write(p/'fields/boundaries.yaml','fields','boundaries',u=bc);write(p/'models/models.yaml','models','registry');write(p/'mesh/mesh.yaml','mesh','registry',files=['mesh/mesh.sfm'])
 lines=['#Information',f'{n} 3 3 1','','#Point',*points,'','#all',' '.join(map(str,range(n*9)))]
 for name,select in [('Left',lambda i,j,k:i==0),('Right',lambda i,j,k:i==n-1),('Bottom',lambda i,j,k:j==0),('Top',lambda i,j,k:j==2),('Front',lambda i,j,k:k==0),('Back',lambda i,j,k:k==2)]:lines+=['','#'+name,' '.join(str((k*3+j)*n+i) for k in range(3) for j in range(3) for i in range(n) if select(i,j,k))]
 lines+=sum([['','#p'+str(i),str(i)] for i in range(n*9)],[]);(p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
 weights=[(1/(n-1))*.25*(.5 if i in [0,n-1] else 1)*(.5 if j in [0,2] else 1)*(.5 if k in [0,2] else 1) for k in range(3) for j in range(3) for i in range(n)]
 return p,dict(n=n,components=nc,initial=initial,D=coefficients,S=source,weights=weights,end=end,dt=dt,constant=constant or legacy,mms=mms,mixed=mixed)
def run(p,exe=binary,expect=None):
 cmd=[str(exe),'run',str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_bytes(r.stdout)
 if expect is None:assert r.returncode==0,r.stdout[-4000:].decode()
 else:assert r.returncode and expect in r.stdout.decode(),r.stdout[-4000:].decode()
 return dict(command=cmd,binary=hashlib.sha256(exe.read_bytes()).hexdigest(),returnCode=r.returncode)
def frames(p):
 for frame in E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet'):
  data=E.parse(p/'result'/frame.get('file'));yield float(frame.get('timestep')),{v.get('Name'):list(map(float,v.text.split())) for v in data.findall('.//PointData/DataArray')}
records=[];last={}
for name,kw in [('constantLegacy',{'legacy':True}),('constantNew',{'constant':True}),('budgetNoSource',{'noSource':True}),('budgetScalar',{}),('budgetVector',{'vector':True}),('mms17',{'n':17,'mms':True}),('mms33',{'n':33,'mms':True}),('mms65',{'n':65,'mms':True}),('mixed',{'n':33,'mixed':True})]:
 p,settings=make(name,**kw);record=run(p,root/'test/linearA8/boundary/baseline/sonicSolver' if kw.get('legacy') else binary);n=settings['n'];nc=settings['components'];maximum=0.;budget=0.;observed=list(frames(p));last[name]=observed[-1][1]['u']
 initialMass=[sum(w*u[c] for w,u in zip(settings['weights'],settings['initial'])) for c in range(nc)];sourceMass=[sum(w*s[c] for w,s in zip(settings['weights'],settings['S'])) for c in range(nc)]
 for time,data in observed:
  u=[data['u'][i*nc:(i+1)*nc] for i in range(n*9)]
  if not settings['mixed']:
   for c in range(nc):budget=max(budget,abs(sum(w*v[c] for w,v in zip(settings['weights'],u))-initialMass[c]-time*sourceMass[c]))
  if 'D' in data:assert data['D']==settings['D'],'Physical coefficient overwritten'
  if settings['mms'] or settings['mixed']:
   maximum=max(maximum,max(abs(v[0]-settings['initial'][i][0]) for i,v in enumerate(u)))
  elif settings['constant']:
   for i,v in enumerate(u):
    x=(i%n)/(n-1);y=(i//n)%3/2;z=(i//(n*3))/2;lam=4*math.sin(math.pi/(2*(n-1)))**2*(n-1)**2
    for c in range(nc):
     ref=1+c+settings['S'][0][c]*time+(.3+.1*c)*math.cos(math.pi*x)*math.exp(-.04*lam*time)+(.2-.03*c)*math.cos(math.pi*y)*math.exp(-.04*8*time)+(.1+.02*c)*math.cos(math.pi*z)*math.exp(-.04*8*time);maximum=max(maximum,abs(v[c]-ref))
 assert abs(observed[-1][0]-settings['end'])<1e-14
 if not settings['mixed']:assert budget<2e-13,(name,budget)
 if settings['constant']:assert maximum<2e-13,(name,maximum)
 record.update(case=name,settings=settings,frames=len(observed),finalTime=observed[-1][0],maxError=maximum if settings['mms'] or settings['mixed'] or settings['constant'] else None,budgetError=budget if not settings['mixed'] else None);records.append(record);print(name,maximum,budget,flush=True)
legacyDifference=max(abs(a-b) for a,b in zip(last['constantLegacy'],last['constantNew']));assert legacyDifference<2e-14
(base/'evidence/legacy-equivalence.json').write_text(json.dumps({'maxFieldDifference':legacyDifference,'baselineBinary':hashlib.sha256((root/'test/linearA8/boundary/baseline/sonicSolver').read_bytes()).hexdigest()},indent=2))
errors=[r['maxError'] for r in records if r['case'].startswith('mms')];assert errors[0]/errors[1]>3.5 and errors[1]/errors[2]>3.5,errors
(base/'evidence/numerical.json').write_text(json.dumps(records,indent=2))
negative=[]
for name,edit,expected in [('missingBoundary',lambda p:change(p,'fields/boundaries.yaml',lambda v:v['u'].pop('Left')),'Missing scalar boundary'),('normalFlux',lambda p:change(p,'fields/boundaries.yaml',lambda v:v['u']['Left'].update(type='normalFlux',value=.1)),'fixedValue/zeroGradient only'),('wrongDimensions',lambda p:change(p,'fields/fields.yaml',lambda v:v['D'].update(dimensions=unit)),'signature')]:
 p=base/'negative'/name
 if p.exists():shutil.rmtree(p,onerror=lambda f,p,e:None if isinstance(e[1],FileNotFoundError) else (_ for _ in ()).throw(e[1]))
 shutil.copytree(base/'cases/budgetScalar',p,ignore=shutil.ignore_patterns('result','*.log','._*'));edit(p);record=run(p,expect=expected);record.update(case=name,expected=expected);negative.append(record);print('rejected',name,flush=True)
(base/'evidence/negative.json').write_text(json.dumps(negative,indent=2))
