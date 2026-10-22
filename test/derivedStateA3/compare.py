from pathlib import Path
import xml.etree.ElementTree as E,hashlib,json,math,re
base=Path(__file__).resolve().parent;root=base.parents[1];records=[]
def data(p):
 tree=E.parse(p);extent=list(map(int,tree.find('.//StructuredGrid').attrib['WholeExtent'].split()));n=[extent[2*i+1]-extent[2*i]+1 for i in range(3)]
 pts=list(map(float,tree.find('.//Points/DataArray').text.split()));coords=[pts[i::3] for i in range(3)]
 axes=[sorted(set(a)) for a in coords];volume=math.prod(a[-1]-a[0] for a in axes);step=math.prod((a[-1]-a[0])/(len(a)-1) for a in axes)
 weights=[step*math.prod(.5 if a==b[0] or a==b[-1] else 1. for a,b in zip(c,axes)) for c in zip(*coords)]
 fields={a.attrib['Name']:(int(a.get('NumberOfComponents','1')),list(map(float,a.text.split()))) for a in tree.findall('.//PointData/DataArray')}
 stats={}
 for key,(components,values) in fields.items():
  if key.startswith('BC_'):continue
  stats[key]=dict(min=min(values),max=max(values),weightedL2=math.sqrt(math.fsum(w*math.fsum(v*v for v in values[i*components:(i+1)*components]) for i,w in enumerate(weights))/volume))
 def integral(values):return math.fsum(x*w for x,w in zip(values,weights))
 totals={}
 if 'rho' in fields:totals['mass']=integral(fields['rho'][1])
 for key in ['MomentumX','MomentumY','MomentumZ','TotalEnergyDensity']:
  if key in fields:totals[key]=integral(fields[key][1])
 if 'rho' in fields and 'MomentumX' not in fields and 'U' in fields:
  totals['momentumFromRhoU']=[integral([r*u for r,u in zip(fields['rho'][1],fields['U'][1][c::3])]) for c in range(3)]
 for key in list(fields):
  if key.startswith('phaseMass.'):
   phase=key.split('.',1)[1];mass=fields[key][1];totals[key]=integral(mass)
   totals['phaseMomentum.'+phase]=[integral([r*u for r,u in zip(mass,fields['U.'+phase][1][c::3])]) for c in range(3)]
 return dict(nodes=n,domain=[[a[0],a[-1]] for a in axes],quadrature='tensor-product trapezoidal nodal physical-volume weights',totals=totals,fields=stats)
for name in ['sod','stageSST','pressure','conservativePressure','eulerianSST','pairRK4','linearTermsRK4']:
 a=base/'before'/name;b=base/'after'/name
 files=sorted(str(p.relative_to(a/'result')) for p in (a/'result').rglob('*') if p.is_file() and not p.name.startswith('._'))
 assert files==sorted(str(p.relative_to(b/'result')) for p in (b/'result').rglob('*') if p.is_file() and not p.name.startswith('._'))
 hashes={p:hashlib.sha256((b/'result'/p).read_bytes()).hexdigest() for p in files}
 assert all((a/'result'/p).read_bytes()==(b/'result'/p).read_bytes() for p in files),name
 traces=lambda p:[line for line in (p/'run.log').read_text().splitlines() if line.startswith(('[SF PLAN]','[SF SCALAR]','[SF] Step time:'))]
 assert traces(a)==traces(b),name
 providers=lambda p:[line for line in (p/'explain.log').read_text().splitlines() if any(x in line for x in ['WHICH:','provider=','occurrence='])]
 assert providers(a)==providers(b),(name,'provider/occurrence differ')
 pvd=next((b/'result').glob('*.pvd'));times=[x.attrib['timestep'] for x in E.parse(pvd).findall('.//DataSet')]
 last=sorted((b/'result').glob('*.vts'))[-1];stats=data(last)
 assert stats==data(a/'result'/last.name)
 record=dict(case=name,files=len(files),frames=len(times),finalTime=times[-1],changed=[],traceEqual=True,traceLines=len(traces(b)),providerOccurrenceEqual=True,resultHashes=hashes,finalOutput=stats)
 records.append(record)
 print(name,'files',len(files),'frames',len(times),'time',times[-1],'all bytes/trace/providers/statistics equal')
(base/'evidence/comparison-final.json').write_text(json.dumps(records,indent=2))
(base/'evidence/final-binary.json').write_text(json.dumps(dict(commit=__import__('subprocess').check_output(['git','rev-parse','HEAD'],text=True).strip(),before=hashlib.sha256((base/'baseline/sonicSolver').read_bytes()).hexdigest(),after=hashlib.sha256((root/'build/sonicSolver').read_bytes()).hexdigest()),indent=2))
