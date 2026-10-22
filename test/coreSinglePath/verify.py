"""Compare recorded A4 payloads, traces, physical-volume statistics and MPI matrices."""
from pathlib import Path
import hashlib,json,re,math,sys,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1]
sys.path.insert(0,str(root/'tools'))
from check_distributed_pressure import fields,matrices
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
def payloads(case):
 return {str(f.relative_to(case/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (case/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix!='.pvd'}
records=[]
for name in ['PISO','SIMPLE','PIMPLE','conservativePressure','sod','stageSST','eulerianSST']:
 a,b=base/'before'/name,base/'after'/name
 assert payloads(a)==payloads(b),name
 trace=lambda c:[s for s in (c/'run.log').read_text().splitlines() if s.startswith(('[SF PLAN]','[SF] Step time:'))]
 assert trace(a)==trace(b),name
 last=sorted(x for x in (b/'result').glob('*.vts') if not x.name.startswith('._'))[-1];stats=data(last);assert stats==data(a/'result'/last.name)
 times=[v.attrib['timestep'] for v in E.parse(next((b/'result').glob('*.pvd'))).findall('.//DataSet')]
 records.append({'case':name,'files':len(payloads(b)),'traceLines':len(trace(b)),'frames':len(times),'finalTime':times[-1],'statistics':stats})
(base/'evidence/statistics-final.json').write_text(json.dumps(records,indent=2))
a,b=base/'before/PISO-MPI2',base/'after/PISO-MPI2'
assert payloads(a)==payloads(b)
cross={}
for label in ['before','after']:
 left,right=fields(base/label/'PISO',2),fields(base/label/'PISO-MPI2',2);assert left.keys()==right.keys()
 u=max(abs(left[k][c]-right[k][c]) for k in left for c in range(3))
 span=max(max(right[k][3]-left[k][3] for k in left if k[2]==z)-min(right[k][3]-left[k][3] for k in left if k[2]==z) for z in {k[2] for k in left})
 assert matrices(base/label/'PISO',0)==matrices(base/label/'PISO-MPI2',0)
 cross[label]={'velocityMaxDifference':u,'pressureShapeDifference':span,'matrixRhsMaxDifference':0,'originalVelocityTolerance':1e-8,'originalPressureTolerance':1e-6,'velocityPass':u<1e-8,'pressurePass':span<1e-6}
assert cross['before']==cross['after']
clock=lambda f:re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',f.read_text())
assert clock(base/'evidence/mpi-before-host.log')==clock(base/'evidence/mpi-after-host-final.log')==clock(base/'after/PISO/run.log')
summary={'sameConfigurationFilesEqual':len(payloads(b)),'differences':[],'serialVsMpi':cross,'ranks':2,'steps':2,'clock':clock(base/'after/PISO/run.log')}
(base/'evidence/mpi-comparison.json').write_text(json.dumps(summary,indent=2))
print('Seven serial pairs and same-configuration MPI payloads/clock/matrices match.')
print('Serial/MPI pressure shape: existing tolerance pass =',cross['after']['pressurePass'])
