from pathlib import Path
import json,math,hashlib,re,ast,sys,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1];sys.path.insert(0,str(root/'tools'))
from check_distributed_pressure import fields,matrices
source=ast.parse((root/'test/resourceA5/verify.py').read_text())
# Physical quadrature source remains in the original helper.
source=ast.parse((root/'test/coreSinglePath/verify.py').read_text())
for n in source.body:
 if isinstance(n,ast.FunctionDef) and n.name=='data':exec(compile(ast.Module(body=[n],type_ignores=[]),'physical-quadrature','exec'))
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix not in ['.pvd','.txt']}
def trace(p):return [s for s in (p/'run.log').read_text().splitlines() if s.startswith(('[SF PLAN]','[SF] Step time:','[SF] Scalar step:'))]
def final(p):
 datasets=E.parse(next((p/'result').glob('*.pvd'))).findall('.//DataSet');return p/'result'/datasets[-1].get('file'),datasets
records=[]
for group,names in [('',['absent','laminar','gravity','kepsilon','kepsilonGravity','disabled','override','sstOnce']),('protection',['flowPair','scalarPair','stageSST','eulerianSST','sstGhost'])]:
 for name in names:
  a=base/group/'before'/name;b=base/group/'after'/name
  assert payload(a)==payload(b),(group,name,'payload difference')
  assert trace(a)==trace(b),(group,name,'execution difference')
  vts,frames=final(b);stats=data(vts);assert stats==data(final(a)[0])
  text=(b/'run.log').read_text();active=name in ['kepsilon','kepsilonGravity','override','sstOnce']
  if not group:
   binds=[line for line in text.splitlines() if line.startswith('[SF] Bound physical resource:')]
   assert bool(binds)==active,(name,binds)
   if not active:assert 'op=turbulence.' not in text and 'flow.turbulence' not in text
   else:assert 'owner=model.turbulence target=k access=write' in (b/'explain.log').read_text(),name
  records.append(dict(group=group or 'modules',case=name,frames=len(frames),finalTime=frames[-1].get('timestep'),payloadFiles=len(payload(b)),traceLines=len(trace(b)),statistics=stats))
# Same physical input with absent/laminar/disabled model yields identical data.
assert payload(base/'after/absent')==payload(base/'after/laminar')==payload(base/'after/disabled')
def centre(p):
 tree=E.parse(final(p)[0]);xyz=list(map(float,tree.find('.//Points/DataArray').text.split()));u=list(map(float,tree.find('.//PointData/DataArray[@Name="U"]').text.split()));at=min(range(len(xyz)//3),key=lambda i:(xyz[3*i]-50)**2+(xyz[3*i+1]-10)**2)
 return u[3*at]
uoff=centre(base/'after/laminar');uon=centre(base/'after/gravity');assert abs(uon-uoff-.01*.001)<2e-12
pressure={}
for label in ['before','after']:
 a=base/'pressure'/label/'original-serial';b=base/'pressure'/label/'original-MPI2'
 left,right=fields(a,2),fields(b,2);assert matrices(a,0)==matrices(b,0)
 span=max(max(right[k][3]-left[k][3] for k in left if k[2]==z)-min(right[k][3]-left[k][3] for k in left if k[2]==z) for z in {k[2] for k in left})
 pressure[label]=dict(shape=span,originalThreshold=1e-6,passOriginalThreshold=span<=1e-6,velocity=max(abs(left[k][c]-right[k][c]) for k in left for c in range(3)))
assert pressure['before']==pressure['after']
for name in ['original-serial','original-MPI2']:assert payload(base/'pressure/before'/name)==payload(base/'pressure/after'/name)
for name in ['serial','mpi2x']:assert payload(base/'levelSetMpi/before'/name)==payload(base/'levelSetMpi/after'/name)
report=dict(pairs=records,gravityBudget=dict(off=uoff,on=uon,expectedIncrement=.00001),pressure=pressure,historicalLevelSetReference='unchanged, unresolved; no replacement of frozen-kernel.json')
(base/'evidence/verification.json').write_text(json.dumps(report,indent=2));print('13 same-configuration pairs and physical statistics identical; gravity budget verified; original pressure failure preserved.')
