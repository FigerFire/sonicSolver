from pathlib import Path
import sys,ast,json,hashlib,re,math,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1]
sys.path.insert(0,str(root/'tools'))
from check_distributed_pressure import fields,matrices
source=ast.parse((base/'verify.py').read_text())
for n in source.body:
 if isinstance(n,ast.FunctionDef) and n.name=='data':exec(compile(ast.Module(body=[n],type_ignores=[]),'physical-quadrature','exec'))
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix!='.pvd'}
clock=lambda p:re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',(p/'run.log').read_text())
checks={}
for name in ['original-serial','original-MPI2']:
 a,b=base/'pressure/before'/name,base/'pressure/after'/name
 assert payload(a)==payload(b) and clock(a)==clock(b),name
 last=sorted(f for f in (b/'result').glob('*.vts') if not f.name.startswith('._'))[-1]
 checks[name]={'sameConfigurationExact':True,'payloadFiles':len(payload(b)),'clock':clock(b),'lastPatchStatistics':data(last)}
cross={}
for label in ['before','after']:
 serial,mpi=base/'pressure'/label/'original-serial',base/'pressure'/label/'original-MPI2'
 left,right=fields(serial,2),fields(mpi,2);assert left.keys()==right.keys()
 velocity=max(abs(left[k][c]-right[k][c]) for k in left for c in range(3))
 pressure=max(max(right[k][3]-left[k][3] for k in left if k[2]==z)-min(right[k][3]-left[k][3] for k in left if k[2]==z) for z in {k[2] for k in left})
 assert matrices(serial,0)==matrices(mpi,0)
 cross[label]={'velocityMaxDifference':velocity,'pressureShapeDifference':pressure,'matrixRhsMaxDifference':0,'originalVelocityTolerance':1e-8,'originalPressureTolerance':1e-6,'velocityPass':velocity<1e-8,'pressurePass':pressure<1e-6}
assert cross['before']==cross['after']
(base/'evidence/pressure-comparison.json').write_text(json.dumps({'sameConfiguration':checks,'serialVsMpi':cross},indent=2))
print('Original pressure serial and MPI-2 before/after payload, clock, first matrix/RHS exact. Cross-configuration:',cross['after'])
