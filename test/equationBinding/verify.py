"""Verify immutable before/after payloads, traces, clocks and physical statistics."""
from pathlib import Path
import json,hashlib,ast,xml.etree.ElementTree as E,math
base=Path(__file__).resolve().parent;root=base.parents[1]
source=ast.parse((root/'test/coreSinglePath/verify.py').read_text())
for n in source.body:
 if isinstance(n,ast.FunctionDef) and n.name=='data':exec(compile(ast.Module(body=[n],type_ignores=[]),'physical-quadrature','exec'))
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix not in ['.pvd','.txt']}
def trace(p):return [line for line in (p/'run.log').read_text().splitlines() if line.startswith(('[SF PLAN]','[SF] Step time:','[SF] Scalar step:','[SF SCALAR stage]'))]
records=[]
for group,names in [('', ['absent','laminar','gravity','kepsilon','kepsilonGravity','disabled','override','sstOnce']),('protection',['flowPair','scalarPair','stageSST','eulerianSST'])]:
 for name in names:
  a,b=base/group/'before'/name,base/group/'after'/name
  assert payload(a)==payload(b),(group,name,'numeric payload differs')
  assert trace(a)==trace(b),(group,name,'trace differs')
  frames=E.parse(next((b/'result').glob('*.pvd'))).findall('.//DataSet');last=b/'result'/frames[-1].get('file');stats=data(last)
  assert stats==data(a/'result'/last.relative_to(b/'result'))
  records.append(dict(group=group or 'module',case=name,frames=len(frames),finalTime=frames[-1].get('timestep'),statistics=stats,traceLines=len(trace(b)),payloadFiles=len(payload(b))))
(base/'evidence/regression.json').write_text(json.dumps(records,indent=2))
assert hashlib.sha256((root/'build/sonicSolver').read_bytes()).digest()==hashlib.sha256((base/'candidate/sonicSolver').read_bytes()).digest()
print('12 before/after numeric payloads, Plan trace, clock and physical-volume statistics identical')
def a7payload(folder):
 return {str(p.relative_to(folder)):hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.rglob('*.vts') if not p.name.startswith('._')}
a7before,a7after=base/'beforeA7/cases',base/'afterA7/cases'
if a7before.exists() and a7after.exists():
 old,new=a7payload(a7before),a7payload(a7after)
 assert old==new and len(old)==39,'A7 frozen numerical payload differs'
 count=len({Path(name).parts[0] for name in old})
 assert count==20,count
 record=dict(cases=count,files=len(old),exact=True,
  beforeTestBinary=hashlib.sha256((base/'baseline/test_coefficientA7').read_bytes()).hexdigest(),
  afterTestBinary=hashlib.sha256((root/'test/coefficientA7/bin/test_coefficientA7').read_bytes()).hexdigest())
 (base/'evidence/a7-regression.json').write_text(json.dumps(record,indent=2))
 print('20 A7 cases / 39 frozen VTS payloads identical')
