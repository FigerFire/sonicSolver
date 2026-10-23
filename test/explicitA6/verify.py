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
