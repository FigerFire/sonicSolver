from pathlib import Path
import ast,json,hashlib,math,re,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1]
source=ast.parse((root/'test/howA9/verify.py').read_text())
for n in source.body:
 if isinstance(n,ast.FunctionDef) and n.name=='data':exec(compile(ast.Module(body=[n],type_ignores=[]),'quadrature','exec'))
def frames(case):
 pvd=E.parse(next((case/'result').glob('*.pvd')));return [(e.get('timestep'),e.get('file')) for e in pvd.findall('.//DataSet')]
def payload(case):
 result={}
 for file in sorted((case/'result').rglob('*')):
  if not file.is_file() or file.name.startswith('._') or file.suffix in ('.pvd','.txt'):continue
  result[str(file.relative_to(case/'result'))]=hashlib.sha256(file.read_bytes()).hexdigest()
 return result
records=[]
for entry in json.loads((base/'evidence/cases-before.json').read_text()):
 name=entry['case'];left,right=base/'before'/name,base/'after'/name
 assert payload(left)==payload(right),(name,'first differing numerical payload')
 assert frames(left)==frames(right),(name,'physical output times')
 def trace(case):
  lines=[line for line in (case/'run.log').read_text().splitlines() if line.startswith(('[SF PLAN]','[SF] Step time:','[SF] Scalar step:','[SF SCALAR stage]'))]
  return [line.replace('stages/stage','recipe.stages/stage') for line in lines]
 a,b=trace(left),trace(right)
 if a!=b:
  differences=[(i,x,y) for i,(x,y) in enumerate(zip(a,b)) if x!=y]
  assert False,(name,'first trace difference',differences[:1],len(a),len(b))
 last=frames(right)[-1];stats=data(right/'result'/last[1]) if last[1].endswith('.vts') else None
 records.append(dict(case=name,finalTime=last[0],frames=len(frames(right)),payloadFiles=len(payload(right)),traceLines=len(b),statistics=stats))
start=json.loads((base/'evidence/start.json').read_text())
for file,digest in start['docs'].items():
 if file=='docs/problem20261020.md':continue # Current document explicitly updated by this task; all historical reports remain frozen.
 assert hashlib.sha256((root/file).read_bytes()).hexdigest()==digest,('Original user document changed',file)
(base/'evidence/regression.json').write_text(json.dumps(dict(cases=records,originalUserDocsPreserved=True),indent=2)+'\n')
print(len(records),'fresh baseline/candidate numerical payloads, output clocks and numerical operation traces are exact')
