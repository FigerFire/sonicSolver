from pathlib import Path
import json,hashlib,re,xml.etree.ElementTree as E,ast,math
base=Path(__file__).resolve().parent;root=base.parents[1]
source=ast.parse((root/'test/howA9/verify.py').read_text())
for n in source.body:
 if isinstance(n,ast.FunctionDef) and n.name=='data':exec(compile(ast.Module(body=[n],type_ignores=[]),'physical-quadrature','exec'))
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix not in ('.pvd','.txt')}
def clock(p):return re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',(p/'run.log').read_text())
records=[]
for left in json.loads((base/'evidence/ibm-start-before.json').read_text()):
 name=left['name'];right=next(r for r in json.loads((base/'evidence/ibm-start-after.json').read_text()) if r['name']==name)
 a,b=base/'ibmStart/before'/name,base/'ibmStart/after'/name;assert left['exit']==right['exit'],name
 if left['exit']==0:
  assert payload(a)==payload(b) and clock(a)==clock(b),(name,'first difference')
  assert len(clock(b))==2,name
  records.append(dict(case=name,status='two-step protection only',clock=clock(b),sameConfigurationExact=True,statistics=data(sorted(f for f in (b/'result').glob('*.vts') if not f.name.startswith('._'))[-1])))
 else:
  reason=lambda p:next(s for s in (p/'run.log').read_text().splitlines() if 'Fatal:' in s).split('Fatal:',1)[1].strip()
  assert reason(a)==reason(b),name
  records.append(dict(case=name,status='configuration blocked before and after',reason=reason(b)))
a,b=base/'ibm/before/cylinderFlowViscousGhost',base/'ibm/after/cylinderFlowViscousGhost'
assert payload(a)==payload(b) and clock(a)==clock(b),'Full no-slip Ghost first difference'
assert len(clock(b))==2000 and float(clock(b)[-1][0])==.02
records.append(dict(case='cylinderFlowViscousGhost',status='configured full .02 seconds',clockSteps=len(clock(b)),finalTime=clock(b)[-1][0],payloadFiles=len(payload(b)),sameConfigurationExact=True,statistics=data(sorted(f for f in (b/'result').glob('*.vts') if not f.name.startswith('._'))[-1])))
(base/'evidence/ibm-comparison.json').write_text(json.dumps(records,indent=2)+'\n')
print('Full no-slip Ghost .02/2000 steps, 7 method-specific two-step protections exact; 3 original configuration conflicts unchanged')
