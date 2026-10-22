"""Compare actual configured-wall runs separately from one-step composition diagnostics."""
from pathlib import Path
import hashlib,json,re,sys,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1]
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix!='.pvd'}
def trace(p):return [l for l in (p/'run.log').read_text().splitlines() if l.startswith(('[SF PLAN]','[SF] Step time:','[SF] State closure:'))]
summary=[]
for mode in ['startup','full']:
 records={label:json.loads((base/'evidence'/('ibm-'+('startup-' if mode=='startup' else '')+label+'.json')).read_text()) for label in ['before','after']}
 before={Path(r['case']).name:r for r in records['before']}
 for r in records['after']:
  name=Path(r['case']).name
  if mode=='full' and name not in {'cylinderFlowGhost','cylinderFlowSSTGhost'}:continue
  old=before[name];a,b=Path(old['case']),Path(r['case'])
  row={'scope':mode,'case':name,'beforeBinary':old['binary'],'afterBinary':r['binary'],'beforeExit':old['runExit'],'afterExit':r['runExit']}
  assert old['runExit']==r['runExit'],name
  if r['runExit']==0:
   assert payload(a)==payload(b),name
   assert trace(a)==trace(b),name
   row.update(status='Matched',payloadFiles=len(payload(b)),traceLines=len(trace(b)),lastClock=r['lastClock'])
   pvd=[f for f in (b/'result').glob('*.pvd') if not f.name.startswith('._')];row['outputTimes']=[x.attrib['timestep'] for f in pvd for x in E.parse(f).findall('.//DataSet')]
   if mode=='full':assert float(row['outputTimes'][0])==0 and float(row['outputTimes'][-1])==(1 if name=='cylinderFlowGhost' else .02),name
  else:
   reason=lambda p:re.findall(r'\[SF\] Fatal: (.*)',(p/'run.log').read_text())
   assert reason(a)==reason(b),name
   row.update(status='ConfigurationInvalid' if 'Conflicting HOW' in (b/'run.log').read_text() else 'Rejected',reason=reason(b))
  summary.append(row)
(base/'evidence/ibm-comparison.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
