"""Copied IBM cases: full/wall preserves horizons; startup explicitly limits to one step."""
from pathlib import Path
import json,os,shutil,subprocess,sys,hashlib,re
from concurrent.futures import ThreadPoolExecutor,as_completed
base=Path(__file__).resolve().parent;root=base.parents[1]
os.chdir(base)
label=sys.argv[1];mode=sys.argv[2] if len(sys.argv)>2 else 'full';binary=base/'baseline/sonicSolver' if label=='before' else root/'build/sonicSolver'
cases=root/'test/IBM/coreSinglePath'/('startup/'+label if mode=='startup' else label);cases.mkdir(parents=True,exist_ok=True)
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
names=['cylinderFlowGhost','cylinderFlowSSTGhost','cylinderFlowPeskin','cylinderFlowVelocityForcing',
 'cylinderFlowVelocityForcingBP','cylinderFlowFictitiousDomain','cylinderFlowDFMExplicitSelfPropelled',
 'cylinderFlowDFMFractionalStepSelfPropelled','cylinderFlowFullyImplicitDLM',
 'cylinderFlowDFMImplicitSelfPropelled','cylinderFlowDFMAugmentedLagrangian']
if mode=='wall':names=names[:2]
manifest=base/'evidence'/('ibm-'+('startup-' if mode=='startup' else '')+label+'.json')
previous=json.loads(manifest.read_text()) if manifest.exists() else []
identity=hashlib.sha256(binary.read_bytes()).hexdigest()
records=[r for r in previous if Path(r['case']).name not in names or r['status']=='Completed' and r['binary']==identity]
completed={Path(r['case']).name for r in records}
def run(name):
 p=cases/name
 if not p.exists():shutil.copytree(root/'test/IBM'/name if label=='before' or mode=='startup' else root/'test/IBM/coreSinglePath/before'/name,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 if mode=='startup':
  f=p/'solvers/runtime.yaml';text=f.read_text()
  text=re.sub(r'endStep:.*','endStep: 1',text) if 'endStep:' in text else text+'\nendStep: 1\n'
  text=re.sub(r'writeControl:.*','writeControl: timeStep',text);text=re.sub(r'writeInterval:.*','writeInterval: 1',text);f.write_text(text)
 record={'scope':'one-step composition diagnostic' if mode=='startup' else 'configured full horizon','case':str(p),'source':str(root/'test/IBM'/name),'binary':hashlib.sha256(binary.read_bytes()).hexdigest(),'models':(p/'models/IBM.yaml').read_text(),'runtime':(p/'solvers/runtime.yaml').read_text()}
 for action in ['explain','run']:
  command=[str(binary),action,str(p)];record[action+'Command']=command
  with (p/(action+'.log')).open('w') as log:
   r=subprocess.Popen(command,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   while True:
    chunk=r.stdout.read(65536)
    if not chunk:break
    log.write(chunk.decode('utf-8',errors='replace'))
   r.wait()
  record[action+'Exit']=r.returncode
 log=(p/'run.log').read_text();record['lastClock']=re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',log)[-1:]
 record['status']='Completed' if record['runExit']==0 else 'FailedOrUnsupported'
 return record
with ThreadPoolExecutor(max_workers=2) as pool:
 for future in as_completed([pool.submit(run,name) for name in names if name not in completed]):
  record=future.result();records.append(record)
  manifest.write_text(json.dumps(sorted(records,key=lambda r:Path(r['case']).name),indent=2))
  print(label,Path(record['case']).name,record['status'],record['lastClock'],flush=True)
