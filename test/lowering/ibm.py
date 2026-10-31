"""Full configured serial IBM protections, isolated native copies; no method substitutions."""
from pathlib import Path
import json,os,shutil,subprocess,hashlib,sys,re
base=Path(__file__).resolve().parent;root=base.parents[1];label=sys.argv[1];mode=sys.argv[2] if len(sys.argv)>2 else "full"
binary=base/('baseline' if label=='before' else 'candidate')/'sonicSolver'
names=['cylinderFlowViscousGhost','cylinderFlowPeskin','cylinderFlowVelocityForcing','cylinderFlowVelocityForcingBP','cylinderFlowFictitiousDomain','cylinderFlowDFMExplicitSelfPropelled','cylinderFlowDFMFractionalStepSelfPropelled','cylinderFlowFullyImplicitDLM','cylinderFlowDFMImplicitSelfPropelled','cylinderFlowDFMAugmentedLagrangian']
if mode=='full':names=names[:1]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'COPYFILE_DISABLE':'1'};records=[]
for name in names:
 p=base/('ibmStart' if mode=='start' else 'ibm')/label/name
 if not p.exists():shutil.copytree(root/'test/IBM'/name,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 if mode=='start':
  runtime=p/'solvers/runtime.yaml';text=runtime.read_text()
  text=re.sub(r'^endStep:.*\n','',text,flags=re.M)+"\nendStep: 2\n"
  text=re.sub(r'^writeControl:.*$', 'writeControl: timeStep',text,flags=re.M)
  text=re.sub(r'^writeInterval:.*$', 'writeInterval: 1',text,flags=re.M);runtime.write_text(text)
 command=[str(binary),'run',str(p)]
 with (p/'run.log').open('w') as log:r=subprocess.run(command,cwd=p,env=env,stdout=log,stderr=subprocess.STDOUT)
 records.append(dict(name=name,command=command,exit=r.returncode,binary=hashlib.sha256(binary.read_bytes()).hexdigest()))
 (base/'evidence'/('ibm-'+mode+'-'+label+'.json')).write_text(json.dumps(records,indent=2));print(label,name,r.returncode,flush=True)
