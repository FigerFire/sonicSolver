from pathlib import Path
import shutil, subprocess, json, os, hashlib
import re
r=Path('test/eulerianTurbulenceOwnership');r.mkdir(exist_ok=True);(r/'evidence').mkdir(exist_ok=True)
b=r/'evidence/sonicSolver-before'
if (r/'evidence/sonicSolver-before').exists():raise RuntimeError('Frozen executable already exists; use compare.py for regression, never overwrite the baseline.')
shutil.copy2('build/sonicSolver',b)
manifest=[]
for model in ['kOmegaSST','kEpsilon','Smagorinsky','laminar']:
 for variant in (['all','subset','reverse','outer3'] if model in ['kOmegaSST','kEpsilon'] else ['all']):
  name=model+'-'+variant;p=r/'baseline'/name
  src=Path('test/turbulence/euleriankEpsilon' if model=='kEpsilon' else 'test/turbulence/euleriankOmegaSST')
  shutil.copytree(src,p,ignore=shutil.ignore_patterns('result','._*','baseline.json','*.log'),dirs_exist_ok=True)
  t=p/'models/turbulence.yaml';v=t.read_text()
  if model=='laminar':v=v.replace('enabled: true','enabled: false')
  if model=='Smagorinsky':v=v.replace('simulationType: RAS','simulationType: LES').replace('RAS:','LES:').replace('model: kOmegaSST','model: Smagorinsky')
  if variant=='subset':v=v.replace('[water, air]','[air]')
  if variant=='reverse':v=v.replace('[water, air]','[air, water]')
  t.write_text(v)
  if variant=='outer3':
   t=p/'solvers/algorithm.yaml';v=t.read_text().replace('outerCorrectors: 1','outerCorrectors: 3').replace('pressureCorrectors: 1','pressureCorrectors: 2').replace('nonOrthogonalCorrectors: 0','nonOrthogonalCorrectors: 1');t.write_text(v)
  cmd=[str(b.resolve()),'run',str(p.resolve())]
  run=subprocess.run(cmd,capture_output=True,text=True,env={**os.environ,'SF_PLAN_TRACE':'1'})
  (p/'run.log').write_text(run.stdout+run.stderr)
  print(name,run.returncode,flush=True)
  manifest.append({'case':str(p),'source':str(src),'command':cmd,'exit':run.returncode,'binary':hashlib.sha256(b.read_bytes()).hexdigest()})
  if run.returncode:print(run.stdout+run.stderr);break
(r/'evidence/baseline-manifest.json').write_text(json.dumps(manifest,indent=2))
