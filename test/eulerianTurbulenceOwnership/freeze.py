from pathlib import Path
import re,shutil,subprocess,os,json
r=Path(__file__).resolve().parent
if (r/'evidence/sonicSolver-instrumented').exists():raise RuntimeError('Frozen executable already exists; use compare.py for regression, never overwrite the baseline.')
for model in ['Smagorinsky','laminar']:
 p=r/'baseline'/(model+'-all')
 for name in ['internalField','boundaries']:
  f=p/'fields'/(name+'.yaml');f.write_text(re.sub(r'(?ms)^(k|omega):\n.*?(?=^[A-Za-z]|\Z)','',f.read_text()))
 run=subprocess.run([str(r/'evidence/sonicSolver-before'),'run',str(p)],text=True,capture_output=True,env={**os.environ,'SF_PLAN_TRACE':'1'});(p/'run.log').write_text(run.stdout+run.stderr);print(model,run.returncode)
shutil.copy2('build/sonicSolver',r/'evidence/sonicSolver-instrumented')
manifest=[]
for src in sorted((r/'baseline').iterdir()):
 if not src.is_dir() or src.name.startswith('.'):continue
 p=r/'instrumented'/src.name;shutil.copytree(src,p,ignore=shutil.ignore_patterns('result','*.log','._*'),dirs_exist_ok=True)
 cmd=[str(r/'evidence/sonicSolver-instrumented'),'run',str(p)]
 run=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'SF_PLAN_TRACE':'1','SF_EULERIAN_TURBULENCE_EVIDENCE':str(p/'numerical.txt')});(p/'run.log').write_text(run.stdout+run.stderr)
 files=[f for f in (src/'result').rglob('*.vts') if not f.name.startswith('.')];same=all(f.read_bytes()==(p/'result'/f.relative_to(src/'result')).read_bytes() for f in files)
 print(src.name,run.returncode,'instrumentation identical',same,flush=True)
 manifest.append(dict(case=src.name,command=cmd,exit=run.returncode,instrumentation_equal=same))
 if run.returncode or not same:raise RuntimeError(src.name+' '+run.stderr)
(r/'evidence/instrumentation-manifest.json').write_text(json.dumps(manifest,indent=2))
