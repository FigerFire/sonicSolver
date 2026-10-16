from pathlib import Path
import shutil,json,subprocess,os
B=Path(__file__).resolve().parent
for name,scheme in [('thermalOne','forwardEuler'),('thermalFive','forwardEuler'),('thermalRkNegative','classicalRK4')]:
 p=B/name
 if not p.exists():shutil.copytree(B/'cliVerified',p,ignore=shutil.ignore_patterns('result','*.log','._*'))
 reg=p/'models/models.yaml';s=reg.read_text();s+='\nmultiPhase:\n  type: multiPhase\n  file: models/multiPhase.yaml\n' if '\nmultiPhase:' not in s else '';reg.write_text(s)
 (p/'models/multiPhase.yaml').write_text('''SonicFile:
  object: models
  type: multiPhase
enabled: true
phaseSystem: thermal
temperature:
  enabled: true
  fieldName: thermalTemperature
  defaultValue: 300
  transport: true
  diffusion: true
  diffusivity: 0.2
  bounded: false
''')
 f=p/'solvers/numerics.yaml';f.write_text(f.read_text().replace('forwardEuler',scheme))
 f=p/'solvers/runtime.yaml';s=f.read_text().replace('endTime: 0.01','endTime: 0.001').replace('endStep: 0','endStep: '+('1' if name=='thermalOne' else '5'));f.write_text(s)
 cmd=[str(B.parents[1]/'build/sonicSolver'),'explain',str(p)]
 r=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'TMPDIR':str(B/'evidence')});(p/'explain.log').write_text(r.stdout+r.stderr);print(name,'explain',r.returncode,flush=True)
 if name!='thermalRkNegative' and r.returncode==0:
  cmd[1]='run';r=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'TMPDIR':str(B/'evidence'),'SF_PLAN_TRACE':'1'});(p/'run.log').write_text(r.stdout+r.stderr);print(name,'run',r.returncode,flush=True)
