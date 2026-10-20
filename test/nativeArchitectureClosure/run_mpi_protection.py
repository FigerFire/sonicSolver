from pathlib import Path
import re,subprocess,os,shutil,json,hashlib
B=Path(__file__).resolve().parent;R=B.parents[1];out=[]
cases={'sod4':('test/levelSetMpiD1/protectionFinal/sod-after',4),
 'Peskin2':('test/IBM/cylinderFlowPeskinMPI2',2),'Peskin4':('test/IBM/cylinderFlowPeskinMPI4',4),
 'FTS2':('test/IBM/cylinderFlowVelocityForcingMPI2',2),'FTS4':('test/IBM/cylinderFlowVelocityForcingMPI4',4),
 'DLM2':('test/IBM/cylinderFlowFictitiousDomainMPI2',2),'DLM4':('test/IBM/cylinderFlowFictitiousDomainMPI4',4)}
for name,(source,n) in cases.items():
 for version,binary in [('before',B/'baseline/sonicSolver'),('after',R/'build/sonicSolver')]:
  p=B/'mpiProtection'/(name+'-'+version)
  if not p.exists():shutil.copytree(R/source,p,ignore=shutil.ignore_patterns('result','*.log','checkpoint','._*'))
  f=p/'mesh/mesh.yaml';s=f.read_text();s=re.sub(r'(?m)^generator:.*\n','',s);s=re.sub(r'(?m)^generation:\n(?:[ \t].*\n)*','',s);f.write_text(s)
  f=p/'solvers/runtime.yaml';s=f.read_text();s=re.sub(r'(?m)^endTime:.*','endTime: 0.001',s);s=re.sub(r'(?m)^endStep:.*','endStep: 0',s);s=re.sub(r'(?m)^writeControl:.*','writeControl: timeStep',s);s=re.sub(r'(?m)^writeInterval:.*','writeInterval: 1',s);s=s.replace('createMesh: true','createMesh: false')
  if name=='sod4':
   s=re.sub(r'(?ms)^parallel:\n(?:[ \t].*\n)*','parallel:\n  enabled: true\n  nProcs: 4\n  split: [4, 1, 1]\n',s)
  f.write_text(s);cmd=['mpirun','-np',str(n),str(binary),'run',str(p)]
  with (p/'run.log').open('w') as log:r=subprocess.run(cmd,cwd=p,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'TMPDIR':str(B/'evidence'),'SF_PLAN_TRACE':'1'},timeout=180)
  print(name,version,r.returncode,flush=True);out.append(dict(name=name,version=version,source=source,cmd=cmd,exit=r.returncode,binary=hashlib.sha256(binary.read_bytes()).hexdigest()))
 a=B/'mpiProtection'/(name+'-before');b=B/'mpiProtection'/(name+'-after')
 def data(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('.') and f.suffix!='.pvd'}
 out[-1]['numerical_output_equal']=bool(data(a)) and data(a)==data(b)
 print('numeric files exact',out[-1]['numerical_output_equal'],flush=True)
 (B/'evidence/mpi-protection-manifest.json').write_text(json.dumps(out,indent=2))
 if not all(v['exit']==0 for v in out[-2:]) or not out[-1]['numerical_output_equal']:raise RuntimeError(name)
