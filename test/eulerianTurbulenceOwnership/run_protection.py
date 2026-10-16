from pathlib import Path
import subprocess,shutil,re,os,json,hashlib
r=Path(__file__).resolve().parent;root=r.parents[1]
cases={'sod':'test/registryIdentity/sod-after','ghost':'test/registryIdentity/ghost-after','stageSST':'test/registryIdentity/stageSST-after','legacyGhostSST':'test/registryIdentity/authoredGhostSST-after','pressure':'test/pressure/constantDensityPiso','legacySST':'test/turbulence/singleFluidSST','Peskin':'test/IBM/cylinderFlowPeskin','forcing':'test/IBM/cylinderFlowVelocityForcing','fractionalDLM':'test/IBM/cylinderFlowFictitiousDomain','DFM':'test/IBM/cylinderFlowDFMExplicitSelfPropelled'}
def ignore_disappeared(function,path,error):
 if not isinstance(error[1],FileNotFoundError):raise error[1]
manifest=[]
for name,source in cases.items():
 for stage,binary in [('before',r/'evidence/sonicSolver-before'),('after',root/'build/sonicSolver')]:
  p=r/'protection'/(name+'-'+stage);shutil.copytree(root/source,p,ignore=shutil.ignore_patterns('result','._*','checkpoints','*.log','baseline.json'),dirs_exist_ok=True)
  if (p/'result').exists():shutil.rmtree(p/'result',onerror=ignore_disappeared)
  f=p/'solvers/runtime.yaml';s=f.read_text();s=re.sub(r'(?m)^endTime:.*','endTime: 0.0002' if name=='pressure' else 'endTime: 0.001',s);s=re.sub(r'(?m)^endStep:.*','endStep: 0',s);s=re.sub(r'(?m)^writeControl:.*','writeControl: timeStep',s);s=re.sub(r'(?m)^writeInterval:.*','writeInterval: 1',s);f.write_text(s)
  if name in ('pressure','legacySST'):
   f=p/'models/thermoDynamics.yaml';s=f.read_text()
   if name=='pressure':s+='  transport:\n    mu: 0.01\n    Pr: 0.72\n'
   else:s=s.replace('properties: {}','properties:\n  equationOfState:\n    gamma: 1.4\n    R: 287.05\n  transport:\n    mu: 0.00001\n    Pr: 0.72')
   f.write_text(s)
  cmd=[str(binary),'run',str(p)]
  if 'createMesh: true' in (p/'solvers/runtime.yaml').read_text():
   mesh=subprocess.run(cmd,text=True,capture_output=True);(p/'mesh-generation.log').write_text(mesh.stdout+mesh.stderr)
   if mesh.returncode:raise RuntimeError('Mesh generation failed: '+name)
   f=p/'solvers/runtime.yaml';f.write_text(f.read_text().replace('createMesh: true','createMesh: false'))
  run=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'SF_PLAN_TRACE':'1'});(p/'run.log').write_text(run.stdout+run.stderr)
  manifest.append(dict(case=name,stage=stage,source=source,command=cmd,exit=run.returncode,binary=hashlib.sha256(binary.read_bytes()).hexdigest()))
  print(name,stage,run.returncode,flush=True)
 if all(v['exit']==0 for v in manifest[-2:]):
  a=r/'protection'/(name+'-before')/'result';b=r/'protection'/(name+'-after')/'result';files=[f for f in a.rglob('*.vts') if not f.name.startswith('.')];other=[f for f in b.rglob('*.vts') if not f.name.startswith('.')]
  same=bool(files) and len(files)==len(other) and all(f.read_bytes()==(b/f.relative_to(a)).read_bytes() for f in files);manifest[-1]['vts_equal']=same;print('strict VTS',same,len(files),flush=True)
  if not same:raise RuntimeError(name)
  traces=lambda p:[line for line in (p/'run.log').read_text().splitlines() if line.startswith('[SF PLAN]')]
  trace_equal=traces(a.parent)==traces(b.parent);manifest[-1]['trace_equal']=trace_equal
  if not trace_equal:raise RuntimeError(name+' operation/time trace')
  outputs=[f for f in a.rglob('*') if f.is_file() and not f.name.startswith('.')];manifest[-1]['output_equal']=all((b/f.relative_to(a)).exists() and f.read_bytes()==(b/f.relative_to(a)).read_bytes() for f in outputs)
  if not manifest[-1]['output_equal']:raise RuntimeError(name+' output')
(r/'evidence/protection-manifest.json').write_text(json.dumps(manifest,indent=2))
