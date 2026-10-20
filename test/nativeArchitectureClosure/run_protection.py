from pathlib import Path
import subprocess,shutil,re,os,json,hashlib
B=Path(__file__).resolve().parent;R=B.parents[1]
cases={'sod':'test/levelSetMpiD1/protectionFinal/sod-after','ghost':'test/levelSetMpiD1/protectionFinal/ghost-after','stageSST':'test/levelSetMpiD1/protectionFinal/stageSST-after','ghostSST':'test/levelSetMpiD1/protectionFinal/legacyGhostSST-after','pressure':'test/levelSetMpiD1/protectionFinal/pressure-after','SST':'test/levelSetMpiD1/protectionFinal/legacySST-after','Peskin':'test/levelSetMpiD1/protectionFinal/Peskin-after','forcing':'test/levelSetMpiD1/protectionFinal/forcing-after','fractionalDLM':'test/levelSetMpiD1/protectionFinal/fractionalDLM-after','DFM':'test/levelSetMpiD1/protectionFinal/DFM-after'}
# Existing supported operators plus explicit negative production KKT; no added physics.
for name in ('cylinderFlowVelocityForcingBP','cylinderFlowDFMFractionalStepSelfPropelled','cylinderFlowFullyImplicitDLM','cylinderFlowDFMImplicitSelfPropelled','cylinderFlowDFMAugmentedLagrangian'):
 cases[name]='test/IBM/'+name
manifest=[]
for name,source in cases.items():
 for stage,binary in [('before',B/'baseline/sonicSolver'),('after',R/'build/sonicSolver')]:
  p=B/'protection'/(name+'-'+stage)
  if not p.exists():shutil.copytree(R/source,p,ignore=shutil.ignore_patterns('result','._*','checkpoints','checkpoint','*.log','baseline.json'))
  # Do not require old input reader to understand new mesh config. Both use the same frozen SFM.
  f=p/'mesh/mesh.yaml';s=f.read_text();s=re.sub(r'(?m)^generator:.*\n','',s);s=re.sub(r'(?m)^generation:\n(?:[ \t].*\n)*','',s);f.write_text(s)
  f=p/'solvers/runtime.yaml';s=f.read_text();s=re.sub(r'(?m)^endTime:.*','endTime: 0.0002' if name=='pressure' else 'endTime: 0.001',s);s=re.sub(r'(?m)^endStep:.*','endStep: 0',s);s=re.sub(r'(?m)^writeControl:.*','writeControl: timeStep',s);s=re.sub(r'(?m)^writeInterval:.*','writeInterval: 1',s);s=s.replace('createMesh: true','createMesh: false');f.write_text(s)
  cmd=[str(binary),'run',str(p)];run=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'SF_PLAN_TRACE':'1','TMPDIR':str(B/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools'},timeout=180);(p/'run.log').write_text(run.stdout+run.stderr)
  item=dict(case=name,stage=stage,source=source,command=cmd,exit=run.returncode,binary=hashlib.sha256(binary.read_bytes()).hexdigest());manifest.append(item)
  print(name,stage,run.returncode,flush=True)
 if all(v['exit']==0 for v in manifest[-2:]):
  a=B/'protection'/(name+'-before');b=B/'protection'/(name+'-after')
  files=lambda p:{str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('.') and f.suffix!='.pvd'}
  same=bool(files(a)) and files(a)==files(b);manifest[-1]['numerical_output_equal']=same
  traces=lambda p:[line for line in (p/'run.log').read_text().splitlines() if line.startswith('[SF PLAN]')]
  trace_equal=traces(a)==traces(b);manifest[-1]['trace_equal']=trace_equal
  print('exact numerical outputs / Plan trace',same,trace_equal,len(files(a)),flush=True)
  if not same or not trace_equal:raise RuntimeError(name+' first difference')
 else:manifest[-1]['expected_unsupported']=manifest[-2]['exit']==manifest[-1]['exit'] and all('Unsupported' in (B/'protection'/(name+'-'+stage)/'run.log').read_text() or 'unsupported' in (B/'protection'/(name+'-'+stage)/'run.log').read_text() for stage in ('before','after'))
 (B/'evidence/protection-manifest.json').write_text(json.dumps(manifest,indent=2))
