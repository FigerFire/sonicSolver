from pathlib import Path
import os,shutil,json,re,subprocess
base=Path(__file__).resolve().parent;binary=base/'candidate/sonicSolver';env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')};records=[]
for name in ['stageKepsilon','initial','whichConflict','unknownMethod','modelConflict']:
 p=base/'negative'/name
 if not p.exists():shutil.copytree(base/'before/override',p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 shutil.copytree(base/'before/override',p,dirs_exist_ok=True,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 if name=='stageKepsilon':
  f=p/'models/turbulence.yaml';d=json.loads(f.read_text());d['transportMethod']='stageConservative';f.write_text(json.dumps(d))
 if name=='initial':
  f=p/'fields/internalField.yaml';s=f.read_text();s=re.sub(r'^epsilon:\n(?:  [^\n]*\n)*','',s,flags=re.M);f.write_text(s)
 if name in ['whichConflict','unknownMethod']:
  f=p/'solvers/which.yaml';d=json.loads(f.read_text())
  if name=='whichConflict':d['bindings'].append(d['bindings'][0])
  else:d['bindings'][0]['method']='MissingNumericalImplementation'
  f.write_text(json.dumps(d))
 if name=='modelConflict':
  f=p/'models/models.yaml';f.write_text(f.read_text()+'\notherTurbulence:\n  type: turbulence\n  file: models/otherTurbulence.yaml\n')
  d=json.loads((p/'models/turbulence.yaml').read_text());d['SonicFile']['type']='otherTurbulence';(p/'models/otherTurbulence.yaml').write_text(json.dumps(d))
 cmd=[str(binary),'explain',str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(p/'explain.log').write_text(r.stdout)
 assert r.returncode and not list((p/'result').glob('*.vts')),(name,r.stdout)
 reason=r.stdout.split('Reason:')[-1].strip()
 expected={'stageKepsilon':'single-fluid SST only','initial':'initialization missing: epsilon','whichConflict':'Conflicting explicit WHICH','unknownMethod':'no registered implementation','modelConflict':'Conflicting explicit target Turbulence'}[name]
 assert expected in reason,(name,reason)
 records.append(dict(case=name,command=cmd,exit=r.returncode,reason=reason));print(name,reason,flush=True)
(base/'evidence/negative.json').write_text(json.dumps(records,indent=2))
