from pathlib import Path
import os,shutil,subprocess,json,hashlib,re,sys
base=Path(__file__).resolve().parent;root=base.parents[1];label=sys.argv[1]
binary=base/'baseline/sonicSolver' if label=='before' else base/'candidate/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'COPYFILE_DISABLE':'1'}
records=[]
for tight in [False]:
 for ranks in [1,2]:
  name=('tight' if tight else 'original')+('-MPI2' if ranks==2 else '-serial');p=base/'pressure'/label/name
  source=root/'test/coreSinglePath/before'/('PISO-MPI2' if ranks==2 else 'PISO')
  if not p.exists():shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
  if tight:
   f=p/'solvers/algorithm.yaml';s=f.read_text();s=re.sub(r'(absoluteTolerance:\s*)[^\s]+',r'\g<1>1e-16',s);f.write_text(s)
  cmd=(['mpirun','-np',str(ranks)] if ranks>1 else [])+[str(binary),'run',str(p)]
  e={**env,'SF_PRESSURE_TRACE_DIR':str(p/'result/assembly')}
  r=subprocess.run(cmd,cwd=p,env=e,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_bytes(r.stdout)
  records.append(dict(name=name,exit=r.returncode,command=cmd,binary=hashlib.sha256(binary.read_bytes()).hexdigest()))
  print(label,name,r.returncode,flush=True)
(base/'evidence'/('pressure-'+label+'.json')).write_text(json.dumps(records,indent=2))
raise SystemExit(any(x['exit'] for x in records))
