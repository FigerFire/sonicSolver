from pathlib import Path
import os,shutil,subprocess,json,hashlib,sys
base=Path(__file__).resolve().parent;root=base.parents[1];label=sys.argv[1]
binary=base/'baseline/sonicSolver' if label=='before' else base/'candidate/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'COPYFILE_DISABLE':'1','SF_PLAN_TRACE':'1'}
records=[]
for name,ranks in [('serial',1),('mpi2x',2)]:
 p=base/'levelSetMpi'/label/name;source=root/'test/levelSetMpiD1/closeout/planeAdvection'/name
 if not p.exists():shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 cmd=(['mpirun','-np',str(ranks)] if ranks>1 else [])+[str(binary),'run',str(p)]
 r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_bytes(r.stdout)
 records.append(dict(name=name,exit=r.returncode,command=cmd,binary=hashlib.sha256(binary.read_bytes()).hexdigest()));print(label,name,r.returncode,flush=True)
if label=='after':
 for name,ranks in [('distributedConstraintLayout',4),('sparseCanonicalCopy',4),('distributedKKT',2)]:
  cmd=['mpirun','-np',str(ranks),str(root/'build'/('test_'+name))]
  r=subprocess.run(cmd,cwd=base,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(base/'evidence'/(name+'.log')).write_bytes(r.stdout)
  records.append(dict(name=name,exit=r.returncode,command=cmd,binary=hashlib.sha256(Path(cmd[3]).read_bytes()).hexdigest()));print(label,name,r.returncode,flush=True)
(base/'evidence'/('mpi-'+label+'.json')).write_text(json.dumps(records,indent=2))
raise SystemExit(any(x['exit'] for x in records))
