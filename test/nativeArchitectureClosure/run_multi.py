from pathlib import Path
import os,shutil,subprocess,json,sys
B=Path(__file__).resolve().parent;R=B.parents[1];phase=sys.argv[1];binary=Path(sys.argv[2]).resolve();out=[]
for problem in ('planeAdvection','circleGeometry','csf'):
 for label,n in [('serial',1),('mpi2',2),('mpi4',4)]:
  p=B/phase/problem/label
  if not p.exists():shutil.copytree(R/'test/levelSetMpiD1/multiPatchCloseout'/problem/label,p,ignore=shutil.ignore_patterns('result','checkpoint','*.log','._*'))
  (p/'checkpoint').mkdir(exist_ok=True);cmd=(['mpirun','-np',str(n)] if n>1 else [])+[str(binary),'run',str(p)]
  with (p/'run.log').open('w') as f:
   run=subprocess.run(cmd,cwd=p,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'TMPDIR':str(B/'evidence'),'SF_PLAN_TRACE':'1','SF_LEVEL_SET_CHECKPOINT_DIR':str(p/'checkpoint')},timeout=180)
  print(problem,label,run.returncode,flush=True);out.append(dict(case=str(p),cmd=cmd,exit=run.returncode))
(B/'evidence'/('manifest-'+phase+'.json')).write_text(json.dumps(out,indent=2));assert all(r['exit']==0 for r in out)
