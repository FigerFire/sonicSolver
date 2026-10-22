"""Host MPI protection with actual exit codes and logs, all files below test/."""
from pathlib import Path
import subprocess,os,json,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1];os.chdir(base)
env={**os.environ,'COPYFILE_DISABLE':'1','DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
records=[]
for name,command in [
 ('mpi-before-host',['mpirun','-np','2',str(base/'baseline/sonicSolver'),'run',str(base/'before/PISO-MPI2')]),
 ('mpi-after-host-final',['mpirun','-np','2',str(root/'build/sonicSolver'),'run',str(base/'after/PISO-MPI2')]),
 ('mpi-layout',['mpirun','-np','4',str(root/'build/test_distributedConstraintLayout')]),
 ('mpi-copy',['mpirun','-np','4',str(root/'build/test_sparseCanonicalCopy')]),
 ('mpi-kkt',['mpirun','-np','2',str(root/'build/test_distributedKKT')])]:
 runenv=dict(env)
 if 'host' in name:
  label='before' if name=='mpi-before-host' else 'after';runenv['SF_PRESSURE_TRACE_DIR']=str(base/label/'PISO-MPI2/result/assembly')
 r=subprocess.run(command,cwd=base,env=runenv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (base/'evidence'/(name+'.log')).write_bytes(r.stdout)
 records.append(dict(name=name,command=command,exit=r.returncode,binarySha256=hashlib.sha256(Path(command[3]).read_bytes()).hexdigest(),appleDoubleCleanupWarning=b'._0' in r.stdout))
 print(name,r.returncode,flush=True)
(base/'evidence/mpi-exits.json').write_text(json.dumps(records,indent=2))
raise SystemExit(any(r['exit'] for r in records))
