"""One-step SSPRK3/TENO5 shared-flux MPI-4 protection; no formal case changes."""
from pathlib import Path
import subprocess,shutil,os,json,hashlib,re
base=Path(__file__).resolve().parent;records=[]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')}
for label,folder in [('before','baseline'),('after','candidate')]:
 p=base/'mpi'/label/'sod4';p.parent.mkdir(parents=True,exist_ok=True)
 if not p.exists():shutil.copytree(base/'before/sod',p,ignore=shutil.ignore_patterns('result','*.log','._*'))
 f=p/'solvers/runtime.yaml';s=f.read_text();s=s.replace('enabled: false','enabled: true').replace('nProcs: 1','nProcs: 4').replace('split: [1, 1, 1]','split: [4, 1, 1]');f.write_text(s)
 binary=base/folder/'sonicSolver';cmd=['mpirun','-np','4',str(binary),'run',str(p)]
 r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_bytes(r.stdout)
 assert not r.returncode,r.stdout[-5000:].decode();assert '4 rank(s)' in r.stdout.decode(),r.stdout[-5000:].decode()
 records.append(dict(phase=label,path=str(p),command=cmd,binary=hashlib.sha256(binary.read_bytes()).hexdigest(),exit=r.returncode))
 print(label,r.returncode,flush=True)
def payload(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._') and f.suffix!='.pvd'}
a,b=base/'mpi/before/sod4',base/'mpi/after/sod4';assert payload(a)==payload(b),'MPI4 first differing numerical payload'
clock=lambda p:re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',(p/'run.log').read_text())
assert clock(a)==clock(b);(base/'evidence/mpi.json').write_text(json.dumps(dict(cases=records,payloadFiles=len(payload(b)),sameConfigurationExact=True,clock=clock(b)),indent=2)+'\n')
print('MPI-4 frozen/current Sod outputs and clocks exact')
