"""A4 frozen pressure output, matrix and operation-order comparison."""
from pathlib import Path
import hashlib, json, os, re, shutil, subprocess, sys
base=Path(__file__).resolve().parent
root=base.parents[1]
os.chdir(base)
sys.path.insert(0,str(root/'tools'))
from check_distributed_pressure import prepare, invoke
label=sys.argv[1]
binary=base/'baseline/sonicSolver' if label=='before' else root/'build/sonicSolver'
os.environ.update(DEVELOPER_DIR='/Library/Developer/CommandLineTools',TMPDIR=str(base/'evidence'),SF_PLAN_TRACE='1')
out=base/label;out.mkdir(exist_ok=True)
for name in ['PISO','SIMPLE','PIMPLE','conservativePressure']:
    p=out/name
    if not p.exists():
        if label=='before' and name!='conservativePressure':
            p=prepare(binary,root/'test/pressure/constantDensityPiso',out,name,(1,1,1),name,2)
        else:
            source=base/'before'/name if label=='after' else root/'test/derivedStateA3/after/conservativePressure'
            shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','._*'))
    env={**os.environ,'SF_PRESSURE_TRACE_DIR':str(p/'result/assembly')}
    with (p/'explain.log').open('w') as log:
        subprocess.run([str(binary),'explain',str(p)],cwd=p,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    with (p/'run.log').open('w') as log:
        subprocess.run([str(binary),'run',str(p)],cwd=p,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(label,name,'passed',flush=True)
if label=='after':
    comparison={}
    for name in ['PISO','SIMPLE','PIMPLE','conservativePressure']:
        a,b=base/'before'/name,base/'after'/name
        files=lambda p:{str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._')}
        af,bf=files(a),files(b)
        # PVD manifests contain paths; VTS, VTM, assembly and diagnostic payloads do not.
        af={k:v for k,v in af.items() if not k.endswith('.pvd')};bf={k:v for k,v in bf.items() if not k.endswith('.pvd')}
        trace=lambda p:[line for line in (p/'run.log').read_text().splitlines() if '[SF PLAN]' in line or 'Step time:' in line]
        difference=[k for k in af.keys()|bf.keys() if af.get(k)!=bf.get(k)]
        comparison[name]={'files':len(af),'differences':difference,'clock_and_operations_identical':trace(a)==trace(b)}
        assert not difference,(name,difference)
        assert trace(a)==trace(b),name
    (base/'evidence/comparison.json').write_text(json.dumps(comparison,indent=2))
    print(json.dumps(comparison,indent=2))
