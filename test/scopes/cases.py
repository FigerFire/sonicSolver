"""Same-input lowering regression; frozen inputs and binaries stay below test/."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
from convert import convert

base=Path(__file__).resolve().parent;root=base.parents[1]
phase=sys.argv[1];binary=base/('baseline' if phase=='before' else 'candidate')/'sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
names=[item['case'] for item in json.loads((root/'test/lowering/evidence/cases-after.json').read_text())]
sources={name:root/'test/lowering/after'/name for name in names}
records=[]
for name,source in sources.items():
    case=base/phase/name
    if not case.exists():
        original=source if phase=='before' else base/'before'/name
        shutil.copytree(original,case,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
    if phase=='after':convert(case)
    command=[str(binary),'run',str(case)]
    run=subprocess.run(command,cwd=case,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (case/'run.log').write_bytes(run.stdout)
    record=dict(case=name,source=str(source),path=str(case),command=command,binary=hashlib.sha256(binary.read_bytes()).hexdigest(),exit=run.returncode)
    records.append(record);print(phase,name,run.returncode,flush=True)
    if run.returncode:raise RuntimeError(run.stdout[-4000:].decode())
(base/'evidence'/('cases-'+phase+'.json')).write_text(json.dumps(records,indent=2)+'\n')
