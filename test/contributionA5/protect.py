from pathlib import Path
import os,shutil,subprocess,sys,json,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1];label=sys.argv[1];binary=base/'baseline/sonicSolver' if label=='before' else base/'candidate/sonicSolver';records=[]
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
for name,source in [('flowPair',root/'test/resourceA5/composition/flowPair'),('scalarPair',root/'test/resourceA5/composition/scalarPair'),('stageSST',root/'test/resourceA5/after/stageSST'),('eulerianSST',root/'test/resourceA5/after/eulerianSST'),('sstGhost',root/'test/IBM/cylinderFlowSSTGhost')]:
 p=base/'protection'/label/name
 if not p.exists():shutil.copytree(source if label=='before' else base/'protection/before'/name,p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store','numerical.txt'))
 identity=hashlib.sha256(binary.read_bytes()).hexdigest()
 cmd=[str(binary),'run',str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_bytes(r.stdout)
 records.append(dict(name=name,path=str(p),binary=identity,command=cmd,exit=r.returncode));print(label,name,r.returncode,flush=True)
 if r.returncode:raise RuntimeError(r.stdout[-2000:].decode())
(base/'evidence'/('protection-'+label+'.json')).write_text(json.dumps(records,indent=2))
