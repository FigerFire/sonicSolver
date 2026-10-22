from pathlib import Path
import subprocess,json,shutil,os,sys,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1]
os.chdir(base)
label=sys.argv[1];binary=base/'baseline/sonicSolver' if label=='before' else root/'build/sonicSolver'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
for name in ['sod','stageSST','eulerianSST']:
 p=base/label/name
 if not p.exists():shutil.copytree(root/'test/mathWhatA2/after'/name,p,ignore=shutil.ignore_patterns('result','*.log','numerical.txt','._*'))
 runtime=p/'solvers/runtime.yaml'
 if name in ['pairRK4','linearTermsRK4']:
  runtime.write_text((root/'test/mathWhatA2/after'/name/'solvers/runtime.yaml').read_text())
 else:
  import re
  text=runtime.read_text()
  text=re.sub(r'endStep:.*','endStep: 1',text) if 'endStep:' in text else text+'\nendStep: 1\n'
  runtime.write_text(text)
 for mode in ['explain','run']:
  r=subprocess.run([str(binary),mode,str(p)],cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(p/(mode+'.log')).write_text(r.stdout)
  print(label,name,mode,r.returncode,r.stdout[-550:] if r.returncode else '',flush=True)
  if r.returncode:break
