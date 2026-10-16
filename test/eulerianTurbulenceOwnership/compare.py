from pathlib import Path
import re,shutil,subprocess,os,json,sys,hashlib
r=Path(__file__).resolve().parent;stage=sys.argv[1];binary=Path('build/sonicSolver').resolve();manifest=[]
for src in sorted((r/'instrumented').iterdir()):
 if not src.is_dir() or src.name.startswith('.'):continue
 p=r/stage/src.name;shutil.copytree(src,p,ignore=shutil.ignore_patterns('result','*.log','._*','numerical.txt'),dirs_exist_ok=True)
 cmd=[str(binary),'run',str(p)]
 run=subprocess.run(cmd,text=True,capture_output=True,env={**os.environ,'SF_PLAN_TRACE':'1','SF_EULERIAN_TURBULENCE_EVIDENCE':str(p/'numerical.txt')});(p/'run.log').write_text(run.stdout+run.stderr)
 files=[f for f in (src/'result').rglob('*.vts') if not f.name.startswith('.')]
 same=run.returncode==0 and len(files)==len([f for f in (p/'result').rglob('*.vts') if not f.name.startswith('.')]) and all(f.read_bytes()==(p/'result'/f.relative_to(src/'result')).read_bytes() for f in files)
 matrix=not (src/'numerical.txt').exists() or (src/'numerical.txt').read_bytes()==(p/'numerical.txt').read_bytes()
 # Expanded local leaves are intentionally visible; original physical trace must remain exact.
 ops=lambda text:re.findall(r'\[PLAN\].*',text)
 physical=lambda text:[line for line in text.splitlines() if ('Time=' in line or 'time=' in line or 'turbulenceIter=' in line) and not line.startswith('[SF PLAN]')]
 a,b=physical((src/'run.log').read_text()),physical(run.stdout+run.stderr)
 print(stage,src.name,run.returncode,'VTS',same,'matrices',matrix,'dt/diagnostics',a==b,flush=True)
 manifest.append(dict(case=src.name,command=cmd,exit=run.returncode,vts_equal=same,matrix_equal=matrix,diagnostics_equal=a==b,binary=hashlib.sha256(binary.read_bytes()).hexdigest()))
 if not (same and matrix and a==b):
  (r/'evidence'/(stage+'-manifest.json')).write_text(json.dumps(manifest,indent=2));raise RuntimeError(src.name+'\n'+run.stderr[-1500:])
(r/'evidence'/(stage+'-manifest.json')).write_text(json.dumps(manifest,indent=2))
