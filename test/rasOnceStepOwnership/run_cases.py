from pathlib import Path
import sys,json,shutil,re,subprocess,os,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1];version=sys.argv[1]
label=sys.argv[2] if len(sys.argv)>2 else version
binary=base/'baseline/sonicSolver' if version=='before' else root/'build/sonicSolver'
specs={'sst':'test/nativeArchitectureClosure/protection/SST-after','kEpsilon':'test/nativeArchitectureClosure/protection/SST-after','les':'test/nativeArchitectureClosure/protection/SST-after','ghostSST':'test/nativeArchitectureClosure/protection/ghostSST-after','stageSST':'test/nativeArchitectureClosure/protection/stageSST-after','eulerianSST':'test/nativeArchitectureClosure/eulerianNative/kOmegaSST-outer3','eulerianKEpsilon':'test/nativeArchitectureClosure/eulerianNative/kEpsilon-all','sod':'test/nativeArchitectureClosure/protection/sod-after'}
records=[]
for name,source in specs.items():
 p=base/label/name
 if p.exists():raise RuntimeError('Refusing to overwrite '+str(p))
 shutil.copytree(root/source,p,ignore=shutil.ignore_patterns('result','checkpoint','*.log','._*','numerical.txt'))
 if name in ('sst','kEpsilon','les'):
  f=p/'fields/internalField.yaml';s=f.read_text();s=re.sub(r'(?m)^k:\n  uniform:.*','k:\n  default: 0.02\n  leftHalf: 0.015\n  rightHalf: 0.025',s)
  s=re.sub(r'(?m)^omega:\n  uniform:.*','omega:\n  default: 2\n  leftHalf: 1.5\n  rightHalf: 2.5',s)
  s=re.sub(r'(?m)^epsilon:\n  uniform:.*','epsilon:\n  default: 0.03\n  leftHalf: 0.025\n  rightHalf: 0.035',s);f.write_text(s)
 if name=='kEpsilon':
  f=p/'models/turbulence.yaml';f.write_text(f.read_text().replace('kOmegaSST','kEpsilon'))
  f=p/'fields/fields.yaml';s=f.read_text();start=s.index('epsilon:');end=s.index('\n',start)+1
  tail=re.search(r'(?m)^\S',s[end:]);stop=end+tail.start() if tail else len(s)
  s=s[:start]+s[start:stop].replace('output: false','output: true')+s[stop:];f.write_text(s)
 if name in ('kEpsilon','les'):
  f=p/'fields/fields.yaml';s=f.read_text()
  for field in (('omega',) if name=='kEpsilon' else ('k','omega','epsilon')):
   s=re.sub(r'(?ms)^'+field+r':.*?(?=^\S|\Z)',lambda m:m.group().replace('output: true','output: false'),s)
  f.write_text(s)
 if name=='les':
  f=p/'models/turbulence.yaml';f.write_text('SonicFile:\n  object: models\n  type: turbulence\nenabled: true\nsimulationType: LES\nLES:\n  model: Smagorinsky\n')
 cmd=[str(binary),'run',str(p)]
 env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','SF_EULERIAN_TURBULENCE_EVIDENCE':str(p/'numerical.txt')}
 run=subprocess.run(cmd,cwd=p,text=True,capture_output=True,env=env);text=run.stdout+run.stderr;(p/'run.log').write_text(text)
 records.append(dict(name=name,source=source,case=str(p),command=cmd,binary=hashlib.sha256(binary.read_bytes()).hexdigest(),exit=run.returncode))
 (base/'evidence'/('cases-'+label+'.json')).write_text(json.dumps(records,indent=2))
 print(name,version,run.returncode,text[-400:] if run.returncode else '',flush=True)
 if run.returncode:raise RuntimeError(name+' failed')
