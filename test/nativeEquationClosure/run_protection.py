from pathlib import Path
import json,subprocess,os,shutil,hashlib,sys,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent;root=base.parents[1]
manifest=[];source=root/'test/stageCoefficientOwnership/verified'
label=sys.argv[1] if len(sys.argv)>1 else 'protection'
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
specs=['sod','eulerianSST','stageSST','ghostSST']
for name in specs:
 records=[]
 for version,binary in [('before',base/'baseline/sonicSolver'),('after',root/'build/sonicSolver')]:
  p=base/label/version/name
  if p.exists():raise RuntimeError('Refusing overwrite '+str(p))
  shutil.copytree(source/name,p,ignore=shutil.ignore_patterns('result','checkpoint','*.log','numerical.txt','._*'))
  cmd=[str(binary),'run',str(p)];actual={**env,'SF_EULERIAN_TURBULENCE_EVIDENCE':str(p/'numerical.txt')}
  r=subprocess.run(cmd,cwd=p,env=actual,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(p/'run.log').write_text(r.stdout)
  assert r.returncode==0,(name,version,r.stdout[-2000:])
  records.append({'case':str(p),'command':cmd,'binary':hashlib.sha256(binary.read_bytes()).hexdigest(),'exit':r.returncode})
 a,b=[Path(x['case']) for x in records]
 files=lambda p:{f.relative_to(p/'result'):f.read_bytes() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._')}
 assert files(a)==files(b),(name,'first output difference')
 if (a/'numerical.txt').exists():assert (a/'numerical.txt').read_bytes()==(b/'numerical.txt').read_bytes()
 trace=lambda p:[x for x in (p/'run.log').read_text().splitlines() if x.startswith('[SF PLAN]') or x.startswith('[SF] Step time:') or x.startswith('[SF] State closure:') or x.startswith('[SF] Eulerian step:')]
 assert trace(a)==trace(b),(name,'first trace difference')
 times=[x.attrib['timestep'] for x in ET.parse(next(f for f in (b/'result').glob('*.pvd') if not f.name.startswith('._'))).findall('.//DataSet')]
 manifest.append({'name':name,'source':str(source/name),'runs':records,'frames':len(times),'times':times,'outputsExact':True,'traceExact':True})
 (base/'evidence'/(label+'.json')).write_text(json.dumps(manifest,indent=2))
 print(name,'exact',len(times),'frames',flush=True)
