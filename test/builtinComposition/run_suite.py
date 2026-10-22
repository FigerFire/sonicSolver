from pathlib import Path
import os,subprocess,shutil,json,hashlib,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1];bin=root/'build/sonicSolver';tool=base/'bin/test_builtinComposition'
runs=base/'evidence/runs';runs.mkdir(parents=True,exist_ok=True)
run=runs/f'{len(list(runs.iterdir()))+1:03d}';run.mkdir()
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1'}
records=[];commands=[]
def command(args,cwd,log,success=True):
 r=subprocess.run([str(a) for a in args],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 log.write_text(r.stdout)
 commands.append({'args':[str(a) for a in args],'cwd':str(cwd),'log':str(log),'exit':r.returncode})
 if success:assert r.returncode==0,(args,r.stdout[-2200:])
 return r
source=root/'test/turbulence/singleFluidSST'
def create(name,source=source):
 p=run/name
 shutil.copytree(source,p,ignore=shutil.ignore_patterns('result','*.log','numerical.txt','._*','reference','inspect'))
 command([tool,'prepare',source,p,name],p,p/'prepare.log')
 return p
fields=lambda p:{f.relative_to(p/'result'):f.read_bytes() for f in (p/'result').glob('*') if f.is_file() and not f.name.startswith('._')}
positives=['laminarAbsent','laminarExplicit','laminarPreset','sstExplicit','sstDefaultOn','sstPreset','sstDisabledPreset','laminarOverridePreset','overridePreset','resolvedDefaultsPreset','gravityOn','gravityOff','gravityDefaultOn','gravityZone']
for name in positives:
 p=create(name)
 command([tool,p,p/'inspect','inspect'],p,p/'inspect.log')
 for mode in ['check','explain','run']:command([bin,mode,p],p,p/(mode+'.log'))
 command([tool,p,p/'reference','reference'],p,p/'reference.log')
 assert fields(p)==fields(p/'reference'),(name,'independent manual numerical result mismatch')
 pvd=next(f for f in (p/'result').glob('*.pvd') if not f.name.startswith('._'))
 times=[float(d.attrib['timestep']) for d in E.parse(pvd).findall('.//DataSet')]
 records.append({'case':str(p),'name':name,'frames':len(times),'finalTime':times[-1],'manualExact':True})
 (run/'manifest.json').write_text(json.dumps({'binary':hashlib.sha256(bin.read_bytes()).hexdigest(),'positive':records},indent=2))
 print(name,len(times),'frames; independent IR/provider/plan/solution exact',flush=True)
for a,b in [('laminarAbsent','laminarExplicit'),('laminarExplicit','laminarPreset'),('sstExplicit','sstDefaultOn'),('sstDefaultOn','sstPreset'),('gravityOn','gravityDefaultOn'),('gravityOn','gravityZone'),('sstDisabledPreset','laminarOverridePreset')]:
 assert fields(run/a)==fields(run/b),(a,b,'equivalent selection output differs')
for name in ['laminarAbsent','laminarExplicit','laminarPreset','sstDisabledPreset','laminarOverridePreset','gravityOff']:
 s=json.loads((run/name/'inspect/signature.json').read_text())
 assert 'k' not in s['math'] and 'mu_t' not in s['state']
 assert all(c['provider']!='flow.turbulence' for c in s['calls'].values())
assert json.loads((run/'sstDisabledPreset/inspect/signature.json').read_text())['state']['mu']==json.loads((run/'sstPreset/inspect/signature.json').read_text())['state']['mu']
for name in ['gravityOn','gravityDefaultOn','gravityZone','gravityOff']:
 p=run/name;pvd=next(f for f in (p/'result').glob('*.pvd') if not f.name.startswith('._'));vts=p/'result'/E.parse(pvd).findall('.//DataSet')[-1].get('file')
 u=next(a for a in E.parse(vts).findall('.//PointData/DataArray') if a.get('Name')=='U');values=list(map(float,u.text.split()))[::3]
 expected=.3+(9.8*.001 if name!='gravityOff' else 0)
 coordinates=list(map(float,E.parse(vts).find('.//Points/DataArray').text.split()))
 interior=[value for i,value in enumerate(values) if coordinates[3*i]==50 and coordinates[3*i+1]==10]
 assert interior and max(abs(x-expected) for x in interior)<1e-10,(name,interior,expected)
 next(row for row in records if row['name']==name)['centralAccelerationCheck']={'coordinates':[50,10],'expectedUx':expected,'actualUx':interior,'expectedDeltaRhoUx':0. if name=='gravityOff' else 1.2*9.8*.001}
override=(run/'overridePreset/explain.log').read_text()
assert 'overrides default' in override and 'solver/numerics/time/default <- explicit' in override
sig=json.loads((run/'overridePreset/inspect/signature.json').read_text())
assert any(item[-1]=='weno7Rusanov' for item in sig['operators'])
assert 'configured recipe: forwardEuler' in override and 'compiled recipe stages, count=1' in override
assert 'object/turbulence/RAS/model' not in (run/'laminarOverridePreset/explain.log').read_text()
for name in ['laminarAbsent','laminarExplicit','laminarPreset','sstDisabledPreset','laminarOverridePreset']:
 p=run/name;pvd=next(f for f in (p/'result').glob('*.pvd') if not f.name.startswith('._'));last=p/'result'/E.parse(pvd).findall('.//DataSet')[-1].get('file')
 names={a.get('Name') for a in E.parse(last).findall('.//PointData/DataArray')}
 assert not ({'omega','TurbulenceOmega','TurbulenceK','TurbulentViscosity'} & names),(name,'inactive output resource exists')
for name in ['sstExplicit','sstDefaultOn','sstPreset','overridePreset','resolvedDefaultsPreset']:
 p=run/name;pvd=next(f for f in (p/'result').glob('*.pvd') if not f.name.startswith('._'))
 final=p/'result'/E.parse(pvd).findall('.//DataSet')[-1].get('file')
 omega=next(a for a in E.parse(final).findall('.//PointData/DataArray') if a.get('Name')=='omega')
 values=list(map(float,omega.text.split()))
 assert min(values)>0 and max(values)<1,(name,'selected SST did not update its original physical storage')
for name,src,expected in [
 ('ambiguity','sstPreset','Ambiguous preset default'),('conflict','sstExplicit','Conflicting explicit target'),
 ('time','sstPreset','time recipe'),('initial','sstPreset','initial'),('boundary','sstPreset','boundary'),
 ('model','sstPreset','unknown composition preset'),('enable','gravityOn',"'enabled'"),
 ('source','gravityOn','requires a finite'),('sourceAlias','gravityOn','Conflicting explicit gravity'),('mpi','sstPreset','MPI'),('binding','sstPreset','explicit nonempty'),('output','laminarAbsent','No runtime field provider'),('field','sstPreset','lacks registered field')]:
 p=create('negative.'+name,run/src)
 r=command([bin,'run',p],p,p/'rejection.log',False)
 assert r.returncode!=0,(name,'invalid input executed')
 assert expected.lower() in r.stdout.lower(),(name,expected,r.stdout[-1800:])
 assert not list((p/'result').glob('*.vts')),(name,'failed case published physical output')
 records.append({'name':'negative.'+name,'case':str(p),'rejected':True,'diagnostic':r.stdout})
print('14 positive / 13 negative cases passed',flush=True)
(run/'manifest.json').write_text(json.dumps({'binary':hashlib.sha256(bin.read_bytes()).hexdigest(),'cases':records,'commands':commands},indent=2))
(base/'evidence/latest.json').write_text(json.dumps({'run':str(run),'manifest':str(run/'manifest.json')},indent=2))
