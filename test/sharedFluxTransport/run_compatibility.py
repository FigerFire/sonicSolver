from pathlib import Path
import shutil,re,subprocess,json,hashlib,sys,os
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from thermophysical_fixture import declare_historical_thermophysical_fixture
root=Path(__file__).resolve().parents[2];name=sys.argv[1] if len(sys.argv)>1 else 'compatibility'
if Path(name).name!=name:raise ValueError('Evidence output must be a child directory name')
out=root/'test/sharedFluxTransport'/name;out.mkdir(exist_ok=True)
executables={'before':root/'test/sharedFluxTransport/evidence/sonicSolver-before','after':root/'build/sonicSolver'}
specs=[('sod','test/Sod/sodCase',.01),('weno7','test/Sod/sodCase_weno7',.01),('ghost','test/IBM/cylinderFlowViscousGhost',.001),('sstGhost','test/IBM/cylinderFlowSSTGhost',.001),('legacySST','test/turbulence/singleFluidSST',.003)]
manifest=[]
for name,source,end in specs:
    pair=[]
    for version,exe in executables.items():
        case=out/(name+'-'+version)
        if case.exists():raise RuntimeError('Do not overwrite evidence: '+str(case))
        shutil.copytree(root/source,case,ignore=shutil.ignore_patterns('result','results','*.log','._*','.DS_Store','__pycache__'))
        declare_historical_thermophysical_fixture(case)
        path=case/'solvers/runtime.yaml';text=path.read_text();text=re.sub(r'^endTime:.*$',f'endTime: {end}',text,flags=re.M);text=re.sub(r'^writeInterval:.*$',f'writeInterval: {end/2}',text,flags=re.M);path.write_text(text)
        command=[str(exe),'run',str(case)];log=case/'run.log'
        with log.open('w') as stream:run=subprocess.run(command,cwd=root,stdout=stream,stderr=subprocess.STDOUT,env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(root/'test/sharedFluxTransport/evidence')})
        entry={'case':name,'version':version,'source':source,'directory':str(case),'endTime':end,'modified':['runtime endTime/writeInterval','explicit historical thermophysical fixture parameters where formerly omitted'],'command':command,'executableSHA256':hashlib.sha256(exe.read_bytes()).hexdigest(),'exit':run.returncode}
        manifest.append(entry);(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
        print(name,version,run.returncode,flush=True);pair.append(case)
    a,b=pair;files=[p for p in (a/'result').rglob('*.vts') if not p.name.startswith('._')];other=[p for p in (b/'result').rglob('*.vts') if not p.name.startswith('._')];equal=len(files)==len(other) and bool(files)
    for path in files:
        peer=b/'result'/path.relative_to(a/'result')
        equal=equal and peer.exists() and path.read_bytes()==peer.read_bytes()
    print(name,'VTS byte equality',equal,'frames',len(files),flush=True)
    manifest[-1]['byteIdenticalVTS']=equal;manifest[-1]['physicalVTSFiles']=len(files);(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
