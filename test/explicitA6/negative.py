from pathlib import Path
import json,os,subprocess,shutil,sys
base=Path(__file__).resolve().parent;root=base.parents[1];binary=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else root/'build/sonicSolver';records=[]
env={**os.environ,'TMPDIR':str(base/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools'}
def missing_metadata(function,path,error):
 if not isinstance(error[1],FileNotFoundError):raise error[1]
def negative(name,source,file,edit,expected,action='explain'):
 p=base/'negative'/name
 if p.exists():shutil.rmtree(p,onerror=missing_metadata)
 shutil.copytree(base/'native'/source,p,ignore=shutil.ignore_patterns('result','*.log','._*'))
 path=p/file;x=json.loads(path.read_text());edit(x);path.write_text(json.dumps(x,indent=2))
 cmd=[str(binary),action,str(p)];r=subprocess.run(cmd,cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);log=r.stdout.decode();(p/'run.log').write_text(log)
 assert r.returncode and expected in log,(name,log[-2000:]);assert '[SF] Scalar step:' not in log,name
 records.append(dict(name=name,expected=expected,command=cmd,exit=r.returncode));print(name,'rejected')
negative('partnerBoundary','cross','fields/boundaries.yaml',lambda x:x['b'].pop('Right'),'Missing scalar boundary backing: b','run')
negative('vectorComponent','vector','fields/internalField.yaml',lambda x:x['Vexperiment'].update(default=[1,2]),'explicit vector components required','run')
negative('vectorResource','vector','fields/fields.yaml',lambda x:x.pop('internalField'),'missing initialization backing','run')
negative('wrongDimension','vector','solvers/equationRegistry.yaml',lambda x:x['add']['equation.Vexperiment']['rhs']['args'][0]['args'][0].update(dimensions=[0,0,-1,0,0,0,0]),'Incompatible mathematical signatures')
negative('variableD','multi','solvers/equationRegistry.yaml',lambda x:x['add']['equation.a']['rhs']['args'][0]['args'].__setitem__(0,{'op':'multiply','args':[{'constant':.03,'dimensions':[0,2,-1,0,0,0,0]},{'symbol':'a'}]}),'Unsupported scalar diffusivity')
negative('unsupportedTime','reaction','solvers/numerics.yaml',lambda x:x['time'].update(default='SSPRK3'),'cannot compile occurrence')
negative('unsupportedBackend','reaction','solvers/runtime.yaml',lambda x:x.update(parallel={'enabled':True,'nProcs':2}),'serial')
(base/'evidence/negative.json').write_text(json.dumps(records,indent=2))
