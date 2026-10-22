from pathlib import Path
import os,json,subprocess,time,shutil,hashlib
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'TMPDIR':str(base/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
binary=root/'build/sonicSolver';records=[]
def run(case,mode,expected=0,reason=''):
 start=time.monotonic();result=subprocess.run([str(binary),mode,str(case)],cwd=case,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (case/(mode+'.log')).write_text(result.stdout)
 entry=dict(case=str(case.relative_to(root)),mode=mode,command=[str(binary),mode,str(case)],seconds=time.monotonic()-start,returncode=result.returncode,expected=expected,reason=reason)
 records.append(entry)
 assert (result.returncode==0)==(expected==0),(entry,result.stdout[-2000:])
 if reason:assert reason in result.stdout,(entry,result.stdout[-2000:])
for name in ['pairRK4','linearTermsRK4']:
 for mode in ['explain','run']:run(base/'after'/name,mode)
for name,reason in [('unspecified','Unknown declared mathematical symbol: dye'),('coefficientUnits','Incompatible mathematical signatures'),('unknownOperator','Unknown mathematical operator')]:
 case=base/'negative'/name
 if not case.exists():shutil.copytree(base/'after/pairRK4',case,ignore=shutil.ignore_patterns('result','*.log','._*'))
 if name=='unspecified':
  p=case/'fields/fields.yaml';d=json.loads(p.read_text());d['dye'].pop('dimensions',None);p.write_text(json.dumps(d,indent=2))
 else:
  p=case/'solvers/equationRegistry.yaml';d=json.loads(p.read_text())
  def edit(x):
   if not isinstance(x,dict):return
   if x.get('op')=='diffusion':
    if name=='coefficientUnits':x['args'][0]['dimensions']=[0]*7
    else:x['op']='unregisteredDifferentialOperator'
   for v in x.values():
    if isinstance(v,list):
     for a in v:edit(a)
    else:edit(v)
  edit(d);p.write_text(json.dumps(d,indent=2))
 run(case,'explain',1,reason)
(base/'evidence/native-final.json').write_text(json.dumps(dict(binary=hashlib.sha256(binary.read_bytes()).hexdigest(),runs=records),indent=2))
print(json.dumps(records,indent=2))
