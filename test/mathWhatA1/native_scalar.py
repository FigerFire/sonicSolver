from pathlib import Path
import subprocess,shutil,json,os,hashlib,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent;root=base.parents[1]
env={**os.environ,'TMPDIR':str(base/'evidence'),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
for label in ['before','after']:
 for name in ['pairRK4','linearTermsRK4']:
  p=base/label/name
  if p.exists():raise RuntimeError('Refusing overwrite '+str(p))
  shutil.copytree(root/'test/nativeEquationClosure/cases'/name,p,ignore=shutil.ignore_patterns('result','reference','*.log','._*'))
  if label=='after':
   field=p/'fields/fields.yaml';doc=json.loads(field.read_text())
   for key,val in doc.items():
    if isinstance(val,dict) and val.get('new'):val['dimensions']=[0]*7
   field.write_text(json.dumps(doc,indent=2))
   eq=p/'solvers/equationRegistry.yaml';doc=json.loads(eq.read_text())
   def annotate(node):
    if not isinstance(node,dict):return
    if node.get('op')=='diffusion':node['args'][0]['dimensions']=[0,2,-1,0,0,0,0]
    if node.get('op')=='multiply' and 'constant' in node['args'][0]:node['args'][0]['dimensions']=[0,0,-1,0,0,0,0]
    if node.get('op')=='add' and 'constant' in node['args'][1]:node['args'][1]['dimensions']=[0,0,-1,0,0,0,0]
    for v in node.values():
     if isinstance(v,list):
      for x in v:annotate(x)
     else:annotate(v)
   annotate(doc);eq.write_text(json.dumps(doc,indent=2))
  binary=base/'baseline/sonicSolver' if label=='before' else root/'build/sonicSolver'
  for mode in ['explain','run']:
   r=subprocess.run([str(binary),mode,str(p)],cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
   (p/(mode+'.log')).write_text(r.stdout)
   print(label,name,mode,r.returncode,r.stdout[-500:] if r.returncode else '',flush=True)
   if r.returncode:break
