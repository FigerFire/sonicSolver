from pathlib import Path
import json,shutil
base=Path(__file__).resolve().parent

def write(p,object,identity,body):
 p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps({'SonicFile':{'object':object,'type':identity},**body},indent=2)+'\n')
def make(name,recipe,ids,eqs):
 p=base/'cases'/name
 if p.exists():raise RuntimeError('Refusing overwrite '+str(p))
 write(p/'case.yaml','case','registry',{'formatVersion':1,'name':name})
 solver={key:{'type':role,'file':'solvers/'+key+'.yaml'} for key,role in [('runtime','runtime'),('numerics','numerics'),('state','stateRegistry'),('equationRegistry','equationRegistry'),('execution','algorithmRegistry'),('how','executionProgram'),('which','providerBindings')]}
 write(p/'solvers/solvers.yaml','solver','registry',solver)
 write(p/'solvers/runtime.yaml','solver','runtime',{'startTime':0,'endTime':.005,'maxDeltaT':.0005,'CFL':.5,'writeControl':'timeStep','writeInterval':1,'writeInitial':True,'parallel':{'enabled':False},'output':{'jobName':'nativePair','outputDir':'result'}})
 write(p/'solvers/numerics.yaml','solver','numerics',{'time':{'default':recipe}})
 write(p/'solvers/state.yaml','solver','state',{'use':ids})
 formula=lambda op,*args:{'op':op,'args':list(args)}
 add={e:{'lhs':formula('ddt',{'symbol':s}),'rhs':formula('diffusion',{'constant':d},{'symbol':s})} for s,e,d in zip(ids,eqs,[.03,.08,.11])}
 if name=='linearTermsRK4':
  add[eqs[0]]['rhs']=formula('add',add[eqs[0]]['rhs'],formula('multiply',{'constant':.2},{'symbol':ids[0]}))
  add[eqs[1]]['rhs']=formula('add',add[eqs[1]]['rhs'],{'constant':.15})
 write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',{'use':eqs,'add':add})
 write(p/'solvers/execution.yaml','solver','execution',{'Explicit':{}})
 calls=[{'kind':'EquationCall','equation':e,'target':s,'occurrence':'advance.'+s} for s,e in zip(ids,eqs)]
 write(p/'solvers/how.yaml','solver','how',{'root':{'kind':'Sequence','id':'native.diffusion','children':calls+[{'kind':'Commit','id':'physical.commit'}]}})
 write(p/'solvers/which.yaml','solver','which',{'bindings':[{'equation':e,'method':'ScalarDiffusionCentral2'} for e in eqs]})
 fields={s:{'new':True,'type':'scalar','output':True} for s in ids}
 fields.update({'internalField':{'type':'internalField','file':'fields/internalField.yaml'},'boundaries':{'type':'boundary','file':'fields/boundaries.yaml'}})
 write(p/'fields/fields.yaml','fields','registry',fields)
 init={ids[0]:{'default':0,'leftHalf':1},ids[1]:{'default':.5,'rightHalf':2}}
 if len(ids)>2:init[ids[2]]={'default':.75,'leftHalf':.2}
 write(p/'fields/internalField.yaml','fields','internalField',init)
 bcs={s:{side:{'type':'zeroGradient'} for side in ['Left','Right','Bottom','Top','Front','Back']} for s in ids}
 bcs[ids[0]].update({'Left':{'type':'fixedValue','value':1},'Right':{'type':'fixedValue','value':0}})
 write(p/'fields/boundaries.yaml','fields','boundaries',bcs)
 write(p/'models/models.yaml','models','registry',{})
 write(p/'mesh/mesh.yaml','mesh','registry',{'files':['mesh/mesh.sfm']})
 n=(33,5,5);point=lambda i,j,k:(k*n[1]+j)*n[0]+i
 lines=['#Information','33 5 5 1','','#Point']
 sets={s:[] for s in ['Left','Right','Bottom','Top','Front','Back','leftHalf','rightHalf']}
 for k in range(n[2]):
  for j in range(n[1]):
   for i in range(n[0]):
    lines.append(f'{i/32:.17g} {j/16:.17g} {k/16:.17g}')
    for s,yes in [('Left',i==0),('Right',i==32),('Bottom',j==0),('Top',j==4),('Front',k==0),('Back',k==4),('leftHalf',i<=16),('rightHalf',i>16)]:
     if yes:sets[s].append(point(i,j,k))
 for s,points in sets.items():lines+=['','#'+s,' '.join(map(str,points))]
 (p/'mesh/mesh.sfm').write_text('\n'.join(lines)+'\n')
 # Unregistered input must neither contribute a third PDE nor be parsed.
 (p/'solvers/unregistered.yaml').write_text('this is deliberately invalid YAML: [\n')
 return {'name':name,'recipe':recipe,'ids':ids,'equations':eqs,'diffusivities':[.03,.08,.11][:len(ids)],'path':str(p)}
CASES=[('pairEuler','forwardEuler',['dye','salt'],['spreadDye','spreadSalt']),('pairRK4','classicalRK4',['dye','salt'],['spreadDye','spreadSalt']),('renamedRK4','classicalRK4',['colourant','dopant'],['diffuseColour','diffuseDopant']),('tripleRK4','classicalRK4',['colourant','dopant','signal'],['diffuseColour','diffuseDopant','diffuseSignal'])]
CASES.append(('linearTermsRK4','classicalRK4',['dye','salt'],['spreadDye','spreadSalt']))
if __name__=='__main__':
 specs=[make(*values) for values in CASES]
 (base/'evidence').mkdir(exist_ok=True)
 (base/'evidence/cases.json').write_text(json.dumps(specs,indent=2))
