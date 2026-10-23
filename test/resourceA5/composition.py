from pathlib import Path
import json,shutil,os,subprocess,sys
base=Path(__file__).resolve().parent;root=base.parents[1]
def load(p):return json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def make(name,flow=True,reverse=False):
 p=base/'composition'/name
 if p.exists():return p
 shutil.copytree(base/'before/pairRK4',p,ignore=shutil.ignore_patterns('result','*.log','._*','.DS_Store'))
 registry=load(p/'solvers/solvers.yaml')
 if flow:
  registry['equations']={'type':'equationSystem','file':'solvers/equations.yaml'}
  save(p/'solvers/equations.yaml',{'SonicFile':{'object':'solver','type':'equations'},'preset':'compressible'})
  save(p/'models/models.yaml',{'SonicFile':{'object':'models','type':'registry'},'thermoDynamics':{'type':'thermoDynamics','file':'models/thermoDynamics.yaml'}})
  save(p/'models/thermoDynamics.yaml',{'SonicFile':{'object':'models','type':'thermoDynamics'},'thermoDynamics':{'equationOfState':'perfectGas','thermo':'hConst','transport':'const'},'properties':{'equationOfState':{'gamma':1.4,'R':287.05},'transport':{'mu':0,'Pr':.72}}})
  f=load(p/'fields/fields.yaml')
  for id,typ,dim in [('rho','scalar',[1,-3,0,0,0,0,0]),('U','vector',[0,1,-1,0,0,0,0]),('p','scalar',[1,-1,-2,0,0,0,0])]:f[id]={'type':typ,'dimensions':dim,'output':True}
  save(p/'fields/fields.yaml',f)
  f=load(p/'fields/internalField.yaml');f.update({'rho':{'default':1},'U':{'default':[.1,0,0]},'p':{'default':1}});save(p/'fields/internalField.yaml',f)
  f=load(p/'fields/boundaries.yaml')
  for id in ['rho','U','p']:f[id]={s:{'type':'zeroGradient'} for s in ['Left','Right','Bottom','Top','Front','Back']}
  save(p/'fields/boundaries.yaml',f)
  f=load(p/'solvers/state.yaml');f['use']=['rho','rhoU','rhoE']+f['use'];save(p/'solvers/state.yaml',f)
  f=load(p/'solvers/equationRegistry.yaml');f['use']=['Continuity','Momentum','Energy']+f['use'];save(p/'solvers/equationRegistry.yaml',f)
  f=load(p/'solvers/how.yaml');calls=[{'kind':'EquationCall','equation':eq,'target':target,'occurrence':'advance.'+target} for eq,target in [('continuity','rho'),('momentum','rhoU'),('energy','rhoE')]]
  f['root']['children']=(f['root']['children'][:-1]+calls if reverse else calls+f['root']['children'][:-1])+[f['root']['children'][-1]];save(p/'solvers/how.yaml',f)
  f=load(p/'solvers/which.yaml');f['bindings']=[{'equation':eq,'method':'ConservativeResidual'} for eq in ['continuity','momentum','energy']]+f['bindings'];save(p/'solvers/which.yaml',f)
  f=load(p/'solvers/numerics.yaml');f['convection']={'default':'WENO5','formulation':'conservativeFluxDifference','reconstruction':'characteristic','flux':'StegerWarming','interfaceFlux':'shared'};save(p/'solvers/numerics.yaml',f)
  f=p/'mesh/mesh.sfm';f.write_text(f.read_text().replace('33 5 5 1','33 5 5 3',1))
 # Two true stage-coupled scalar equations, with original D=0.03/0.08 and cross reaction.
 f=load(p/'solvers/equationRegistry.yaml')
 for eq,other in [('spreadDye','salt'),('spreadSalt','dye')]:
  rhs=f['add'][eq]['rhs'];f['add'][eq]['rhs']={'op':'add','args':[rhs,{'op':'multiply','args':[{'constant':.7,'dimensions':[0,0,-1,0,0,0,0]},{'symbol':other}]}]}
 save(p/'solvers/equationRegistry.yaml',f);save(p/'solvers/solvers.yaml',registry)
 return p
if __name__=='__main__':
 env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence'),'SF_PLAN_TRACE':'1','SF_SCALAR_TRACE':'1'}
 binary=base/'baseline/sonicSolver' if len(sys.argv)>1 and sys.argv[1]=='before' else root/'build/sonicSolver'
 for name,flow,reverse in [('flowPair',True,False),('reversedPair',True,True),('scalarPair',False,False)]:
  p=make(name,flow,reverse)
  for mode in ['explain','run']:
   r=subprocess.run([str(binary),mode,str(p)],cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   (p/(mode+('-before' if 'baseline' in str(binary) else '')+'.log')).write_bytes(r.stdout)
   print(name,mode,r.returncode,r.stdout[-1800:].decode() if r.returncode else '',flush=True)
   if r.returncode:break
