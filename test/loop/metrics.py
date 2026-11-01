"""Read only the four authorized loop cases. No solver or other regression launches."""
from pathlib import Path
import xml.etree.ElementTree as ET
import json,math,csv,re,hashlib
base=Path(__file__).resolve().parent

def vtk(path):
 t=ET.parse(path);a={x.get('Name'):list(map(float,x.text.split())) for x in t.findall('.//PointData/DataArray')}
 q=list(map(float,t.find('.//Points/DataArray').text.split()));return list(zip(q[::3],q[1::3],q[2::3])),a

def frames(folder):
 return [(float(x.get('timestep')),folder/x.get('file')) for x in ET.parse(next(folder.glob('*.pvd'))).findall('.//DataSet')]

def weights(points):
 axes=[sorted(set(p[a] for p in points)) for a in range(3)];w=[]
 for axis in axes:
  w.append({x:((axis[min(i+1,len(axis)-1)]-axis[max(i-1,0)])/2 if len(axis)>1 else 1) for i,x in enumerate(axis)})
 return [w[0][x]*w[1][y]*w[2][z] for x,y,z in points],list(map(len,axes))

def csvtrace(p,context=False):
 with p.open() as h:
  meta=next(h).strip().split(',') if context else None
  return meta,[{k:(int(v) if k in {'row','column','cell','axis','lower','upper','i','j','k','id','owner','solved'} else float(v)) for k,v in r.items()} for r in csv.DictReader(h)]

def traces(path,word,context=True):
 out=[]
 for p in path.glob('*.csv'):
  if p.name.startswith('._') or word not in p.name:continue
  meta,rows=csvtrace(p,context);out.append((p,meta,rows))
 return sorted(out,key=lambda r:int(r[1][3]) if context else int(r[0].name.split('-')[2]))

result={'density':{},'pressure':{},'scope':'four authorized serial cases, only their output/trace comparisons'}
for name in ['baseline','feedback']:
 case=base/'density'/name;folder=case/'result/trace-04';series=frames(folder);points,a=vtk(series[-1][1]);w,shape=weights(points);V=math.fsum(w);beta=.2 if name=='feedback' else 0
 integrals=[];mxErr=pErr=rhoErr=energyErr=0
 for t,p in series:
  _,v=vtk(p);expectedU=.25*math.exp(-beta*t);expectedE=2.5+.0375*math.exp(-2*beta*t)
  mxErr=max(mxErr,max(abs(x/1.2-expectedU) for x in v['MomentumX']));pErr=max(pErr,max(abs(x-1) for x in v['p']));rhoErr=max(rhoErr,max(abs(x-1.2) for x in v['rho']));energyErr=max(energyErr,max(abs(x-expectedE) for x in v['TotalEnergyDensity']))
  integrals.append({'time':t,'mass':math.fsum(x*y for x,y in zip(w,v['rho'])),'momentumX':math.fsum(x*y for x,y in zip(w,v['MomentumX'])),'energy':math.fsum(x*y for x,y in zip(w,v['TotalEnergyDensity']))})
 log=(case/'run-03.log').read_text();lines=[x for x in log.splitlines() if x.startswith('[SF PLAN]')];commits=[x for x in lines if 'op=time.commit' in x];dt=[float(re.search(r'dt=([^ ]+)',x)[1]) for x in commits]
 cp=folder/'checkpoints';q=traces(cp,'stage boundary-ready Q');r=traces(cp,'local residual before GlobalDof SUM');g=traces(cp,'global residual after GlobalDof SUM')
 perStage=[];rates=[]
 for (qf,qmeta,qr),(rf,rmeta,rr),(gf,gmeta,gr) in zip(q,r,g):
  qm={(x['x'],x['y'],x['z']):x for x in qr if x['interior']==1};errs=[];es=[]
  for x in rr:
   v=qm[(x['x'],x['y'],x['z'])];expected=[0,beta*v['Q1'],0,0,beta*v['Q1']**2/v['Q0']]
   errs.append(max(abs(x['RL'+str(c)]-expected[c]) for c in range(5)))
  globalCopy=max(abs(x['RL'+str(c)]-(y['RG'+str(c)] if y['hasGlobal'] else y['RL'+str(c)])) for x,y in zip(rr,gr) for c in range(5))
  rates.append([-math.fsum(x['RL'+str(c)] for x in rr)/len(rr)*V for c in range(5)])
  perStage.append({'stageIndex':int(qmeta[3]),'stageTime':float(qmeta[5]),'maxResidualVsDeclaredSource':max(errs),'localVsEffectiveGlobalResidual':globalCopy})
 # First-step actual RK quadrature, volume integral of the uniform manufactured flow.
 delta=[integrals[1][key]-integrals[0][key] for key in ['mass','momentumX','energy']];rk=[dt[0]/6*sum(b*rates[i][c] for i,b in enumerate([1,2,2,1])) for c in [0,1,4]]
 fluxBefore=traces(cp,'candidate face flux before canonical COPY');fluxAfter=traces(cp,'face flux after canonical COPY')
 faceDifference=max(abs(x['F'+str(c)]-y['F'+str(c)]) for (_,_,br),(_,_,ar) in zip(fluxBefore,fluxAfter) for x,y in zip(br,ar) for c in range(5))
 stageOps=[x for x in lines if ('op=flow.stage.rhs' in x or 'op=explicit.stage.execute' in x)]
 expectedRhs=40*4
 # Underlying fields must agree with the earlier full run, despite read-only tracing.
 _,old=vtk(frames(case/'result')[-1][1]);traceInvariant=old==a
 impulses=[];works=[];amplitude=.3
 for h in dt:
  k1=-beta*amplitude;m2=amplitude+h*k1/2;k2=-beta*m2;m3=amplitude+h*k2/2;k3=-beta*m3;m4=amplitude+h*k3;k4=-beta*m4
  impulses.append(V*h*(k1+2*k2+2*k3+k4)/6)
  works.append(-V*h*beta*(amplitude**2+2*m2**2+2*m3**2+m4**2)/1.2/6)
  amplitude+=h*(k1+2*k2+2*k3+k4)/6
 fullBudget={'impulse':math.fsum(impulses),'mechanicalWork':math.fsum(works),'momentumError':abs(integrals[-1]['momentumX']-integrals[0]['momentumX']-math.fsum(impulses)),'energyError':abs(integrals[-1]['energy']-integrals[0]['energy']-math.fsum(works)),'reference':'independent RK4 of declared uniform source; first step additionally checked against actual residual snapshots'}
 data={'fullTimeSourceBudget':fullBudget,'shapeNodes':shape,'volume':V,'physicalFrames':len(series),'finalTime':series[-1][0],'physicalSteps':len(commits),'dtMin':min(dt),'dtMax':max(dt),'RhsEvaluations':len(stageOps),'uMaxAnalyticError':mxErr,'pMaxChange':pErr,'rhoMaxChange':rhoErr,'energyDensityMaxAnalyticError':energyErr,'integralsInitial':integrals[0],'integralsFinal':integrals[-1],'firstStepRkBudgetError':max(abs(x-y) for x,y in zip(delta,rk)),'firstStepActualIntegratedRhsQuadrature':rk,'firstStepStageTrace':perStage,'canonicalBeforeAfterMaxDifference':faceDifference,'traceVsRun02FieldsEqual':traceInvariant,'finalRhoRange':[min(a['rho']),max(a['rho'])],'finalPressureRange':[min(a['p']),max(a['p'])]}
 checks={'end':series[-1][0]==1 and len(commits)==40,'analytic':mxErr<=1e-9 and pErr<=1e-9 and rhoErr<=1e-9,'rhsOnce':len(stageOps)==expectedRhs,'stageSource':max(x['maxResidualVsDeclaredSource'] for x in perStage)<1e-11,'budget':data['firstStepRkBudgetError']<1e-10 and max(fullBudget['momentumError'],fullBudget['energyError'])<1e-9,'canonical':faceDifference==0,'traceInvariant':traceInvariant}
 grouped={}
 for line in lines:
  step=int(re.search(r'step=(\d+)',line)[1]);grouped.setdefault(step,[]).append(line.split(' op=')[1])
 if beta:
  stage=['time.stage.open','flow.stage.prepare','time.stage.ready','affine.solve@11:drag.update','affine.solve@11:work.update','flow.stage.rhs','time.stage.rhsReady','flow.stage.advance','time.stage.close']
  expected=['flow.step.prepare','time.group.dt','flow.step.begin']+stage*4+['time.group.publishReady','flow.step.commit','affine.publish@11:drag.update','affine.publish@11:work.update','time.commit']
  data['traceOrder']=all(ops==expected for ops in grouped.values())
  checks['traceOrder']=data['traceOrder']
 checks['stageTimes']=[x['stageTime'] for x in perStage]==[0,.0125,.0125,.025]
 data['checks']=checks;result['density'][name]=data

for name in ['baseline','feedback']:
 case=base/'pressure'/name;folder=case/'result/trace-04';series=frames(folder);points,a=vtk(series[-1][1]);w,shape=weights(points);beta=.2 if name=='feedback' else 0
 ref=[.5*(.25-y*y) if not beta else .1/beta*(1-math.cosh(math.sqrt(beta/.1)*y)/math.cosh(math.sqrt(beta/.1)*.5)) for x,y,z in points]
 u=a['U'][::3];l2=math.sqrt(sum(ww*(x-y)**2 for ww,x,y in zip(w,u,ref))/sum(ww*y*y for ww,y in zip(w,ref)))
 log=(case/'run-03.log').read_text();lines=[x for x in log.splitlines() if x.startswith('[SF PLAN]')];commit=[x for x in lines if 'op=time.commit' in x];dt=[float(re.search(r'dt=([^ ]+)',x)[1]) for x in commit]
 cp=folder/'assembly';paths=[p for p in cp.glob('state-2047-*.csv') if not p.name.startswith('._') and not p.name.endswith(('-faces.csv','-matrix.csv'))];paths.sort(key=lambda p:int(p.name.split('-')[2]))
 snap={p.name.rsplit('-',1)[-1].removesuffix('.csv'):p for p in paths} # select accepted/momentum explicitly below
 def read(p):return csvtrace(p)[1]
 momentum=next(p for p in paths if p.name.endswith('-momentum.csv'));accepted=next(p for p in paths if p.name.endswith('-accepted.csv'))
 m={int(v['cell']):v for v in read(momentum)};acc={int(v['cell']):v for v in read(accepted)}
 faces=lambda p:{(int(x['axis']),int(x['lower'])):x for x in read(p.with_name(p.stem+'-faces.csv'))}
 mf=faces(momentum);af=faces(accepted);ijk={tuple(int(v[c]) for c in ['i','j','k']):int(v['cell']) for v in m.values()}
 def neighbor(v,axis,sign):
  at=[int(v[c]) for c in ['i','j','k']];at[axis]+=sign;return ijk[tuple(at)]
 sourceError=momentumError=div=0;momentumSamples=[]
 for cell,v in m.items():
  if not v['solved']:continue
  z=acc[cell];net=0
  for axis in range(2):net+=af[(axis,cell)]['flux']-af[(axis,neighbor(v,axis,-1))]['flux']
  div=max(div,abs(net*z['inverseVolume']))
  for c in range(3):
   adv=lap=0
   for axis in range(2):
    for sign in [-1,1]:
     lower=cell if sign==1 else neighbor(v,axis,-1);f=mf[(axis,lower)];lo=m[lower];up=m[int(f['upper'])];area=[f[k] for k in ['Ax','Ay','Az']]
     flow=sum(.5*(lo['U'+str(k)]+up['U'+str(k)])*area[k] for k in range(3));adv+=sign*flow*(lo if flow>=0 else up)['U'+str(c)]
     nb=m[neighbor(v,axis,sign)];h=math.sqrt(sum((up[k]-lo[k])**2 for k in ['x','y','z']));lap+=math.sqrt(sum(x*x for x in area))*(nb['U'+str(c)]-v['U'+str(c)])/h
   spatialRhs=v['inverseVolume']*(-adv+.1*lap);source=-beta*v['U0'] if c==0 else 0
   inferred=(v['H'+str(c)]-v['U'+str(c)])/v['dt']-spatialRhs;sourceError=max(sourceError,abs(inferred-source))
   grad=0
   if c<2:
    lo=acc[neighbor(v,c,-1)];up=acc[neighbor(v,c,1)];grad=(up['p']-lo['p'])/(up[['x','y','z'][c]]-lo[['x','y','z'][c]])
   residual=(z['U'+str(c)]-v['U'+str(c)])/v['dt']-spatialRhs-source+grad
   momentumError=max(momentumError,abs(residual));momentumSamples.append(residual)
 consistency=0;velocityConsistency=0;previous=None
 for i,p in enumerate(paths):
  if not p.name.endswith('-flux-before.csv'):continue
  after=paths[i+1];before=faces(p);post=faces(after)
  now={int(x['cell']):x for x in read(p)}
  for cell,v in now.items():
   if not v['solved']:continue
   for c in range(3):
    deltaGrad=0
    if c<2:
     lo=now[neighbor(v,c,-1)];up=now[neighbor(v,c,1)];deltaGrad=(up['deltaP']-lo['deltaP'])/(up[['x','y','z'][c]]-lo[['x','y','z'][c]])
    if previous is None:
     pressureGrad=0
     if c<2:
      lo=m[neighbor(v,c,-1)];up=m[neighbor(v,c,1)];pressureGrad=(up['p']-lo['p'])/(up[['x','y','z'][c]]-lo[['x','y','z'][c]])
     initial=m[cell]['H'+str(c)]-m[cell]['rAU']*pressureGrad
    else:initial=previous[cell]['U'+str(c)]
    velocityConsistency=max(velocityConsistency,abs(v['U'+str(c)]-(initial-v['rAU']*deltaGrad)))
  previous={int(x['cell']):x for x in read(after)}
  for k,f in before.items():consistency=max(consistency,abs(post[k]['flux']-(f['flux']-f['response']*(f['deltaR']-f['deltaL']))))
 matrix=read(accepted.with_name(accepted.stem+'-matrix.csv'));solution={int(x['row']):x['solution'] for x in matrix};rhs={int(x['row']):x['rhs'] for x in matrix};Ax={row:0 for row in rhs}
 for x in matrix:Ax[int(x['row'])]+=x['coefficient']*solution[int(x['column'])]
 residual=[Ax[r]-b for r,b in rhs.items()];rnorm=math.sqrt(sum(x*x for x in residual));bnorm=math.sqrt(sum(x*x for x in rhs.values()))
 affine=sum('op=affine.solve@' in x for x in lines);continuity=[tuple(map(float,x)) for x in re.findall(r'maxDiv\(before/after\)=([^/\s]+)/([^\s]+)',log)]
 data={'shapeNodes':shape,'volume':sum(w),'finalTime':series[-1][0],'physicalSteps':len(commit),'dtMin':min(dt),'dtMax':max(dt),'profileRelativeWeightedL2':l2,'uMax':max(u),'uReferenceMax':max(ref),'transverseVelocityMax':max(abs(x) for x in a['U'][1::3]),'pressureRange':[min(a['p']),max(a['p'])],'continuityFinalFromCanonicalFaces':div,'continuityFinalLog':continuity[-1],'momentumResidualMax':momentumError,'momentumResidualL2':math.sqrt(sum(x*x for x in momentumSamples)/len(momentumSamples)),'momentumResidualMaxNormalizedByG':momentumError/.1,'sourceVsAssembledHMaxError':sourceError,'pressureMatrixResidualL2':rnorm,'pressureMatrixRhsL2':bnorm,'pressureMatrixRelativeResidual':rnorm/bnorm if bnorm else None,'velocityCorrectionIdentityMaxError':velocityConsistency,'lastFluxCorrectionIdentityMaxError':consistency,'affineEvaluations':affine,'matrixRows':len(rhs),'endReason':'endTime=10 reached at step 2048; endStep cap 2200 not reached'}
 grouped={}
 for line in lines:
  step=int(re.search(r'step=(\d+)',line)[1]);grouped.setdefault(step,[]).append(line.split(' op=')[1])
 passOps=['pressure.boundary.prepare','pressure.assemble','pressure.solve','velocity.correct','pressure.update.prepare','flux.correct','pressure.correction.commit']
 expected=['pressure.prepare']+(['affine.solve@14:damping.update'] if beta else [])+['momentum.assemble','momentum.solve']+passOps*2+['pressure.step.commit','time.commit']
 data['traceOrder']=all(ops==expected for ops in grouped.values())
 _,prior=vtk(frames(case/'result/trace-03')[-1][1]);data['finalFieldsVsTrace03Equal']=prior==a
 data['checks']={'traceInvariant':prior==a,'traceOrder':data['traceOrder'],'end':series[-1][0]==10 and len(commit)==2048,'profile':l2<=.05,'continuity':div<=1e-7,'momentum':momentumError<=1e-8,'source':sourceError<=1e-11,'flux':consistency<=1e-14 and velocityConsistency<=1e-12,'matrixAbsolute':rnorm<=1e-12,'frequency':affine==(2048 if beta else 0)}
 result['pressure'][name]=data

result['allPassed']=all(all(v['checks'].values()) for group in ['density','pressure'] for v in result[group].values())
(base/'evidence/metrics.json').write_text(json.dumps(result,indent=2)+'\n')
for group in ['density','pressure']:
 for name,data in result[group].items():print(group,name,json.dumps(data,indent=2))
print('ALL METRIC CHECKS:',result['allPassed'])
raise SystemExit(0 if result['allPassed'] else 1)
