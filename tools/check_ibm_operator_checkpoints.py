#!/usr/bin/env python3
"""Validate opt-in IBM checkpoints for the canonical stationary cylinder methods."""
from pathlib import Path
import argparse,csv,json,math
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory",type=Path)
parser.add_argument("--center",type=float,nargs=3,required=True)
parser.add_argument("--json",type=Path,required=True)
args=parser.parse_args()
base=args.directory
def read(p):
 with p.open() as f:return float(next(f).split(',')[1]),list(csv.DictReader(f))
def entries(case,label,t):
 return [rows for p in (case/'checkpoints').glob('*-'+label+'.csv') for time,rows in [read(p)] if time==t]
def state(case,label,t):
 data={}
 for rows in entries(case,label,t):
  for r in rows:
   key=tuple(float(r[a]) for a in ['x','y','z'])
   if key not in data or r['owner']=='1':data[key]=r
 return data
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
summary=[]
for case in sorted(base.iterdir()):
 if not case.is_dir():continue
 tag='peskin' if 'Peskin' in case.name else 'fts' if 'Velocity' in case.name else 'fractional';steps=[]
 times=sorted({t for p in (case/'checkpoints').glob('*-constraint-result.csv') for t,_ in [read(p)]})
 for t in times:
  result=[r for rows in entries(case,'constraint-result',t) for r in rows];r=result[0]
  if any(x!=r for x in result):raise RuntimeError('global result replicas differ')
  before=state(case,tag+'-predictor',t);after=state(case,tag+'-corrected',t)
  label='peskin-consumed-force' if tag=='peskin' else tag+'-force'
  forces=[r for rows in entries(case,label,t) for r in rows if r['owner']=='1']
  bodyForce=[0.]*3;torque=[0.]*3;power=0.;energy=0.
  for f in forces:
   point=tuple(float(f[a]) for a in ['x','y','z']);vol=float(f['volume']);fv=[float(f[a]) for a in ['Fx','Fy','Fz']]
   qa,qb=before[point],after[point];va=[float(qa['Q'+str(i+1)])/float(qa['Q0']) for i in range(3)];vb=[float(qb['Q'+str(i+1)])/float(qb['Q0']) for i in range(3)]
   force=[-vol*v for v in fv];moment=cross([point[i]-args.center[i] for i in range(3)],force)
   for i in range(3):bodyForce[i]+=force[i];torque[i]+=moment[i]
   power+=vol*sum(fv[i]*(va[i]+vb[i])*.5 for i in range(3))
   energy+=vol*(float(qb['Q4'])-float(qa['Q4']))
  dt=t-(steps[-1]['time'] if steps else 0.)
  metrics=dict(time=t,nodalDomainVolume=math.fsum(float(f['volume']) for f in forces),forceIntegralError=max(abs(bodyForce[i]-float(r['force'+a])) for i,a in enumerate(['X','Y','Z'])),torqueIntegralError=max(abs(torque[i]-float(r['torque'+a])) for i,a in enumerate(['X','Y','Z'])),powerIntegralError=abs(power-float(r['power'])),mechanicalEnergyIncrementError=abs(energy-dt*power),constraintResidual=float(r['residual']),force=bodyForce,torque=torque,power=power)
  if tag!='fractional':
   edges=[r for rows in entries(case,'edges',t) for r in rows];keys=[(r['marker'],r['eulerian']) for r in edges]
   if len(keys)!=len(set(keys)):raise RuntimeError('duplicate global edge')
   adj=[r for rows in entries(case,'weighted-adjoint',t) for r in rows];metrics.update(markers=len({r['marker'] for r in edges}),uniqueEdges=len(edges),partitionOfUnityError=max(float(r['unityError']) for r in adj),weightedAdjointError=max(float(r['absoluteError']) for r in adj))
   lambdas={};lambdaConflicts=0
   for rows in entries(case,tag+'-lambda',t):
    for row in rows:
     marker=int(row['marker']);lv=tuple(float(row[a]) for a in ['xValue','yValue','zValue'])
     if marker in lambdas and lambdas[marker]!=lv:lambdaConflicts+=1
     lambdas[marker]=lv
   spread={}
   for row in edges:
    p=tuple(float(row[a]) for a in ['x','y','z']);scale=float(row['weight'])
    # surface measure belongs to marker, available in lambda record.
    marker=int(row['marker']);measure=float(next(x['measure'] for rows in entries(case,tag+'-lambda',t) for x in rows if int(x['marker'])==marker));scale=float(row['weight'])*measure/float(row['volume'])
    value=spread.setdefault(p,[0.,0.,0.])
    for i in range(3):value[i]+=scale*lambdas[marker][i]
   actualLabel='peskin-next-force' if tag=='peskin' else tag+'-force'
   actual=state(case,actualLabel,t)
   metrics['lambdaReplicaConflicts']=lambdaConflicts
   metrics['actualAdjointSpreadingError']=max(abs(float(actual[p][a])-v[i]) for p,v in spread.items() for i,a in enumerate(['Fx','Fy','Fz']))
   if lambdaConflicts:raise RuntimeError('canonical lambda COPY mismatch')
   local={};global_values={}
   for rows in entries(case,'local-mass-response',t):
    for row in rows:local.setdefault(int(row['marker']),[]).append(float(row['xValue']))
   for rows in entries(case,'global-mass-response',t):
    for row in rows:
     value=float(row['xValue']);marker=int(row['marker'])
     if not math.isfinite(value) or value<=0:raise RuntimeError('invalid reduced mass response')
     if marker in global_values and global_values[marker]!=value:raise RuntimeError('reduced mass response replicas differ')
     global_values[marker]=value
   if any(not math.isfinite(v) or v<0 for values in local.values() for v in values):raise RuntimeError('invalid local contribution')
   metrics['zeroLocalMassContributions']=sum(v==0 for values in local.values() for v in values)
   metrics['massResponseSumError']=max(abs(math.fsum(values)-global_values[marker]) for marker,values in local.items())
   moments={}
   for row in edges:
    marker=int(row['marker']);value=moments.setdefault(marker,[0.,0.,0.])
    for i,a in enumerate(['x','y','z']):value[i]+=float(row['weight'])*float(row[a])
   coordinates={int(row['marker']):[float(row[a]) for a in ['x','y','z']] for rows in entries(case,tag+'-lambda',t) for row in rows}
   metrics['interpolationFirstMomentError']=max(abs(moments[marker][i]-coordinates[marker][i]) for marker in moments for i in range(3))
  steps.append(metrics)
 summary.append(dict(case=case.name,steps=steps));print(case.name,steps[-1])
args.json.write_text(json.dumps(summary,indent=2))
