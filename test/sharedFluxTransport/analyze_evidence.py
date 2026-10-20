from pathlib import Path
import json,re,hashlib,math,xml.etree.ElementTree as ET
root=Path(__file__).resolve().parents[2];base=root/'test/sharedFluxTransport'
def files(case):return sorted(p for p in (case/'result').rglob('*.vts') if not p.name.startswith('._'))
def times(case):
 log=(case/'run.log').read_text();return [(float(t),float(dt)) for t,dt in re.findall(r'Step time: time=([\deE+.-]+), dt=([\deE+.-]+)',log)]
def arrays(path):
 x=ET.parse(path);coords=[float(v) for v in x.find('.//Points/DataArray').text.split()];fields={}
 for a in x.findall('.//PointData/DataArray'):
  values=[float(v) for v in (a.text or '').split()];components=int(a.get('NumberOfComponents','1'));fields[a.get('Name')]=(components,values)
 return [tuple(coords[i:i+3]) for i in range(0,len(coords),3)],fields
def field_metrics(paths):
 points={};replica_error=0
 for path in paths:
  coords,data=arrays(path)
  for node,xyz in enumerate(coords):
   row={name:tuple(values[node*n:(node+1)*n]) for name,(n,values) in data.items()}
   if xyz in points:
    for name in ('rho','p','MomentumX','MomentumY','MomentumZ','TotalEnergyDensity'):
     if name in row:replica_error=max(replica_error,max(abs(a-b) for a,b in zip(row[name],points[xyz][name])))
   else:points[xyz]=row
 axes=[sorted(set(x[a] for x in points)) for a in range(3)];steps=[(x[-1]-x[0])/(len(x)-1) if len(x)>1 else 1 for x in axes]
 result={'physicalPoints':len(points),'shape':list(map(len,axes)),'bounds':[[x[0],x[-1]] for x in axes],'replicaMaxDifference':replica_error,'fields':{}}
 for name in ('rho','p','MomentumX','MomentumY','MomentumZ','TotalEnergyDensity','TurbulenceK','omega','TurbulentViscosity'):
  if name not in next(iter(points.values())):continue
  values=[v[name][0] for v in points.values()];total=0
  for xyz,row in points.items():
   weight=math.prod(h*(.5 if len(axis)>1 and pos in (axis[0],axis[-1]) else 1) for pos,axis,h in zip(xyz,axes,steps));total+=weight*row[name][0]
  result['fields'][name]={'min':min(values),'max':max(values),'physicalTrapezoidalIntegral':total,'L2Nodal':math.sqrt(sum(v*v for v in values)/len(values))}
 return result
results=[]
pairs=[(name,base/'compatibility-closeout'/(name+'-before'),base/'compatibility-closeout'/(name+'-after')) for name in ('sod','weno7','ghost','sstGhost','legacySST')]
pairs.append(('sod-mpi4',base/'compatibility/sod-mpi4-before',base/'compatibility/sod-mpi4-closeout'))
for name,a,b in pairs:
 fa,fb=files(a),files(b);same=len(fa)==len(fb) and bool(fa) and all((b/'result'/f.relative_to(a/'result')).read_bytes()==f.read_bytes() for f in fa);ta,tb=times(a),times(b)
 final=max(float(re.search(r'_t([\deE+.-]+)(?:_block\d+)?\.vts',f.name).group(1)) for f in fb)
 paths=[f for f in fb if float(re.search(r'_t([\deE+.-]+)(?:_block\d+)?\.vts',f.name).group(1))==final]
 result={'name':name,'before':str(a),'after':str(b),'byteIdenticalPhysicalVTS':same,'physicalVTSFiles':len(fb),'dtSequenceIdentical':ta==tb,'actualSteps':len(tb),'finalTime':final,'dtSequence':tb,'finalMetrics':field_metrics(paths)}
 if not same or ta!=tb:raise RuntimeError('Frozen regression failed '+name)
 results.append(result)
native=base/'sstTransport/nativeSST-closeout';ts=times(native);paths=[f for f in files(native) if '_t0.001.' in f.name];m=field_metrics(paths)
if len(ts)!=2 or abs(ts[-1][0]-.001)>1e-12:raise RuntimeError('Native production clock incorrect')
for name in ('rho','p','TurbulenceK','omega','TurbulentViscosity'):
 if not math.isfinite(m['fields'][name]['min']) or m['fields'][name]['min']<=0:raise RuntimeError('Invalid native production field '+name)
results.append({'name':'nativeSST','directory':str(native),'actualSteps':len(ts),'finalTime':ts[-1][0],'dtSequence':ts,'finalMetrics':m})
(base/'evidence/regression-metrics.json').write_text(json.dumps(results,indent=2))
print([(x['name'],x['actualSteps'],x['finalTime'],x.get('byteIdenticalPhysicalVTS','new mathematics')) for x in results])
# Valid binary hashes remain local; ignore macOS AppleDouble metadata in frame counts.
(base/'evidence/final-binary-hashes.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'build/sonicSolver',root/'build/test_stageSST',root/'build/test_sharedFluxTransport',base/'evidence/sonicSolver-before']},indent=2))
