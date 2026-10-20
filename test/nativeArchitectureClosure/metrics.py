from pathlib import Path
import xml.etree.ElementTree as ET,json,re,math
B=Path(__file__).resolve().parent;out=[]
records=json.loads((B/'evidence/protection-manifest.json').read_text())+json.loads((B/'evidence/mpi-protection-manifest.json').read_text())
for record in records:
 if record.get('stage',record.get('version'))!='after' or record['exit']!=0:continue
 name=record.get('case',record.get('name'));parent='protection' if 'case' in record else 'mpiProtection';p=B/parent/(name+'-after');before=p.with_name(name+'-before')
 def frames(case):
  all={}
  for f in sorted((case/'result').rglob('*')):
   if f.name.startswith('.') or f.suffix not in ('.vts','.vtu') or '_patches' in f.name:continue
   key=re.search(r'_(\d{6})(?:_|\.)',f.name)
   if not key:continue
   piece=ET.parse(f).find('.//Piece');point=piece.find('./Points/DataArray')
   if point is None or point.attrib.get('format','ascii')!='ascii':raise AssertionError(('non-ascii',f))
   xyz=list(map(float,point.text.split()));data={a.attrib['Name']:(int(a.attrib.get('NumberOfComponents',1)),list(map(float,a.text.split()))) for a in piece.findall('./PointData/DataArray')}
   group=all.setdefault(key.group(1),{})
   for n in range(len(xyz)//3):
    coord=tuple(round(v,12) for v in xyz[3*n:3*n+3]);values={k:arr[c*n:c*(n+1)] for k,(c,arr) in data.items() if not k.startswith('BC_') and k not in ('vtkGhostType','GlobalPointId')}
    if coord in group:
     if group[coord]!=values:raise AssertionError(('replica data conflict',f,coord))
    else:group[coord]=values
  return all
 a,c=frames(before),frames(p);assert a==c,(name,'numerical first difference')
 for frame,rows in c.items():
  fields=next(iter(rows.values())).keys();stats={}
  axis=[sorted(set(pt[i] for pt in rows)) for i in range(3)];cartesian=math.prod(len(s) for s in axis)==len(rows) and all(len(s)>1 for s in axis)
  weights=[]
  for coords in axis:
   weights.append({x:((coords[i+1]-x)/2 if i==0 else (x-coords[i-1])/2 if i==len(coords)-1 else (coords[i+1]-coords[i-1])/2) for i,x in enumerate(coords)} if len(coords)>1 else {})
  for field in fields:
   vals=[v for row in rows.values() for v in row[field]]
   assert all(math.isfinite(v) for v in vals),(name,field)
   stats[field]={'min':min(vals),'max':max(vals),'nodalL1':sum(abs(v) for v in vals),'nodalL2':math.sqrt(sum(v*v for v in vals))}
   if cartesian:stats[field]['wholeDomainTrapezoidalIntegral']=[sum(row[field][i]*math.prod(weights[d][coord[d]] for d in range(3)) for coord,row in rows.items()) for i in range(len(next(iter(rows.values()))[field]))]
  out.append(dict(case=name,frame=frame,uniquePhysicalPoints=len(rows),exact=True,stats=stats,wholeDomainIntegralAvailable=cartesian))
(B/'evidence/numerical-metrics.json').write_text(json.dumps(out,indent=2))
print(len(out),'frames: exact replica-aware fields and all derived metrics; whole-domain Cartesian quadrature only, not a body-subtracted conservation proof')
