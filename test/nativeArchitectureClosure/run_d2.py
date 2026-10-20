from pathlib import Path
import sys,shutil,json
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B.parents[1]/'test/levelSetMpiD1'));import cases
cases.BASE=B
phase=sys.argv[1];binary=Path(sys.argv[2]).resolve();cases.generate(phase)
for original,graded in [('planeAdvection','gradedPlaneAdvection'),('circleGeometry','gradedCircle'),('csf','gradedCsf')]:
 for label in ('serial','mpi4xy'):
  p=B/phase/graded/label
  if p.exists():continue
  shutil.copytree(B/phase/original/label,p,ignore=shutil.ignore_patterns('result','checkpoint','*.log','._*'))
  f=p/'mesh/mesh.sfm';s=f.read_text();pre,raw=s.split('#Point\n',1);points,post=raw.split('#Cell',1);updated=[]
  for line in points.splitlines():
   if not line.strip():continue
   x,y,z=map(float,line.split());updated.append('%.17g %.17g %.17g'%(100*(x/100)**1.2,20*(y/20)**1.1,z))
  f.write_text(pre+'#Point\n'+'\n'.join(updated)+'\n#Cell'+post)
p=B/phase/'serialReinitialization'/'serial'
if not p.exists():
 shutil.copytree(B/phase/'circleGeometry/serial',p,ignore=shutil.ignore_patterns('result','checkpoint','*.log','._*'))
 f=p/'models/multiPhase.yaml';f.write_text(f.read_text().replace('reinitializationSteps: 0','reinitializationSteps: 3'))
cases.run(phase,binary,[])
