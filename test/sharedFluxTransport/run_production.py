from pathlib import Path
import shutil,re,subprocess,json,hashlib,os
root=Path(__file__).resolve().parents[2];base=root/'test/sharedFluxTransport';exe=root/'build/sonicSolver';before=base/'evidence/sonicSolver-before';specs=[]
ignore=shutil.ignore_patterns('result','results','*.log','._*','.DS_Store','__pycache__')
# Frozen before/after MPI equivalence uses the same physical input and partition.
for version,binary in [('before',before),('after',exe)]:
 dest=base/'compatibility'/('sod-mpi4-'+version)
 if dest.exists():raise RuntimeError('Refuse to overwrite '+str(dest))
 shutil.copytree(base/'compatibility/sod-before',dest,ignore=ignore)
 p=dest/'solvers/runtime.yaml';s=p.read_text().replace('enabled: false','enabled: true').replace('nProcs: 1','nProcs: 4').replace('split: [1, 1, 1]','split: [4, 1, 1]');p.write_text(s)
 specs.append({'name':'sod-mpi4-'+version,'directory':str(dest),'source':'compatibility/sod-before','modified':['parallel enabled, nProcs=4 split=[4,1,1]'],'command':['mpirun','-np','4',str(binary),'run',str(dest)],'hash':hashlib.sha256(binary.read_bytes()).hexdigest(),'endTime':.01})
case=base/'sstTransport/nativeSST';case.parent.mkdir(parents=True,exist_ok=True)
if case.exists():raise RuntimeError('Refuse to overwrite '+str(case))
shutil.copytree(base/'compatibility/legacySST-before',case,ignore=ignore)
p=case/'models/turbulence.yaml';p.write_text(p.read_text().replace('enabled: true','enabled: true\ntransportMethod: stageConservative'))
p=case/'solvers/runtime.yaml';s=re.sub(r'^endTime:.*$', 'endTime: 0.001',p.read_text(),flags=re.M);s=re.sub(r'^writeInterval:.*$','writeInterval: 0.0005',s,flags=re.M);p.write_text(s)
p=case/'fields/internalField.yaml';p.write_text('SonicFile:\n  object: fields\n  type: internalField\nrho:\n  uniform: 1.2\nU:\n  uniform: [0.3, 0, 0]\np:\n  uniform: 100000\nk:\n  uniform: 0.01\nomega:\n  uniform: 1.0\nphi:\n  uniform: 0\nepsilon:\n  uniform: 0\n')
p=case/'fields/boundaries.yaml';s=p.read_text().replace('k:\n','k:\n  Left: zeroGradient\n').replace('omega:\n','omega:\n  Left: zeroGradient\n');s=s.replace('  Left:','  nativeLeft:').replace('  Right:','  nativeRight:').replace('  sideWalls:','  nativeSides:');p.write_text(s)
p=case/'mesh/mesh.sfm';s=p.read_text();nx,ny,nz=101,21,2
sets={'nativeLeft':[i+nx*(j+ny*k) for k in range(nz) for j in range(ny) for i in [0]],'nativeRight':[i+nx*(j+ny*k) for k in range(nz) for j in range(ny) for i in [nx-1]],'nativeSides':[i+nx*(j+ny*k) for k in range(nz) for j in [0,ny-1] for i in range(nx)]}
for name,ids in sets.items():s+='\n'+name+' '+str(len(ids))+'\n'+' '.join(map(str,ids))+'\n'
p.write_text(s)
specs.append({'name':'nativeSST','directory':str(case),'source':'test/turbulence/singleFluidSST','modified':['explicit transportMethod=stageConservative','uniform initial rho=1.2,U=(.3,0,0),p=100000,k=.01,omega=1','explicit physical face point sets nativeLeft/nativeRight/nativeSides; zeroGradient k/omega and existing flow face tags switched; frontBack empty retained','endTime=.001, writeInterval=.0005'],'command':[str(exe),'run',str(case)],'hash':hashlib.sha256(exe.read_bytes()).hexdigest(),'endTime':.001})
(base/'evidence/production-manifest.json').write_text(json.dumps(specs,indent=2));print('prepared',len(specs),'cases')
