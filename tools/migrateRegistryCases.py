#!/usr/bin/env python3
"""Losslessly migrate supported Sonic OpenFOAM case dictionaries to registry YAML (JSON subset).
Generated results are never touched. Run with --apply after reviewing --check output.
"""
import argparse,json,re,shutil
from pathlib import Path

def tokens(text):
    return re.findall(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/|[{}();\[\]]|[^\s{}();\[\]]+',text)
def parse(path):
    ts=[x for x in tokens(path.read_text()) if not x.startswith(('//','/*'))];i=0
    def atom(t):
        if t.startswith(('"',"'")):return t[1:-1]
        if t in ('true','false'):return t=='true'
        try:return int(t)
        except ValueError:pass
        try:return float(t)
        except ValueError:return t
    def value():
        nonlocal i
        t=ts[i];i+=1
        if t in ('(','['):
            end=')' if t=='(' else ']';a=[]
            while i<len(ts) and ts[i]!=end:a.append(value())
            if i==len(ts):raise ValueError(f'{path}: unclosed list')
            i+=1;return a
        if t=='{':return mapping('}')
        return atom(t)
    def mapping(end=None):
        nonlocal i
        out={}
        while i<len(ts) and ts[i]!=end:
            if ts[i]==';':i+=1;continue
            key=str(atom(ts[i]));i+=1
            if key in out:raise ValueError(f'{path}: duplicate key {key}')
            if ts[i]=='{':i+=1;out[key]=mapping('}');continue
            a=[]
            while i<len(ts) and ts[i]!=';':a.append(value())
            if i==len(ts):raise ValueError(f'{path}: missing semicolon {key}')
            i+=1
            if len(a)==2 and a[0]=='uniform':a=a[1:]
            out[key]=a[0] if len(a)==1 else a
        if end:
            if i==len(ts):raise ValueError(f'{path}: missing {end}')
            i+=1
        return out
    result=mapping();result.pop('FoamFile',None);return result

def write(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def migrate(case,apply):
    control=parse(case/'system/controlDict');control.pop('application');solution=parse(case/'system/fvSolution') if (case/'system/fvSolution').exists() else {}
    control['output']=parse(case/'system/outputDict');control['parallel']=solution.get('parallel',{})
    # Keep output locations and numerical settings unchanged.
    schemes=parse(case/'system/fvSchemes') if (case/'system/fvSchemes').exists() else {}
    names={'ddtSchemes':'time','divSchemes':'convection','laplacianSchemes':'diffusion','sourceSchemes':'sources'}
    aliases={'timeScheme':('time','default'),'convectionScheme':('convection','default'),'fluxSplitter':('convection','flux'),'viscousScheme':('diffusion','default'),'sourceScheme':('sources','default')}
    flat={k:schemes.pop(k) for k in list(schemes) if k in aliases}
    # These old spellings were not consumed by the old reader. Record them for review.
    ignored={k:schemes.pop(k) for k in ('convection','flux','viscous') if k in schemes}
    unknown=set(schemes)-set(names)
    if unknown:raise ValueError(f'{case}: unsupported legacy scheme keys {unknown}')
    numerics={names[k]:v for k,v in schemes.items()};numerics['transport']=solution.get('viscous',{})
    for key,value in flat.items():
        group,option=aliases[key];numerics.setdefault(group,{}).setdefault(option,value)
    if 'pressureBased' in solution:numerics['pressureCorrection']=solution['pressureBased']
    if 'solverFormulation' in solution:ignored['fvSolution.solverFormulation']=solution['solverFormulation']
    if set(solution)-{'parallel','viscous','pressureBased','solverFormulation'}:raise ValueError(f'{case}: unsupported fvSolution keys')
    fields={};boundaries={}
    start=str(control.get('startTime',0));start=str(int(float(start))) if float(start).is_integer() else start
    for p in sorted((case/start).glob('*')):
        if p.name.startswith('.') or not p.is_file():continue
        data=parse(p);initial=data.pop('internalField',None)
        typ='vector' if isinstance(initial,list) and len(initial)==3 else 'scalar'
        if initial is None:raise ValueError(f'{p}: no internalField')
        boundaries[p.name]=data.pop('boundaryField',{})
        f={'type':typ,'domain':'fluid','location':'point','dimensions':data.pop('dimensions','unspecified'),
           'storage':'primary','initial':initial,'sets':data.pop('internalSets',{}),'output':True,'restart':False}
        aliases={'rho':'Density','U':'Velocity','p':'Pressure','T':'Temperature'}
        f['outputName']=aliases.get(p.name,p.name)
        # Checkpoint ownership is bound separately; primitive initial inputs are not canonical restart state.
        if data:raise ValueError(f'{p}: unmigrated field properties {data.keys()}')
        fields[p.name]=f
    closures={};constraints={};geometry={};moves=[];configfiles=[]
    types={'thermophysicalProperties':'Thermophysical','phaseProperties':'PhaseSystem','phaseChange':'PhaseChange',
           'turbulenceProperties':'Turbulence','IBMProperties':'ImmersedBoundary','ILWProperties':'ILW',
           'g':'Gravity','MRFProperties':'MRF','wallHeatSourceProperties':'WallHeat'}
    for name,typ in types.items():
        p=case/'constant'/name
        if not p.exists():continue
        params=parse(p)
        if name=='IBMProperties':
            block=params.get('IBM',params)
            for key in ('geometryFiles','files'):
                if key in block:
                    values=block[key];values=values if isinstance(values,list) else [values]
                    dest=[]
                    for v in values:
                        source=case/v if '/' in v else case/'constant/triSurface'/v
                        if not source.exists():raise ValueError(f'{case}: missing geometry {source}')
                        target=case/'geometry'/source.name
                        moves.append((source,target));dest.append('geometry/'+source.name)
                        geometry[source.stem]={'type':'File','parameters':{'file':dest[-1]}}
                    block[key]=dest
        target=constraints if typ in ('ImmersedBoundary','ILW') else closures
        target[name]={'type':typ,'parameters':params};configfiles.append(p)
    phase = closures.get('phaseProperties',{}).get('parameters',{})
    phase_kind = phase.get('phaseSystem')
    if phase_kind is None:
        phase_kind = phase.get('multiPhase',{}).get('type')
    turbulence = closures.get('turbulenceProperties',{}).get('parameters',{})
    turbulence_output = (turbulence.get('enabled',False)
                         and str(turbulence.get('simulationType','')).lower() == 'ras')
    turbulence_model = str(turbulence.get('RAS',{}).get('model','')).lower()
    turbulence_fields = {'k'}
    if turbulence_model == 'kepsilon': turbulence_fields.add('epsilon')
    if turbulence_model == 'komegasst': turbulence_fields.add('omega')
    for name,field in fields.items():
        field['output'] = (
            name in ('rho','U','p')
            or (name == 'T' and phase_kind != 'eulerianEulerian')
            or (name == 'phi' and str(phase_kind).lower() == 'levelset')
            or (name in turbulence_fields and turbulence_output)
            or ('.' in name and phase_kind == 'eulerianEulerian'))
    meshfiles=['mesh/mesh.sfm']
    if (case/'constant/mesh_sets.sfm').exists():meshfiles.append('mesh/mesh_sets.sfm')
    for p in (case/'constant').glob('*'):
        if (p.is_file() and not p.name.startswith('.')
                and (p.suffix.lower() in ('.sfm','.txt','.vts','.vtm'))):
            moves.append((p,case/'mesh'/p.name))
    generator=case/'system/blockMeshDict'
    if generator.exists():moves.append((generator,case/'mesh/blockMeshDict'))
    preset='eulerianEulerian' if any(o['parameters'].get('phaseSystem')=='eulerianEulerian' for o in closures.values()) else 'compressible'
    documents={'case.yaml':{'formatVersion':1,'name':case.name},'runtime.yaml':control,'numerics.yaml':numerics,
        'solver.yaml':parse(case/'system/solverProperties') if (case/'system/solverProperties').exists() else {'type':'densityBase'},
        'mesh/mesh.yaml':{'files':meshfiles,'generator':'mesh/blockMeshDict'},
        'model/fields.yaml':{'fields':fields},'model/boundaries.yaml':{'boundaries':boundaries},
        'model/closures.yaml':{'closures':closures},'model/constraints.yaml':{'constraints':constraints},
        'model/geometry.yaml':{'geometry':geometry},'model/equations.yaml':{'equations':{'flow':{'type':'BuiltinEquationSystem','parameters':{'preset':preset}}}}}
    if apply:
        if (case/'case.yaml').exists():raise ValueError(f'{case}: already migrated, refusing overwrite')
        for rel,obj in documents.items():write(case/rel,obj)
        for src,dst in dict(moves).items():
            if src==dst:continue
            dst.parent.mkdir(parents=True,exist_ok=True)
            if dst.exists():raise ValueError(f'{dst} already exists')
            shutil.move(str(src),str(dst))
        # Old inputs are backed up by the task snapshot; outputs and other assets stay in place.
        for p in configfiles:p.unlink()
        for p in (case/'system').glob('*'):
            if p.is_file() and not p.name.startswith('.'):
                if p.name=='sonicDict' and parse(p):raise ValueError(f'{p}: nonempty legacy extension')
                p.unlink()
        for p in (case/start).glob('*'):
            if p.is_file() and not p.name.startswith('.'):p.unlink()
        for legacy in (case/start,case/'system',case/'constant'):
            if legacy.exists() and not any(
                    item for item in legacy.rglob('*')
                    if not item.name.startswith('._')
                    and (item.is_file() or item.is_symlink())):
                shutil.rmtree(legacy)
    return {'case':str(case),'fields':len(fields),'objects':len(closures)+len(constraints),'files':len(documents),'ignoredLegacySchemeKeys':ignored}

def main():
    p=argparse.ArgumentParser();p.add_argument('root',nargs='?',default='test');p.add_argument('--apply',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    cases=sorted(x.parent.parent for x in Path(a.root).rglob('controlDict') if x.parent.name=='system' and not (x.parent.parent/'case.yaml').exists())
    results=[migrate(c,False) for c in cases] # complete preflight before mutations
    if a.apply:results=[migrate(c,True) for c in cases]
    print(json.dumps(results,indent=2))
if __name__=='__main__':main()
