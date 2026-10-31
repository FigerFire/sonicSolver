"""Create isolated native inputs only. Does not run any solver or check."""
from pathlib import Path
import json
import shutil

root = Path(__file__).resolve().parents[2]
base = Path(__file__).resolve().parent

def write(path, obj, kind, **body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'SonicFile': {'object': obj, 'type': kind}, **body}, indent=2)+'\n')

def sym(name): return {'symbol': name}
def number(value, dimensions=None):
    return {'constant': value, **({'dimensions': dimensions} if dimensions is not None else {})}
def op(name, *args, occurrence=None):
    return {'op': name, 'args': list(args), **({'occurrence': occurrence} if occurrence else {})}
def call(equation, target, identity):
    return {'equation': equation, 'target': target, 'id': identity}
def read(identity, version, *symbols):
    return {'call': identity, 'reads': [{'symbol': s, 'version': version} for s in symbols]}

rate = [0, 0, -1, 0, 0, 0, 0]
force = [1, -2, -2, 0, 0, 0, 0]
power = [1, -1, -3, 0, 0, 0, 0]
acceleration = [0, 1, -2, 0, 0, 0, 0]
manifest = []
for mode in ['density', 'pressure']:
    source = root/('test/Sod/sodCase' if mode == 'density' else 'test/pressure/constantDensityPoiseuille')
    for feedback in [False, True]:
        p = base/mode/('feedback' if feedback else 'baseline')
        if p.exists():
            raise RuntimeError(f'Refusing to overwrite an existing input/evidence directory: {p}')
        for folder in ['mesh', 'solvers', 'models', 'fields']:
            (p/folder).mkdir(parents=True)
        # Preserve supplied geometry inputs, never source results or symlinks.
        for name in ['generation.yaml', 'mesh.sfm']:
            path = source/'mesh'/name
            if path.is_file() and not path.is_symlink(): shutil.copy2(path, p/'mesh'/name)
        write(p/'mesh/mesh.yaml', 'mesh', 'registry', files=['mesh/mesh.sfm'], generation={'type': 'structuredMesh', 'file': 'mesh/generation.yaml'})
        write(p/'case.yaml', 'case', 'registry', formatVersion=2, name=f'{mode}-{p.name}')
        roles = [('runtime','runtime'), ('numerics','numerics'), ('algorithm','algorithm'), ('state','stateRegistry'), ('equationRegistry','equationRegistry'), ('execution','algorithmRegistry'), ('how','executionProgram'), ('which','providerBindings')]
        write(p/'solvers/solvers.yaml', 'solver', 'registry', **{k: {'type': v, 'file': f'solvers/{k}.yaml'} for k,v in roles})
        write(p/'solvers/execution.yaml', 'solver', 'execution', **({'Explicit': {}} if mode == 'density' else {'PISO': {'correctors': 2, 'outerCorrectors': 1, 'nonOrthogonalCorrectors': 0}}))
        write(p/'models/models.yaml', 'models', 'registry', thermoDynamics={'type':'thermoDynamics','file':'models/thermoDynamics.yaml'})
        fields = {s: {'type': 'vector' if s=='U' else 'scalar', 'domain':'fluid','location':'point','dimensions': d,'storage':'primary','output':True} for s,d in [('U',[0,1,-1,0,0,0,0]), ('p',[1,-1,-2,0,0,0,0]), ('rho',[1,-3,0,0,0,0,0])]}
        fields.update(internalField={'type':'internalField','file':'fields/internalField.yaml'}, boundaries={'type':'boundary','file':'fields/boundaries.yaml'})
        additions, extensions, state_reads, bindings = {}, {}, [], []
        if mode == 'density':
            write(p/'solvers/runtime.yaml','solver','runtime',startTime=0,endTime=1.0,endStep=10000,CFL=.5,maxDeltaT=.025,createMesh=True,writeControl='timeStep',writeInterval=1,writeInitial=True,output={'jobName':'DensityLoop','outputDir':'result'},parallel={'enabled':False,'nProcs':1})
            write(p/'solvers/numerics.yaml','solver','numerics',time={'default':'classicalRK4'},convection={'default':'TENO5','formulation':'conservativeFluxDifference','reconstruction':'characteristic','flux':'StegerWarming','interfaceFlux':'shared'})
            write(p/'solvers/algorithm.yaml','solver','algorithm',type='densityBase',algorithm='Explicit')
            write(p/'models/thermoDynamics.yaml','models','thermoDynamics',thermoDynamics={'equationOfState':'perfectGas','thermo':'hConst','transport':'const'},properties={'equationOfState':{'gamma':1.4,'R':287.05},'transport':{'mu':0.0,'Pr':.72}})
            write(p/'fields/internalField.yaml','fields','internalField',U={'uniform':[.25,0,0]},p={'uniform':1.0},rho={'uniform':1.2})
            sides = ['left','right','bottom','top','front','back']
            # The copy uses six named mesh patches rather than the Sod IC sets.
            g = json.loads((p/'mesh/generation.yaml').read_text())
            g['patches'] = [{'name': n, 'faces': [f]} for n,f in zip(sides, [[0,4,7,3],[1,2,6,5],[0,1,5,4],[3,7,6,2],[0,3,2,1],[4,5,6,7]])]
            (p/'mesh/generation.yaml').write_text(json.dumps(g,indent=2)+'\n')
            write(p/'fields/boundaries.yaml','fields','boundaries', **{s:{side: 'empty' if side in ['front','back'] else 'zeroGradient' for side in sides} for s in ['U','p','rho']})
            equations = ['Continuity','Momentum','Energy']; solution = ['rho','rhoU','rhoE']
            body = []
            if feedback:
                additions = {'drag': {'lhs': sym('F'), 'rhs': op('multiply',number(-.2,rate),op('component',sym('rhoU'),number(0)))}, 'power': {'lhs':sym('P'),'rhs':op('divide',op('multiply',sym('F'),op('component',sym('rhoU'),number(0))),sym('rho'))}}
                for name,d in [('F',force),('P',power)]: fields[name]={'new':True,'type':'scalar','domain':'fluid','location':'point','dimensions':d,'storage':'derived','output':False}
                extensions = {'momentum': {'rhs': op('source',op('vector',sym('F'),number(0,force),number(0,force)),occurrence='drag.source')}, 'energy': {'rhs':op('source',sym('P'),occurrence='work.source')}}
                body += [call('drag','F','drag.update'),call('power','P','work.update')]
                state_reads += [read('drag.update','Stage','rhoU'),read('work.update','Stage','F','rhoU','rho'),read('momentum.rhs','Stage','F'),read('energy.rhs','Stage','P')]
            body += [call(e,t,i) for e,t,i in [('continuity','rho','mass.rhs'),('momentum','rhoU','momentum.rhs'),('energy','rhoE','energy.rhs')]]
            how = [{'kind':'Sequence','id':'explicit.stage','frequency':'Stage','children':body},{'kind':'Commit','id':'physical.accept'}]
            bindings += [{'equation':e,'method':'ConservativeResidual'} for e in ['continuity','momentum','energy']]
        else:
            shutil.copy2(source/'solvers/algorithm.yaml',p/'solvers/algorithm.yaml')
            shutil.copy2(source/'solvers/runtime.yaml',p/'solvers/runtime.yaml')
            shutil.copy2(source/'fields/internalField.yaml',p/'fields/internalField.yaml')
            shutil.copy2(source/'fields/boundaries.yaml',p/'fields/boundaries.yaml')
            write(p/'solvers/numerics.yaml','solver','numerics',time={'default':'forwardEuler'},terms={'convection':'primitiveUpwind1','diffusion':'central2Explicit'})
            write(p/'models/thermoDynamics.yaml','models','thermoDynamics',thermoDynamics={'equationOfState':'rhoConst','transport':'const'},properties={'equationOfState':{'rho':1.0},'transport':{'mu':.1,'Pr':.72}})
            equations = ['Momentum','Continuity']; solution = ['U','p']; how=[]
            if feedback:
                fields['acceleration']={'new':True,'type':'scalar','domain':'fluid','location':'point','dimensions':acceleration,'storage':'derived','output':False}
                additions={'damping':{'lhs':sym('acceleration'),'rhs':op('multiply',number(-.2,rate),op('component',sym('U'),number(0)))}}
                extensions={'momentum':{'rhs':op('source',op('vector',sym('acceleration'),number(0,acceleration),number(0,acceleration)),occurrence='damping.source')}}
                how.append({**call('damping','acceleration','damping.update'),'frequency':'Step'})
                state_reads += [read('damping.update','Current','U'),read('predictor','Current','acceleration')]
            how += [call('momentum','U','predictor'),{'kind':'Loop','id':'pressure','repetitions':2,'children':[{'kind':'Loop','id':'nonOrthogonal','repetitions':1,'children':[{**call('pSimple','deltaP','pressure.solve'),'targetKind':'Correction'}]},call('correctU','U','velocity.correct'),call('correctP','p','pressure.correct'),call('correctFluxp','phi','flux.correct')]},{'kind':'Commit','id':'physical.accept'}]
            bindings += [{'equation':e,'method':m} for e,m in [('momentum','PressureMomentum'),('pSimple','PressureCorrection'),('correctU','VelocityCorrection'),('correctP','PressureUpdate'),('correctFluxp','FluxCorrection')]]
        for equation in additions: bindings.append({'equation':equation,'method':'LocalAffine'})
        write(p/'fields/fields.yaml','fields','registry',**fields)
        write(p/'solvers/state.yaml','solver','state',use=solution,calls=state_reads)
        write(p/'solvers/equationRegistry.yaml','solver','equationRegistry',use=equations+list(additions),**({'add': additions, 'extend':extensions} if additions else {}))
        write(p/'solvers/how.yaml','solver','how',root={'kind':'Sequence','id':'HOW','children':how})
        write(p/'solvers/which.yaml','solver','which',bindings=bindings)
        manifest.append({'path':str(p.relative_to(root)),'source':str(source.relative_to(root)),'preparedOnly':True,'ran':False,'command':[str(root/'build/sonicSolver'),'run',str(p)]})
(base/'evidence/run-plan.json').write_text(json.dumps(manifest,indent=2)+'\n')
