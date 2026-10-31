"""Compare A9's freshly executed baseline and candidate; never use A8 output as truth."""
from pathlib import Path
import hashlib, json, math, xml.etree.ElementTree as E

base = Path(__file__).resolve().parent
root = base.parents[1]
def data(path):
    tree=E.parse(path);extent=list(map(int,tree.find('.//StructuredGrid').get('WholeExtent').split()))
    n=[extent[2*i+1]-extent[2*i]+1 for i in range(3)]
    points=list(map(float,tree.find('.//Points/DataArray').text.split()));coords=[points[i::3] for i in range(3)]
    # Structured indices preserve an axis even when coordinate serialization differs by roundoff.
    axes=[coords[0][:n[0]],coords[1][::n[0]][:n[1]],coords[2][::n[0]*n[1]]]
    weights=[]
    for axis in axes:
        assert len(axis)>1 and all(b>a for a,b in zip(axis,axis[1:])), 'Non-monotone diagnostic axis'
        weights.append([.5*(axis[min(i+1,len(axis)-1)]-axis[max(i-1,0)]) for i in range(len(axis))])
    volume=math.prod(a[-1]-a[0] for a in axes);cellWeights=[]
    for k in range(n[2]):
        for j in range(n[1]):
            for i in range(n[0]):
                idx=(k*n[1]+j)*n[0]+i
                assert all(abs(coords[c][idx]-axes[c][v])<1e-12*max(1.,abs(axes[c][-1]-axes[c][0])) for c,v in enumerate((i,j,k))), 'Non-rectilinear diagnostic mesh'
                cellWeights.append(weights[0][i]*weights[1][j]*weights[2][k])
    assert abs(math.fsum(cellWeights)-volume)<1e-13*volume
    fields={a.get('Name'):(int(a.get('NumberOfComponents','1')),list(map(float,a.text.split()))) for a in tree.findall('.//PointData/DataArray') if not a.get('Name').startswith('BC_')}
    integral=lambda values:math.fsum(w*v for w,v in zip(cellWeights,values))
    statistics={key:dict(min=min(values),max=max(values),weightedL2=math.sqrt(math.fsum(w*math.fsum(v*v for v in values[i*c:(i+1)*c]) for i,w in enumerate(cellWeights))/volume)) for key,(c,values) in fields.items()}
    totals={}
    if 'rho' in fields:totals['mass']=integral(fields['rho'][1])
    for key in ('MomentumX','MomentumY','MomentumZ','TotalEnergyDensity'):
        if key in fields:totals[key]=integral(fields[key][1])
    if 'rho' in fields and 'MomentumX' not in fields and 'U' in fields:
        totals['momentumFromRhoU']=[integral([r*u for r,u in zip(fields['rho'][1],fields['U'][1][c::3])]) for c in range(3)]
    for key,(_,values) in fields.items():
        if key.startswith('phaseMass.'):
            phase=key.split('.',1)[1];totals[key]=integral(values)
            totals['phaseMomentum.'+phase]=[integral([r*u for r,u in zip(values,fields['U.'+phase][1][c::3])]) for c in range(3)]
    return dict(nodes=n,domain=[[a[0],a[-1]] for a in axes],quadrature='structured-index tensor-product trapezoidal physical-volume weights',totals=totals,fields=statistics)

def payload(case):
    return {str(f.relative_to(case / 'result')): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in (case / 'result').rglob('*')
            if f.is_file() and not f.name.startswith('._') and f.suffix not in ('.pvd', '.txt')}

def trace(case):
    return [line for line in (case / 'run.log').read_text().splitlines()
            if line.startswith(('[SF PLAN]', '[SF] Step time:', '[SF] Scalar step:', '[SF SCALAR stage]'))]

records = []
for group, names in [('', ['absent', 'laminar', 'gravity', 'kepsilon', 'kepsilonGravity', 'disabled', 'override', 'sstOnce']),
                     ('protection', ['flowPair', 'scalarPair', 'stageSST', 'eulerianSST'])]:
    for name in names:
        before, after = base / group / 'before' / name, base / group / 'after' / name
        assert payload(before) == payload(after), (group, name, 'numerical payload first difference')
        assert trace(before) == trace(after), (group, name, 'operation/clock first difference')
        frames = E.parse(next((after / 'result').glob('*.pvd'))).findall('.//DataSet')
        last = after / 'result' / frames[-1].get('file')
        stats = data(last)
        assert stats == data(before / 'result' / last.relative_to(after / 'result'))
        records.append(dict(group=group or 'module', case=name, frames=len(frames),
                            finalTime=frames[-1].get('timestep'), statistics=stats,
                            traceLines=len(trace(after)), payloadFiles=len(payload(after))))

def vts(folder):
    return {str(f.relative_to(folder)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in folder.rglob('*.vts') if not f.name.startswith('._')}

groups = []
for before, after, expected in [('targetsBefore/cases', 'targetsAfter/cases', 11),
                                ('before/blocks/cases', 'after/blocks/cases', 13),
                                ('before/boundaries/cases', 'after/boundaries/cases', 9)]:
    left, right = vts(base / before), vts(base / after)
    assert left and left == right, (before, 'VTS payload first difference')
    assert len({Path(name).parts[0] for name in left}) == expected
    groups.append(dict(before=before, after=after, cases=expected, files=len(left), exact=True))

start = json.loads((base / 'evidence/start.json').read_text())
for name, digest in start['userDocsSha256'].items():
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest, ('preexisting user document changed', name)

(base / 'evidence/regression.json').write_text(json.dumps(dict(cases=records, nativeGroups=groups,
    originalUserDocsPreserved=True), indent=2))
print('12 CFD payloads, operation/clock traces and volume statistics; 11 targets + 13 blocks + 9 boundary cases match freshly executed A9 baseline.')
