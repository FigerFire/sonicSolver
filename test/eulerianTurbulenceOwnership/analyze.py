"""Quantify the preserved fields, row systems and physical-time trace.

Only reads Phase C evidence inside test/; does not execute or change cases.
Physical integrals use Cartesian nodal dual volumes, including half volumes
at both EMPTY end planes, rather than an unweighted nodal sum.
"""
from pathlib import Path
from collections import Counter, defaultdict
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parent


def arrays(path):
    piece = ET.parse(path).find('.//Piece')
    extent = list(map(int, piece.attrib['Extent'].split()))
    shape = [extent[2*a+1]-extent[2*a]+1 for a in range(3)]
    coordinates = list(map(float, piece.find('Points/DataArray').text.split()))
    axes = [[coordinates[3*(i*(1 if a == 0 else shape[0] if a == 1 else shape[0]*shape[1]))+a]
             for i in range(shape[a])] for a in range(3)]
    widths = []
    for axis in axes:
        widths.append([((axis[min(i+1,len(axis)-1)]-axis[max(i-1,0)])/2)
                       if len(axis)>1 else 1.0 for i in range(len(axis))])
    weights = [widths[0][i]*widths[1][j]*widths[2][k]
               for k in range(shape[2]) for j in range(shape[1]) for i in range(shape[0])]
    values = {a.attrib['Name']: (int(a.attrib.get('NumberOfComponents',1)),list(map(float,a.text.split())))
              for a in piece.findall('PointData/DataArray') if a.attrib['type'].startswith('Float')}
    return shape, weights, values


def statistics(values, weights):
    return {'min':min(values),'max':max(values),
            'L1_mean':math.fsum(abs(v) for v in values)/len(values),
            'L2_rms':math.sqrt(math.fsum(v*v for v in values)/len(values)),
            'Linf':max(abs(v) for v in values),
            'volume_integral':math.fsum(v*w for v,w in zip(values,weights))}


summary = []
for case in sorted((BASE/'final').iterdir()):
    if not case.is_dir() or case.name.startswith('.'): continue
    reference = BASE/'instrumented'/case.name
    log = (case/'run.log').read_text()
    old_log = (reference/'run.log').read_text()
    trace = [l for l in log.splitlines() if l.startswith('[SF PLAN]')]
    old_trace = [l for l in old_log.splitlines() if l.startswith('[SF PLAN]')]
    new_leaves = ('ee.turbulence.boundary','ee.turbulence.validate','ee.turbulence.refresh')
    retained = [l for l in trace if not l.split('op=')[-1] in new_leaves]
    assert retained == old_trace, 'First old-operation divergence: '+case.name
    counts = Counter(l.split('op=')[-1] for l in trace)
    dts = [float(re.search(r'dt=([^ ]+)',l).group(1)) for l in trace if l.endswith('op=ee.time.commit')]
    steps = len(dts)
    frames = []
    for file in sorted((case/'result').glob('*.vts')):
        if file.name.startswith('.'): continue
        shape,weights,fields = arrays(file)
        _,_,before = arrays(reference/'result'/file.name)
        assert fields == before
        metrics = {}
        for name,(components,values) in fields.items():
            for component in range(components):
                v=values[component::components]
                metrics[name if components==1 else name+f'[{component}]'] = statistics(v,weights)
        phases=['water','air']
        mass=math.fsum(metrics['phaseMass.'+p]['volume_integral'] for p in phases)
        enthalpy=math.fsum(metrics['phaseEnthalpy.'+p]['volume_integral'] for p in phases)
        momentum=[math.fsum(w*fields['phaseMass.'+p][1][i]*fields['U.'+p][1][3*i+c]
                    for p in phases for i,w in enumerate(weights)) for c in range(3)]
        frames.append({'file':file.name,'volume':math.fsum(weights),'fields':metrics,
                       'conserved':{'total_phase_mass':mass,'total_phase_momentum':momentum,
                                    'total_phase_enthalpy':enthalpy},
                       'all_before_after_metric_deltas':0.0})
    evidence = case/'numerical.txt'
    diagnostics = evidence.read_text().splitlines() if evidence.exists() else []
    old_values = defaultdict(list)
    for line in diagnostics:
        if line.startswith('old '):
            parts=line.split();old_values[parts[1]].append(parts[2:])
    outer = int(re.search(r'(?m)^  outerCorrectors: (\d+)',(case/'solvers/algorithm.yaml').read_text()).group(1))
    frozen_comparisons = 0
    for variable,snapshots in old_values.items():
        assert len(snapshots)==steps*outer
        for step in range(steps):
            for repeat in range(1,outer):
                assert snapshots[step*outer] == snapshots[step*outer+repeat], 'OldTime advanced in outer'
                frozen_comparisons += 1
    # Physical values and old time must refer to the same phase/storage slot.
    # Compare actual next-step snapshots at every physical node against the
    # preceding committed output, not just a selected-array ordinal.
    commit_checks=0
    for variable,snapshots in old_values.items():
        name=variable.removeprefix('previous.') if hasattr(variable,'removeprefix') else variable[9:]
        phase=name.split('.')[-1]
        for step in range(1,steps):
            shape,_,fields=arrays(case/'result'/frames[step]['file'])
            mass=fields['phaseMass.'+phase][1]
            output_name=name if name in fields else name.split('.')[0]
            assert output_name in fields and (output_name==name or len(old_values)==2)
            phi=fields[output_name][1]
            ng=3;mx=shape[0]+2*ng;my=shape[1]+2*ng
            current=list(map(float,snapshots[step*outer]))
            for k in range(shape[2]):
                for j in range(shape[1]):
                    for i in range(shape[0]):
                        point=(k*shape[1]+j)*shape[0]+i
                        cell=((k+ng)*my+j+ng)*mx+i+ng
                        assert current[cell]==mass[point]*phi[point], 'Committed old-time product differs'
                        commit_checks+=1
    solves = [l.split() for l in diagnostics if l.startswith('solve ')]
    row_count = sum(l.startswith('row ') for l in diagnostics)
    summary.append({'case':case.name,'steps':steps,'final_time':math.fsum(dts),'dt_sequence':dts,
                    'trace_before':len(old_trace),'trace_after':len(trace),'operations':dict(counts),
                    'retained_trace_equal':True,'matrix_rows':row_count,'linear_solves':len(solves),
                    'linear_iterations':sum(int(s[1]) for s in solves),
                    'max_linear_residual':max((float(s[2]) for s in solves),default=0),
                    'frozen_outer_comparisons':frozen_comparisons,'physical_commit_product_checks':commit_checks,
                    'oldtime_changes_between_steps':sum(a!=b for snapshots in old_values.values()
                        for a,b in zip(snapshots[::outer],snapshots[outer::outer])),
                    'matrix_and_solution_bytes_equal':not evidence.exists() or evidence.read_bytes()==(reference/'numerical.txt').read_bytes(),
                    'frames':frames})
(BASE/'evidence/metrics.json').write_text(json.dumps(summary,indent=2))
for item in summary:
    print(item['case'],'steps',item['steps'],'old-time frozen checks',item['frozen_outer_comparisons'],
          'commit products',item['physical_commit_product_checks'],'rows',item['matrix_rows'],
          'solves',item['linear_solves'],'VTS frames',len(item['frames']),'all deltas 0')
