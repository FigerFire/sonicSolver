"""Write actual short-case settings and executable identities for the report."""
from pathlib import Path
import hashlib
import json
import re
import xml.etree.ElementTree as ET

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def registered_inputs(case):
    queue=['case.yaml','solvers/solvers.yaml','models/models.yaml','fields/fields.yaml','mesh/mesh.yaml']
    seen=set()
    while queue:
        name=queue.pop(0)
        if name in seen:continue
        seen.add(name);file=case/name
        if not file.exists():continue
        contents=file.read_text()
        yield name,contents
        queue.extend(re.findall(r'(?m)^\s*file:\s*([^\s#]+\.yaml)',contents))


cases=[]
for item in json.loads((BASE/'evidence/final-manifest.json').read_text()):
    cases.append(('Eulerian '+item['case'],BASE/'final'/item['case'],
                  'test/turbulence/euleriankEpsilon' if item['case'].startswith('kEpsilon') else 'test/turbulence/euleriankOmegaSST'))
for item in json.loads((BASE/'evidence/protection-manifest.json').read_text()):
    if item['stage']=='after':cases.append(('保护 '+item['case'],BASE/'protection'/(item['case']+'-after'),item['source']))

lines=['# Phase C 实际算例设置与证据索引',
       '', '本文由 `settings.py` 从最终实际运行副本生成。所有 CFD 均为 serial、单块、startTime=0、endTime≤0.001 s。',
       '完整输入仅展示母 registry 及其显式注册的 YAML 子文件；不将目录内未注册文件宣称为已启用。',
       '网格范围来自数值结果中的真实物理坐标；采用输入的 SI 参数解释，求解器不自动判定 mm/m。',
       '原正式 case 未改写。每个对照只缩短时间、逐步输出；pressure/legacySST 另补显式原物性参数，pressure 先生成网格再求解。',
       'Eulerian variants 的额外差异见中文主报告；before/after 数学输入相同。未列值为实现缺省/继承，而非新设定。','']
manifest=[]
for title,case,source in cases:
    frames=sorted(f for f in (case/'result').glob('*.vts') if not f.name.startswith('.'))
    assert frames
    xml=ET.parse(frames[-1]);piece=xml.find('.//Piece')
    ext=list(map(int,piece.attrib['Extent'].split()));shape=[ext[2*a+1]-ext[2*a]+1 for a in range(3)]
    xyz=list(map(float,piece.find('Points/DataArray').text.split()))
    bounds=[(min(xyz[a::3]),max(xyz[a::3])) for a in range(3)]
    log=(case/'run.log').read_text()
    trace=[l for l in log.splitlines() if l.startswith('[SF PLAN]')]
    commits=[l for l in trace if l.endswith('op=ee.time.commit') or l.endswith('op=time.commit')]
    timeline=[(int(re.search(r'step=(\d+)',l).group(1)),float(re.search(r'time=([^ ]+)',l).group(1)),float(re.search(r'dt=([^ ]+)',l).group(1))) for l in commits]
    final=timeline[-1][1]+timeline[-1][2]
    contents=dict(registered_inputs(case))
    lines+=['## '+title,'',f'- 来源：`{source}`；实际运行：`{case.relative_to(ROOT)}`。',
            f'- 网格物理范围 x/y/z：{bounds}；节点数 {shape}；单元数 {[n-1 for n in shape]}；单 patch，serial 1 rank。',
            f'- 实际 {len(commits)} 个物理步，最终 time={final:.17g}；{len(frames)} 个 VTS（含初值）。',
            '- dt 序列：`'+', '.join(f'{x[2]:.17g}' for x in timeline)+'`。',
            f'- 命令：`build/sonicSolver run {case.relative_to(ROOT)}`；日志 `run.log`、数值输出 `result/`。before 使用冻结的 `evidence/sonicSolver-before`。',
            '- 验证：同配置新旧逐帧输出、操作/time/dt trace；Eulerian 另比对全矩阵和 old-time。短时 preservation 不代表长期物理验证。','']
    if 'WHAT / EQUATIONS' in log:
        section=log.split('WHAT / EQUATIONS',1)[1].split('\n\n',1)[0].strip()
        lines+=['实际编译 WHAT（不是根据目录名推测）：','', '```text',section,'```','']
    for name,text in contents.items():
        lines += ['### '+name,'','```yaml',text.rstrip(),'```','']
    generator=case/'mesh/blockMeshDict'
    if generator.exists():lines+=['### 背景网格生成输入（数值网格已由上面的 VTS 核验）','','```text',generator.read_text().rstrip(),'```','']
    manifest.append({'name':title,'source':source,'case':str(case.relative_to(ROOT)),
                     'shape':shape,'bounds':bounds,'steps':len(commits),'final_time':final,'frames':len(frames),
                     'dt_sequence':[x[2] for x in timeline],
                     'inputs':{name:hashlib.sha256(text.encode()).hexdigest() for name,text in contents.items()}})
(BASE/'evidence/case-settings.md').write_text('\n'.join(lines)+'\n')
(BASE/'evidence/case-settings.json').write_text(json.dumps(manifest,indent=2))
(BASE/'evidence/executable-identities.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in
    [BASE/'evidence/sonicSolver-before',BASE/'evidence/sonicSolver-instrumented',ROOT/'build/sonicSolver']},indent=2))
for entry in manifest:
    print(entry['name'],entry['shape'],entry['bounds'],entry['steps'],entry['final_time'])
