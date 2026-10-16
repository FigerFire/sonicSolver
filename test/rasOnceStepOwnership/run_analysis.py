# -*- coding: utf-8 -*-
"""Compare frozen pre-migration output without tolerances; retain real case settings."""
from pathlib import Path
from collections import Counter,defaultdict
import hashlib,json,math,re,subprocess,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent;root=base.parents[1];out=base/'evidence'
before=json.loads((out/'cases-beforeFinal.json').read_text());after=json.loads((out/'cases-final.json').read_text())
assert len(before)==len(after)==8 and len({x['name'] for x in after})==8, 'Incomplete CFD manifest'
assert all(x['binary']==hashlib.sha256((root/'build/sonicSolver').read_bytes()).hexdigest() for x in after), 'Final executable changed since CFD run'
assert (out/'sweep-before.txt').read_bytes()==(out/'sweep-after.txt').read_bytes(),'First difference in kernel array trace'
summary=[];settings=[]
for a,b in zip(before,after):
 assert a['name']==b['name'] and a['exit']==b['exit']==0
 p,q=Path(a['case']),Path(b['case'])
 for f in p.rglob('*'):
  if f.is_file() and f.suffix in ('.yaml','.sfm','.STL') and 'result' not in f.parts:
   assert f.read_bytes()==(q/f.relative_to(p)).read_bytes(),('input changed',f)
 files=sorted(f for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('.'))
 assert [str(f.relative_to(p/'result')) for f in files]==sorted(str(f.relative_to(q/'result')) for f in (q/'result').rglob('*') if f.is_file() and not f.name.startswith('.'))
 for f in files:assert f.read_bytes()==(q/'result'/f.relative_to(p/'result')).read_bytes(),('first output difference',f)
 la,lb=(p/'run.log').read_text(),(q/'run.log').read_text()
 traces=lambda s:[l for l in s.splitlines() if l.startswith('[SF PLAN]') or l.startswith('[SF] Step time:') or l.startswith('[SF] State closure:') or l.startswith('[SF] Eulerian step:')]
 assert traces(la)==traces(lb),('first execution/clock difference',a['name'])
 if (p/'numerical.txt').exists():assert (p/'numerical.txt').read_bytes()==(q/'numerical.txt').read_bytes(),('Eulerian matrix/solution difference',a['name'])
 groups=defaultdict(list)
 for step,op in re.findall(r'\[SF PLAN\] step=(\d+).*?op=(\S+)',lb):groups[int(step)].append(op)
 assert groups, (a['name'],'missing Plan trace')
 for step,ops in groups.items():
  counts=Counter(ops);clock='ee.time.commit' if a['name'].startswith('eulerian') else 'time.commit'
  assert counts[clock]==1
  if a['name'] in ('sst','kEpsilon','ghostSST','les'):
   closure='turbulence.closure.refresh' if a['name']=='les' else 'turbulence.advance'
   # Operation spelling comes from the frozen production Plan, not a new schedule.
   found=[i for i,x in enumerate(ops) if x.startswith('turbulence.')]
   assert len(found)==1 and ops.index('flow.dt.compute')<found[0]<ops.index('flow.step.begin')
 frames=[]
 for f in files:
  if f.suffix!='.vts':continue
  x=ET.parse(f);piece=x.find('.//Piece');xyz=list(map(float,piece.find('./Points/DataArray').text.split()))
  data={v.attrib['Name']:list(map(float,v.text.split())) for v in piece.findall('./PointData/DataArray')}
  assert all(math.isfinite(v) for values in data.values() for v in values)
  stat={k:dict(min=min(v),max=max(v),nodalL1=math.fsum(abs(z) for z in v),nodalL2=math.sqrt(math.fsum(z*z for z in v))) for k,v in data.items() if not k.startswith('BC_')}
  frames.append(dict(file=f.name,stats=stat))
 settingsJson=json.loads(subprocess.check_output([str(root/'build/native_caseSettings'),str(q)],text=True,cwd=base))
 times=[d.attrib['timestep'] for d in ET.parse(next((q/'result').glob('*.pvd'))).findall('.//DataSet')]
 item=dict(name=a['name'],source=a['source'],before=a,after=b,steps=len(groups),frames=len(frames),times=times,outputsExact=True,planClockClosureExact=True,
  eulerianMatrixExact=(q/'numerical.txt').exists(),extent=piece.attrib['Extent'],bounds=[[min(xyz[d::3]),max(xyz[d::3])] for d in range(3)],final=frames[-1],settings=settingsJson)
 settings.append(item);summary.append({k:v for k,v in item.items() if k not in ('settings','before','after','final')})
 (out/(a['name']+'-frame-statistics.json')).write_text(json.dumps(frames,indent=2))
(out/'comparison.json').write_text(json.dumps(summary,indent=2));(out/'case-settings.json').write_text(json.dumps(settings,indent=2))
lines=['# Once-step RAS 逐例设置与运行证据','','2026-10-16。所有8对输入、全部输出文件逐字节比较，Plan/clock/closure trace全等，未使用容差。字段初值/边界及模型块由正式 CaseIO 工具读取；磁盘未注册文件不参与执行。以下保存注册参数，不将文件中的 inactive 字段或模型当作已测试功能。','','共同命令：`baseline/sonicSolver run <beforeFinal/name>` 与 `build/sonicSolver run <afterFinal/name>`，脚本以各case为cwd，SF_PLAN_TRACE=1。每例完整命令、binary SHA、实际终点、网格范围、所有字段IC/BC、数值参数及注册模型见下。坐标按当前CFD输入SI解释（m/s、Pa、kg/m³、Pa·s）；mesh本身不携带长度单位标签。均serial单source zone/单patch，parallel.enabled=false时nProcs不生效。','','短时回归不等于长期湍流壁面物理验收；Eulerian采用实际PIMPLE而非配置中保留的density RK；非IBM案例的wallClosure、热壁面、ILW不适用。Ghost采用stationaryNoSlipAdiabatic、Euler ILW关闭、Central2黏性及原壁距/SST壁值；未测试forcing/DLM/KKT。统计为节点extrema/L1/L2，包含输出点，Ghost含标记点；不是物理体积分或守恒预算。','']
for item in settings:
 lines+=['## '+item['name'],'','```json',json.dumps(item,ensure_ascii=False,indent=2),'```','']
(out/'case-settings.md').write_text('\n'.join(lines))
print(json.dumps(summary,indent=2))
