# -*- coding: utf-8 -*-
"""Compare frozen pre-migration output without tolerances; retain real case settings."""
from pathlib import Path
from collections import Counter,defaultdict
import hashlib,json,math,re,subprocess,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent;root=base.parents[1];out=base/'evidence'
before=json.loads((out/'cases-before.json').read_text());after=json.loads((out/'cases-verified.json').read_text())
assert len(before)==len(after)==10 and len({x['name'] for x in after})==10, 'Incomplete CFD manifest'
assert all(x['binary']==hashlib.sha256((root/'build/sonicSolver').read_bytes()).hexdigest() for x in after), 'Final executable changed since CFD run'
for recipe in ('forwardEuler','classicalRK4'):
 assert (out/'beforeKernel'/('checkpoint-'+recipe+'.txt')).read_bytes()==(out/'verifiedKernel'/('checkpoint-'+recipe+'.txt')).read_bytes(), ('first kernel coefficient/RHS/flux/stage difference',recipe)
 assert (out/'beforeKernel'/('checkpoint-'+recipe+'.txt')).read_bytes()==(out/'unit/stageCoefficientEvidence'/('checkpoint-'+recipe+'.txt')).read_bytes(), ('final CTest kernel difference',recipe)

summary=[];settings=[]
for a,b in zip(before,after):
 assert a['name']==b['name'] and a['exit']==b['exit']==0
 p,q=Path(a['case']),Path(b['case'])
 for f in p.rglob('*'):
  if f.is_file() and not f.name.startswith('._') and f.suffix in ('.yaml','.sfm','.STL') and 'result' not in f.parts:
   assert f.read_bytes()==(q/f.relative_to(p)).read_bytes(),('input changed',f)
 files=sorted(f for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('.'))
 assert [str(f.relative_to(p/'result')) for f in files]==sorted(str(f.relative_to(q/'result')) for f in (q/'result').rglob('*') if f.is_file() and not f.name.startswith('.'))
 for f in files:assert f.read_bytes()==(q/'result'/f.relative_to(p/'result')).read_bytes(),('first output difference',f)
 if (p/'checkpoint').exists():
  def checkpoints(case):
   return {re.sub(r'^flow-\d+-','',f.name):f.read_bytes() for f in (case/'checkpoint').glob('*.csv') if not f.name.startswith('._')}
  assert checkpoints(p) and checkpoints(p)==checkpoints(q),('first flow checkpoint difference',a['name'])
 la,lb=(p/'run.log').read_text(),(q/'run.log').read_text()
 traces=lambda s:[re.sub(r'0x[0-9a-fA-F]+','<address>',l) for l in s.splitlines() if l.startswith('[SF PLAN]') or l.startswith('[SF TRACE]') or l.startswith('[SF] Step time:') or l.startswith('[SF] State closure:') or l.startswith('[SF] Eulerian step:')]
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
  axes=[sorted(set(xyz[d::3])) for d in range(3)]
  weights=[{v:((a[i+1]-v)/2 if i==0 else (v-a[i-1])/2 if i==len(a)-1 else (a[i+1]-a[i-1])/2) for i,v in enumerate(a)} for a in axes]
  quadrature=[math.prod(weights[d][xyz[3*i+d]] for d in range(3)) for i in range(len(xyz)//3)]
  integrals={k:math.fsum(w*z for w,z in zip(quadrature,data[k])) for k in ('rho','MomentumX','MomentumY','MomentumZ','TotalEnergyDensity') if k in data}
  frames.append(dict(file=f.name,stats=stat,rectangularTrapezoidalIntegrals=integrals))
 settingsJson=json.loads(subprocess.check_output([str(root/'build/native_caseSettings'),str(q)],text=True,cwd=base))
 times=[d.attrib['timestep'] for d in ET.parse(next((q/'result').glob('*.pvd'))).findall('.//DataSet')]
 item=dict(name=a['name'],source=a['source'],before=a,after=b,steps=len(groups),frames=len(frames),times=times,outputsExact=True,planClockClosureExact=True,
  eulerianMatrixExact=(q/'numerical.txt').exists(),extent=piece.attrib['Extent'],bounds=[[min(xyz[d::3]),max(xyz[d::3])] for d in range(3)],final=frames[-1],settings=settingsJson)
 settings.append(item);summary.append({k:v for k,v in item.items() if k not in ('settings','before','after','final')})
 (out/(a['name']+'-frame-statistics.json')).write_text(json.dumps(frames,indent=2))
(out/'comparison.json').write_text(json.dumps(summary,indent=2));(out/'case-settings.json').write_text(json.dumps(settings,indent=2))
lines=['# Stage coefficient 逐例设置与运行证据','','2026-10-19。所有10对输入、全部输出文件逐字节比较，Plan/clock/closure trace全等，未使用容差。字段初值/边界及模型块由正式 CaseIO 工具读取；磁盘未注册文件不参与执行。以下保存注册参数，不将文件中的 inactive 字段或模型当作已测试功能。','','共同命令：`baseline/sonicSolver run <before/name>` 与 `build/sonicSolver run <verified/name>`，脚本以各case为cwd，SF_PLAN_TRACE=1。每例完整命令、binary SHA、实际终点、网格范围、所有字段IC/BC、数值参数及注册模型见下。坐标按当前CFD输入SI解释（m/s、Pa、kg/m³、Pa·s）；mesh本身不携带长度单位标签。均serial单source zone/单patch，parallel.enabled=false时nProcs不生效。','','短时回归不等于长期湍流壁面物理验收；Eulerian采用实际PIMPLE而非配置中保留的density RK；非IBM案例的wallClosure、热壁面、ILW不适用。Ghost采用stationaryNoSlipAdiabatic、Euler ILW关闭、Central2黏性及原壁距/SST壁值；未测试forcing/DLM/KKT。统计包含节点extrema/L1/L2，以及可用rho/动量/总能量密度的Cartesian矩形域梯形积分；Ghost积分包含solid输出点，不是cut-cell流体体积分。Eulerian缺少总能量输出时不以相焓代替。均不宣称完成守恒预算物理验收。','']
for item in settings:
 lines+=['## '+item['name'],'','```json',json.dumps(item,ensure_ascii=False,indent=2),'```','']
(out/'case-settings.md').write_text('\n'.join(lines))
print(json.dumps(summary,indent=2))
