from pathlib import Path
import subprocess,json,os,hashlib,math,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent;root=base.parents[1];e=base/'evidence'
latest=json.loads((e/'latest.json').read_text());m=json.loads(Path(latest['manifest']).read_text());protection=json.loads((e/'protection-delivery.json').read_text())
assert m['sha256']==hashlib.sha256((root/'build/sonicSolver').read_bytes()).hexdigest()
env={**os.environ,'TMPDIR':str(e),'DEVELOPER_DIR':'/Library/Developer/CommandLineTools'}
items=[]
for record in m['positive']:
 items.append({'category':'native/reference','record':record,'case':record['case']})
for record in protection:items.append({'category':'frozen-production','record':record,'case':record['runs'][1]['case']})
for item in items:
 p=Path(item['case']);utility=root/'build/native_caseSettings'
 item['registeredSettings']=json.loads(subprocess.check_output([str(utility),str(p)],cwd=base,env=env,text=True))
 frames=sorted(f for f in (p/'result').glob('*.vts') if not f.name.startswith('._'))
 doc=ET.parse(frames[-1]);piece=doc.find('.//Piece');xyz=list(map(float,piece.find('./Points/DataArray').text.split()))
 axes=[sorted(set(xyz[d::3])) for d in range(3)]
 weights=[{v:((a[i+1]-v)/2 if i==0 else (v-a[i-1])/2 if i==len(a)-1 else (a[i+1]-a[i-1])/2) for i,v in enumerate(a)} for a in axes]
 w=[math.prod(weights[d][xyz[3*i+d]] for d in range(3)) for i in range(len(xyz)//3)]
 item['bounds']=[[a[0],a[-1]] for a in axes];item['nodes']=[len(a) for a in axes]
 item['finalMetrics']={}
 for a in piece.findall('./PointData/DataArray'):
  v=list(map(float,a.text.split()));assert all(math.isfinite(x) for x in v)
  item['finalMetrics'][a.attrib['Name']]={'min':min(v),'max':max(v),'nodalL1':math.fsum(abs(x) for x in v),'nodalL2':math.sqrt(math.fsum(x*x for x in v))}
  if len(v)==len(w):item['finalMetrics'][a.attrib['Name']]['rectangularTrapezoidalIntegral']=math.fsum(a*b for a,b in zip(w,v))
(e/'case-settings.json').write_text(json.dumps(items,ensure_ascii=False,indent=2))
lines=['# Native Equation Closure 逐例实际配置与证据','', '2026-10-19。生产CLI与C++参考分别列明；不是将新增native功能冒充修改前可运行基线。所有输入路径、命令、binary SHA、frame/time及registeredSettings见逐例JSON。fields输出未激活的配置不等同于执行。', '', '全新标量例无EOS/NS/IBM/ILW/mu/Pr。Frozen CFD按原已注册模型和solver配置执行；Eulerian实际是PIMPLE，不能由配置RK4字样宣称RK4。全部串行单patch，无多rank MPI回归。坐标与D按一致SI解释，SFM没有单位自动转换。', '', 'finalMetrics为输出节点统计；积分是矩形Cartesian域梯形积分，Ghost包含solid输出点，不作为cut-cell流体体积/守恒预算。Eulerian相焓不是总能量。未展示自动缩放云图。', '', 'Final CLI SHA256: `'+m['sha256']+'`','']
for item in items:lines+=['## '+item['category']+' — '+Path(item['case']).name,'','```json',json.dumps(item,ensure_ascii=False,indent=2),'```','']
(e/'case-settings.md').write_text('\n'.join(lines))
print(len(items),'actual positive settings collected; native frames',sum(x['frames'] for x in m['positive']),'protection frames',sum(x['frames'] for x in protection))
