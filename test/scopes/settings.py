"""Retain actual per-case inputs and outputs for the report, never infer settings from names."""
from pathlib import Path
import json,re,xml.etree.ElementTree as E
base=Path(__file__).resolve().parent
families=[('数值对照',base/'before'),('时间作用域',base/'native/cases'),('原压力串并行',base/'pressure/before'),('Sod四rank',base/'mpi/before')]
lines=['# 本轮逐例实际输入与完成状态','', '本附件由 test/scopes/settings.py 从本轮实际配置、网格和输出读取。before/after 的具体差异另见 scope-migration.json；无单位转换，几何数值按 case 的 SI 输入约定解释。source/typed 输入并不意味着未活动模块参与执行。','']
def compact(file):
 text=file.read_text()
 try:v=json.loads(text)
 except ValueError:return text.strip()
 v.pop('SonicFile',None)
 if file.name=='internalField.yaml':
  for key,values in v.items():
   if not isinstance(values,dict):continue
   points={k:x for k,x in values.items() if k.startswith('point') and isinstance(x,(int,float))}
   if points:
    rest={k:x for k,x in values.items() if k not in points};rest['explicitPointValues']={'count':len(points),'min':min(points.values()),'max':max(points.values()),'source':'actual fields/internalField.yaml (retained)'};v[key]=rest
 return json.dumps(v,ensure_ascii=False,indent=2)
for group,folder in families:
 for case in sorted(folder.iterdir()) if folder.exists() else []:
  if not (case/'case.yaml').is_file():continue
  if group=='IBM 完整 no-slip Ghost' and case.name!='cylinderFlowViscousGhost':continue
  lines+=['## '+group+' / '+case.name,'', '路径：`'+str(case)+'`；日志 `run.log`；数值目录 `result/`。','']
  mesh=case/'mesh/mesh.sfm'
  if mesh.exists():
   content=mesh.read_text().splitlines();at=content.index('#Information');shape=content[at+1];at=content.index('#Point')+1;points=[]
   while at<len(content) and content[at].strip() and not content[at].startswith('#'):
    points.append(list(map(float,content[at].split())));at+=1
   lines+=['网格文件头（nodes + NG）：`'+shape+'`；实际坐标范围：`'+json.dumps([[min(p[i] for p in points),max(p[i] for p in points)] for i in range(3)])+'`。','']
  pvds=list((case/'result').glob('*.pvd'))
  if pvds:
   frames=E.parse(pvds[0]).findall('.//DataSet');lines+=['实际输出：'+str(len(frames))+' 帧；末帧 physical time=`'+frames[-1].get('timestep')+'`。','']
  log=case/'run.log'
  if log.exists():
   text=log.read_text(errors='replace');clocks=re.findall(r'Step time: time=([^,\n]+), dt=([^\n]+)',text)
   if clocks:lines+=['实际 flow clock：'+str(len(clocks))+' 次提交；最后 time/dt=`'+str(clocks[-1])+'`。','']
   if 'Fatal' in text:lines+=['未完成/拒绝：`'+next(line for line in reversed(text.splitlines()) if 'Fatal' in line)+'`。','']
  for rel in ['case.yaml','solvers/runtime.yaml','solvers/numerics.yaml','solvers/equations.yaml','solvers/algorithm.yaml','solvers/state.yaml','solvers/equationRegistry.yaml','solvers/how.yaml','solvers/which.yaml','fields/fields.yaml','fields/internalField.yaml','fields/boundaries.yaml','models/models.yaml']:
   file=case/rel
   if file.is_file():lines+=['### '+rel,'','```yaml',compact(file),'```','']
  for file in sorted((case/'models').glob('*.yaml')):
   if file.name!='models.yaml' and not file.name.startswith('._'):lines+=['### models/'+file.name,'','```yaml',compact(file),'```','']
(base/'evidence/settings.md').write_text('\n'.join(lines)+'\n')
print('Per-case input appendix generated from current baseline/native copies')
