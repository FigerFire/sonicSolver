from pathlib import Path
import subprocess,json,re,hashlib,struct
B=Path(__file__).resolve().parent;R=B.parents[1];rows=[]
def add(p,source,purpose,command=None):
 if p.is_dir() and (p/'case.yaml').exists():rows.append(dict(path=p,source=source,purpose=purpose,command=command))
for group in ['d2After','d2MultiAfter']:
 for p in sorted((B/group).glob('*/*')):add(p,'test/levelSetMpiD1 cases.py / multiPatchCloseout','D2 纯 ownership 对照；包括 exact C2/C6/C12，MPI composite 与 serial Unsupported。')
for rec in json.loads((B/'evidence/protection-manifest.json').read_text()):
 if rec['stage']=='after':add(B/'protection'/(rec['case']+'-after'),rec['source'],'同数学／同输入短时冻结前后对照；unsupported KKT 是能力负例。',rec['command'])
for rec in json.loads((B/'evidence/mpi-protection-manifest.json').read_text()):
 if rec['version']=='after':add(B/'mpiProtection'/(rec['name']+'-after'),rec['source'],'同 MPI 分区配置前后数值字段对照；不扩展物理能力。',rec['cmd'])
for rec in json.loads((B/'evidence/eulerianNative-manifest.json').read_text()):add(B/'eulerianNative'/rec['case'],'test/levelSetMpiD1/eulerianFinal/'+rec['case'],'移除副本旧 generator 注册，复用冻结 SFM；比较矩阵、dt、诊断和 VTS。',rec['command'])
for name in ['guiMeshFinal','thermalOneVerified','thermalFiveVerified','thermalBoundaryRegisteredOneVerified','thermalBoundaryRegisteredFiveVerified','thermalFamilyNegative','thermalOneFinal','cliDefault','cliVerified','thermalOne','thermalFive','thermalRkNegative','thermalBoundaryOne','thermalBoundaryFive','thermalBoundaryRegisteredOne','thermalBoundaryRegisteredFive','guiSerial','guiMpi','guiSerialFinal','guiMpiFinal']:
 add(B/name,'Application::createNativeCase / 本轮生成的明确副本','native IO／独立温度／GUI 测试；原 thermalCase 未注册实际 thermal 模块，不能当作旧 CFD baseline。')
for p in sorted((R/'test/guiNativePreview').glob('20261014*')):
 if (p/'case.yaml').exists():add(p,'GUI native template test case','独立 GUI preview：sonicSolver steps 1；不是完整时域。')
manifest=[];catalog={};sections=['# 全模块收口：实际运行算例与设置附录','', '日期：2026-10-14。仅已注册文档由生产 CaseIO 解析。以下逐例标明来源、命令、网格、终点和配置索引；配置表包含全部初值、边界、EOS/物性、模型参数、WHAT/STATE/HOW/WHICH 与并行设置。`enabled:false` 是已配置但未执行，不以文件存在推断功能。','', 'SFM 本身无独立单位标签，坐标数值不自动换单位；热模板按 SI（T 为 K、速度 m/s、D m²/s、时间 s）解释。Sod 是标准归一化激波管。IBM 原始算例保留其坐标／物性数值，本轮不作真实圆柱物理标度、长时尾流或网格收敛验收。网格前三个 Information 值为物理节点维数，结构单元数各减一；halo depth 为第四值。','', 'before/after 对照使用同一份配置，保护组统一把停止时间改为 0.001 s（pressure 0.0002 s），按步每步输出；其他数学、IC/BC/模型/数值参数不变。D2 endStep=3、endTime=0.001；重初始化仅 serial，steps=3、pseudoDt=0.01，不增加 physical commit。所有输出为各 case/result/，checkpoint 为 case/checkpoint/。旧 PVD step 标签与新 physical-time 标签的差异单独说明，不能算作字段误差。','', '热温度常量版本仅是接线 smoke；thermalBoundaryRegistered* 另加显式注册的 custom thermalTemperature，左310 K fixedValue、右300 K fixedValue、上下zeroGradient、前后empty。EOS T 从 rhoE 派生，独立温度不写 rhoE。kernel harness 才是原 Forward Euler 算术的独立等价证据。','', 'GUI offscreen 仅证明 Qt/Controller/QProcess，OpenGL 不可用；host Cocoa 图像另行检查。开发中失败的 cliDefault 网格、thermalBoundary 未注册字段及历史 generator 输入均保留；修正用独立目录，未将它们算作成功数值回归。','']
for index,r in enumerate(rows,1):
 p=r['path'];relative=str(p.relative_to(R));log=p/'run.log'
 if not log.exists():log=p/'explain.log'
 text=log.read_text(errors='replace') if log.exists() else ''
 cmd=r['command'] or [str(R/'build/sonicSolver'),'run',str(p)]
 result=subprocess.run([str(R/'build/native_caseSettings'),str(p)],capture_output=True,text=True)
 config=json.loads(result.stdout) if result.returncode==0 else {'readerError':result.stderr.strip()}
 if not r['command'] and config.get('runtime',{}).get('parallel',{}).get('enabled',False):
  cmd=['mpirun','-np',str(config['runtime']['parallel']['nProcs'])]+cmd
 # Case names/output job labels do not describe different physics.
 normalized=json.loads(json.dumps(config));normalized.pop('name',None)
 normalized.get('runtime',{}).get('output',{}).pop('jobName',None)
 key=hashlib.sha256(json.dumps(normalized,sort_keys=True).encode()).hexdigest()[:12];catalog.setdefault(key,normalized)
 mesh=[]
 files=config.get('mesh',{}).get('files',[])
 for name in files:
  f=p/name
  if not f.exists():mesh.append(dict(file=name,missing=True));continue
  s=f.read_text(errors='replace');info=re.search(r'#Information\s*\n([^\n]+)',s)
  point=re.search(r'#Point\s*\n(.*?)\n#',s,re.S)
  if point:
   xyz=[tuple(map(float,l.split())) for l in point.group(1).splitlines() if len(l.split())==3]
   if xyz:mesh.append(dict(file=name,information=info.group(1) if info else '?',points=len(xyz),bounds=[(min(c[a] for c in xyz),max(c[a] for c in xyz)) for a in range(3)],sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 stl=[]
 for f in (p/'mesh').rglob('*.stl'):
  if f.name.startswith('._'):continue
  data=f.read_bytes();vertices=[]
  if len(data)>=84 and 84+50*struct.unpack_from('<I',data,80)[0]==len(data):
   for n in range(struct.unpack_from('<I',data,80)[0]):v=struct.unpack_from('<12f',data,84+n*50);vertices.extend((v[3:6],v[6:9],v[9:12]))
  else:vertices=[tuple(map(float,v)) for v in re.findall(rb'vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)',data)]
  if vertices:stl.append(dict(file=str(f.relative_to(p)),bounds=[(min(v[a] for v in vertices),max(v[a] for v in vertices)) for a in range(3)]))
 times=re.findall(r'\[SF\] Step time:.*',text);steps=set(re.findall(r'\[SF PLAN\] step=(\d+)',text));ee=re.findall(r'\[SF\].*(?:Time=|time=).*',text)
 final=times[-1] if times else (ee[-1] if ee else '没有成功 advance 的日志；GUI 日志见 evidence 对应文件。')
 item=dict(case=relative,source=r['source'],purpose=r['purpose'],settings=key,command=cmd,mesh=mesh,stl=stl,actualStepCount=len(times) or len(steps),actualFinal=final,log=str(log.relative_to(R)))
 manifest.append(item)
 sections+=['## '+str(index)+'. '+relative,'','来源：`'+r['source']+'`。目的：'+r['purpose'],'','命令：`'+' '.join(cmd)+'`。GUI 组实际 QProcess 命令／steps limit 以 evidence/gui-*.log 和 GUI controller 测试配置为准。','', '实际终点：`'+final+'`；成功 Step time 记录 '+str(len(times))+' 条，trace step 集合 '+str(len(steps))+' 项。配置的 endTime/endStep 另见配置 **'+key+'**，不能用帧数冒充步数。','', '网格：`'+json.dumps(mesh,ensure_ascii=False)+'`。','', '物体位置／尺寸：`'+json.dumps(stl,ensure_ascii=False)+'`；无 STL 时不适用，IBM target、center 和壁面参数以配置及实际启用模块为准。','', '完整 WHAT/STATE、IC/BC、EOS/transport、HOW/WHICH、MPI/rank、CFL、maxDeltaT、输出间隔：配置 **'+key+'**。日志 `'+item['log']+'`；output `'+relative+'/result/`，checkpoint 若启用则 `'+relative+'/checkpoint/`。','']
sections+=['# 注册配置表','', '按全部数学／状态／执行／数值配置的 canonical JSON 去重。注释和 YAML 风格未作为数学参数；原文件与 reader 数据快照均保留。文件中未执行的 enabled:false 模块也列出；本表不将其计为测试能力。','']
for key,c in catalog.items():sections+=['## 配置 '+key,'','```json',json.dumps(c,ensure_ascii=False,indent=2),'```','']
(B/'evidence/case-settings.md').write_text('\n'.join(sections));(B/'evidence/case-settings.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));(B/'evidence/config-catalog.json').write_text(json.dumps(catalog,ensure_ascii=False,indent=2))
print('documented',len(manifest),'cases with',len(catalog),'unique registered configurations')
