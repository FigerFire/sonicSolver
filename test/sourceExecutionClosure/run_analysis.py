from pathlib import Path
from collections import Counter, defaultdict
import hashlib
import json
import math
import re
import subprocess
import xml.etree.ElementTree as ET

base=Path(__file__).resolve().parent
root=base.parents[1]
evidence=base/"evidence"
records=json.loads((evidence/"cases.json").read_text())
reports=[]
catalog={}
settings=[]
for record in records:
    case=Path(record["case"])
    parsed=subprocess.run([str(root/"build/native_caseSettings"),str(case)],capture_output=True,text=True,check=True)
    config=json.loads(parsed.stdout)
    key=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()[:12]
    catalog[key]=config
    item=dict(case=str(case.relative_to(root)),settings=key,command=record["command"],
        binary=record["binary"],exit=record["exit"],actual=record["physical_lines"][-1:] or record["error"],
        log=str((case/"run.log").relative_to(root)),outputs=str((case/"result").relative_to(root)))
    settings.append(item)
    if record["exit"]:
        continue
    text=(case/"run.log").read_text()
    grouped=defaultdict(list)
    for step,operation in re.findall(r"\[SF PLAN\] step=(\d+).*?op=([^\s]+)",text):
        grouped[int(step)].append(operation)
    coupling=config["algorithm"]["PIMPLE"]
    outer=int(coupling["outerCorrectors"])
    corrections=int(coupling["pressureCorrectors"])*(int(coupling["nonOrthogonalCorrectors"])+1)
    for step,operations in grouped.items():
        counts=Counter(operations)
        for op in ("ee.dt.compute","ee.step.begin","ee.step.commit","ee.time.commit"):
            assert counts[op]==1,(record["name"],step,op,counts[op])
        for op in ("ee.interphase.compute","ee.sources.assemble","ee.sources.validate","ee.energy.solve"):
            assert counts[op]==outer,(record["name"],step,op,counts[op])
        assert counts["ee.pressure.solve"]==outer*corrections
        for index,op in enumerate(operations):
            if op=="ee.sources.assemble":assert operations[index-1]=="ee.interphase.compute"
        assert operations[-2:]==["ee.step.commit","ee.time.commit"]
    frames=[]
    for file in sorted((case/"result").glob("*.vts")):
        if file.name.startswith("."):continue
        piece=ET.parse(file).find(".//Piece")
        xyz=list(map(float,piece.find("./Points/DataArray").text.split()))
        data={a.attrib["Name"]:list(map(float,a.text.split())) for a in piece.findall("./PointData/DataArray")}
        assert all(math.isfinite(v) for values in data.values() for v in values)
        stats={name:dict(min=min(values),max=max(values),nodalL1=math.fsum(abs(v) for v in values),
            nodalL2=math.sqrt(math.fsum(v*v for v in values))) for name,values in data.items() if not name.startswith("BC_")}
        errors={};scaledErrors={};integrals={}
        if "SourceOnly" in record["name"]:
            n=len(data["p"])
            for kind,names in {
                "mass":["massSource.water","massSource.air"],
                "momentumX":["momentumSourceX.water","momentumSourceX.air"],
                "momentumY":["momentumSourceY.water","momentumSourceY.air"],
                "momentumZ":["momentumSourceZ.water","momentumSourceZ.air"],
                "energy":["energySource.water","energySource.air"],
            }.items():
                error=0;scaled=0
                for cell in range(n):
                    expected=0
                    if kind=="energy":
                        expected=sum(data[key][cell] for key in ("interphaseMechanicalHeating","wallBoilingQconv","wallBoilingQquench","wallBoilingQevap"))
                    actual=sum(data[key][cell] for key in names)
                    scale=1+abs(expected)+sum(abs(data[key][cell]) for key in names)
                    error=max(error,abs(actual-expected))
                    scaled=max(scaled,abs(actual-expected)/scale)
                    assert abs(actual-expected)<=2e-12*scale,(file,kind,cell,actual,expected)
                errors[kind]=error
                scaledErrors[kind]=scaled
            axes=[sorted(set(xyz[d::3])) for d in range(3)]
            assert math.prod(len(a) for a in axes)==n
            weights=[{x:((a[i+1]-x)/2 if i==0 else (x-a[i-1])/2 if i==len(a)-1
                else (a[i+1]-a[i-1])/2) for i,x in enumerate(a)} for a in axes]
            pointWeights=[math.prod(weights[d][xyz[3*i+d]] for d in range(3)) for i in range(n)]
            for label,keys in {
                "massRate":["massSource.water","massSource.air"],
                "forceX":["momentumSourceX.water","momentumSourceX.air"],
                "forceY":["momentumSourceY.water","momentumSourceY.air"],
                "forceZ":["momentumSourceZ.water","momentumSourceZ.air"],
                "enthalpySourcePower":["energySource.water","energySource.air"],
                "mechanicalHeatingPower":["interphaseMechanicalHeating"],
                "wallInputPower":["wallBoilingQconv","wallBoilingQquench","wallBoilingQevap"],
            }.items():
                integrals[label]=math.fsum(pointWeights[i]*sum(data[key][i] for key in keys) for i in range(n))
            expected=0 if file.name.endswith("_000000.vts") or "rpi" not in record["name"] else 70
            assert abs(integrals["wallInputPower"]-expected)<=1e-10*max(1,expected),(file,integrals)
        frames.append(dict(file=str(file.relative_to(root)),stats=stats,sourceBudgetMaxAbsoluteErrors=errors,
            sourceBudgetMaxScaledErrors=scaledErrors,cartesianTrapezoidalSourceIntegrals=integrals))
        item["mesh"]=dict(extent=piece.attrib["Extent"],nodes=len(xyz)//3,
            bounds=[[min(xyz[d::3]),max(xyz[d::3])] for d in range(3)])
    item["actualCommittedSteps"]=len(grouped)
    reports.append(dict(name=record["name"],version=record["version"],steps=len(grouped),outer=outer,
        pressurePassesPerOuter=corrections,traceValid=True,frames=frames))
(evidence/"trace-and-budgets.json").write_text(json.dumps(reports,indent=2))
(evidence/"settings.json").write_text(json.dumps(settings,ensure_ascii=False,indent=2))
(evidence/"config-catalog.json").write_text(json.dumps(catalog,ensure_ascii=False,indent=2))
lines=["# Source execution closure：逐例设置与证据", "", "日期：2026-10-15。配置由生产 CaseIO 读取；before/after 使用同一输入；失败日志也保留。", "",
    "所有案例来源 test/nativeArchitectureClosure/eulerianNative/laminar-all。IC/BC、EOS/物性和 runtime 的完整值见各配置表。坐标按 SI 使用；无 IBM/STL，wallClosure 与 ILW 不适用。每个 case 的 result/ 独立，不覆盖正式算例。", "",
    "共同物理系统：逐相 phaseMass/momentum/enthalpy + shared pressure，N−1 相质量闭合；PIMPLE、Upwind、已有隐式动量/焓矩阵与 HYPRE。numerics 中保留 classicalRK4，但此 Eulerian provider 的实际 HOW 是 PIMPLE 物理步，不执行 density RK stages。未验证 RK4 相变。", "",
    "共同 IC：p=101325 Pa，T.water=T.air=300 K，U.water=U.air=0，alpha.water=0.9、alpha.air=0.1。水 perfectLiquid：rho=998.2、mu=0.001002、Cp=Cv=4182、k=0.598；空气 perfectGas：配置rho=1.1766、R=287.05、mu=1.846e−5、Cp=1007、Cv=718、k=0.0262，实际密度由当前相闭合恢复。", "",
    "共同 BC：Left 固定 U=0、T=300、alpha=各相初值、p=101325；Right U/T/alpha/p zeroGradient；bottom/top U=0，T/alpha/p zeroGradient；frontBack empty。rho 各相边界值及 fields registry 完整声明见 JSON。RPI 的 wallHeatFlux 是独立 source setting，不改 T 的 zeroGradient BC。", "",
    "共同 runtime：startTime=0、endTime=0.0002 s、无显式endStep、CFL=0.5、maxDeltaT=0.0001 s、初始输出及每步输出。serial，单 source zone；实际dt和终点逐例取run.log。", "",
    "实验差异：disabled* 仅显式注册停用模型，缺少参数不得触发active校验；embeddedDisabledRpi 使用multiPhase内嵌块。Lee采用Ce=Cc=0.001、L=2e6，Tsat=290为蒸发、310为冷凝。RPI fixed 对照为qConv=1000、qQuench=2000、固定气泡闭式和壁温390，随流场变化触发原有预算拒绝；Balanced改用已有hConv=50的coefficient闭式与heatFluxBalance壁温求根，qWall=7000、qQuench=2000，未改容差。Outer3仅outer改为3。SourceOnly额外移除MRF注册和对应source selection；磁盘未注册的MRF文件不执行。其他实验保留原MRF：all zone、中心(0.5,0.1,0)、z轴、omega=1。", "",
    "source-only验收检查瞬时源率预算，并非包含全部边界通量、压力功的整体流体总能量守恒证明。普通sum/L1/L2只是节点统计；不把它当作物理积分。", ""]
for item in settings:
    lines += ["## "+item["case"],"","```json",json.dumps(item,ensure_ascii=False,indent=2),"```",""]
lines += ["# 完整已注册配置","", "各case的settings键索引下表；未注册文件不参与读取。", ""]
for key,value in catalog.items():
    lines += ["## "+key,"","```json",json.dumps(value,ensure_ascii=False,indent=2),"```",""]
(evidence/"case-settings.md").write_text("\n".join(lines))
print(len(reports),"successful launches; trace counts and",sum(len(r["frames"]) for r in reports),"frames checked")
