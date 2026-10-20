"""Short native cases: preserve inputs, binaries, commands and every launch result."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

base = Path(__file__).resolve().parent
root = base.parents[1]
source = root / "test/nativeArchitectureClosure/eulerianNative/laminar-all"
evidence = base / "evidence"
manifest = json.loads((evidence / "cases.json").read_text()) if len(sys.argv)>1 else []

def prepare(name, phase=None, embedded=False):
    case = base / name
    if case.exists():
        raise RuntimeError(f"Refusing to replace existing evidence: {case}")
    shutil.copytree(source, case, ignore=shutil.ignore_patterns("result", "checkpoint", "*.log", "._*", "numerical.txt"))
    if phase:
        if embedded:
            f = case / "models/multiPhase.yaml"
            f.write_text(f.read_text()+"\nphaseChange:\n"+"\n".join("  "+s for s in phase.splitlines())+"\n")
        else:
            f = case / "models/models.yaml"
            f.write_text(f.read_text()+"\nphaseChange:\n  type: phaseChange\n  file: models/phaseChange.yaml\n")
            (case / "models/phaseChange.yaml").write_text("SonicFile:\n  object: models\n  type: phaseChange\n"+phase+"\n")
    if "Outer3" in name:
        f=case / "solvers/algorithm.yaml"
        f.write_text(f.read_text().replace("outerCorrectors: 1", "outerCorrectors: 3"))
    if "SourceOnly" in name:
        f=case / "models/models.yaml"
        f.write_text(f.read_text().replace("MRF:\n  type: MRF\n  file: models/MRF.yaml\n", ""))
        f=case / "solvers/numerics.yaml"
        f.write_text(f.read_text().replace("sources:\n  default: MRF\n", ""))
    return case

disabled = "enabled: false\nphaseChangeModel: RPI\nphases: [water, air]"
lee = "enabled: true\nphaseChangeModel: Lee\nphases: [water, air]\nevaporationCoefficient: 0.001\ncondensationCoefficient: 0.001\nsaturation:\n  Tsat: {tsat}\n  latentHeat: 2000000"
rho = 101325/(287.05*300)
nucleation = .002 / (math.pi*.002**3/6*rho*20)
rpi = f"""enabled: true
phaseChangeModel: RPI
phases: [water, air]
saturation:
  Tsat: 290
  latentHeat: 2000000
RPI:
  departureDiameterModel: explicit
  departureDiameter: 0.002
  departureFrequencyModel: explicit
  departureFrequency: 20
  nucleationSiteDensityModel: explicit
  nucleationSiteDensity: {nucleation:.17g}
  convectiveHeatFluxModel: fixed
  convectiveHeatFlux: 1000
  quenchingHeatFluxModel: fixed
  quenchingHeatFlux: 2000
  wallTemperatureModel: explicit
  maximumWallVaporFraction: 0.9
  heatFluxBalanceTolerance: 0.0001
wallBoiling:
  bottomHeating:
    patch: bottom
    heatFlux: 7000
    wallTemperature: 390"""
experiments = [
    ("laminar", None, False),
    ("disabledRpi", disabled, False),
    ("disabledLee", disabled.replace("RPI", "Lee"), False),
    ("embeddedDisabledRpi", "enabled: false\nmodel: RPI", True),
    ("leeEvaporation", lee.format(tsat=290), False),
    ("leeCondensation", lee.format(tsat=310), False),
    ("rpiActive", rpi, False),
    ("rpiBalanced", rpi.replace("convectiveHeatFluxModel: fixed\n  convectiveHeatFlux: 1000", "convectiveHeatFluxModel: coefficient\n  convectiveHeatTransferCoefficient: 50").replace("wallTemperatureModel: explicit", "wallTemperatureModel: heatFluxBalance"), False),
    ("rpiBalancedOuter3", rpi.replace("convectiveHeatFluxModel: fixed\n  convectiveHeatFlux: 1000", "convectiveHeatFluxModel: coefficient\n  convectiveHeatTransferCoefficient: 50").replace("wallTemperatureModel: explicit", "wallTemperatureModel: heatFluxBalance"), False),
    ("leeEvaporationOuter3", lee.format(tsat=290), False),
    ("leeEvaporationSourceOnly", lee.format(tsat=290), False),
    ("leeCondensationSourceOnly", lee.format(tsat=310), False),
    ("rpiBalancedSourceOnly", rpi.replace("convectiveHeatFluxModel: fixed\n  convectiveHeatFlux: 1000", "convectiveHeatFluxModel: coefficient\n  convectiveHeatTransferCoefficient: 50").replace("wallTemperatureModel: explicit", "wallTemperatureModel: heatFluxBalance"), False),
    ("rpiBalancedOuter3SourceOnly", rpi.replace("convectiveHeatFluxModel: fixed\n  convectiveHeatFlux: 1000", "convectiveHeatFluxModel: coefficient\n  convectiveHeatTransferCoefficient: 50").replace("wallTemperatureModel: explicit", "wallTemperatureModel: heatFluxBalance"), False),
]
for name, phase, embedded in experiments:
    if len(sys.argv)>1 and name not in sys.argv[1:]:
        continue
    for version in ("before", "after"):
        case = prepare(f"{name}-{version}-verified", phase, embedded)
        binary = evidence / "sonicSolver-before" if version == "before" else root / "build/sonicSolver"
        command = [str(binary), "run", str(case)]
        run = subprocess.run(command, cwd=case, text=True, capture_output=True,
            env={**os.environ, "TMPDIR": str(evidence), "DEVELOPER_DIR": "/Library/Developer/CommandLineTools", "SF_PLAN_TRACE": "1"})
        log = run.stdout+run.stderr
        (case / "run.log").write_text(log)
        files = {str(p.relative_to(case / "result")):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (case / "result").rglob("*") if p.is_file() and not p.name.startswith(".") and p.suffix != ".pvd"}
        record = dict(name=name, version=version, case=str(case), source=str(source), embedded=embedded,
            command=command, binary=hashlib.sha256(binary.read_bytes()).hexdigest(), exit=run.returncode,
            files=files, physical_lines=[line for line in log.splitlines() if "Time=" in line or "time=" in line or "totalWallHeat=" in line],
            error=log[-1000:] if run.returncode else "")
        manifest.append(record)
        (evidence / "cases.json").write_text(json.dumps(manifest, indent=2))
        print(name, version, run.returncode, len(files), record["error"][-400:], flush=True)

checks = []
for name, _, _ in experiments:
    if not any(x["name"]==name for x in manifest):
        continue
    before, after = [x for x in manifest if x["name"] == name]
    checks.append(dict(name=name, before_exit=before["exit"], after_exit=after["exit"],
        output_exact=before["exit"] == after["exit"] == 0 and before["files"] == after["files"],
        physical_log_exact=before["physical_lines"] == after["physical_lines"]))
laminar = next(x for x in manifest if x["name"] == "laminar" and x["version"] == "after")
for name in ("disabledRpi", "disabledLee", "embeddedDisabledRpi"):
    after = next(x for x in manifest if x["name"] == name and x["version"] == "after")
    checks.append(dict(name=name+"-vs-absent", output_exact=after["exit"] == 0 and after["files"] == laminar["files"],
        physical_log_exact=after["physical_lines"] == laminar["physical_lines"]))
(evidence / "comparisons.json").write_text(json.dumps(checks, indent=2))
print(json.dumps(checks, indent=2))
for check in checks:
    if check["name"]=="rpiActive":
        assert check["before_exit"]!=0 and check["after_exit"]!=0
        assert all("RPI wall ledger heat imbalance" in x["error"] for x in manifest if x["name"]=="rpiActive")
    elif check["name"] in ("disabledRpi","disabledLee","embeddedDisabledRpi"):
        assert check["before_exit"]!=0 and check["after_exit"]==0
    else:
        assert check["output_exact"] and check["physical_log_exact"],check
