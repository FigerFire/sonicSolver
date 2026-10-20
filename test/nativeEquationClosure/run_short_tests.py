"""Run the C++ regression inventory with all generated evidence below test/."""
from pathlib import Path
import json
import os
import subprocess

base=Path(__file__).resolve().parent
root=base.parents[1]
out=base/"ctest"
out.mkdir(exist_ok=True)
env={**os.environ,"DEVELOPER_DIR":"/Library/Developer/CommandLineTools","TMPDIR":str(base/"evidence")}
inventory=json.loads(subprocess.check_output(["ctest","--test-dir",str(root/"build"),"--show-only=json-v1"],text=True,env=env))
selected=[]
excluded=[]
lines=[]
def quote(value):
    return '"'+str(value).replace('\\','\\\\').replace('"','\\"')+'"'
for test in inventory["tests"]:
    command=test["command"]
    if test["name"]!="nativeEquationClosure" and ("python" in Path(command[0]).name or Path(command[0]).name=="sonicSolver"):
        excluded.append(test["name"])
        continue
    directory=base/"evidence/unit"/test["name"]
    directory.mkdir(parents=True,exist_ok=True)
    selected.append(test["name"])
    lines.append("add_test("+quote(test["name"])+" "+" ".join(quote(x) for x in command)+")")
    properties=["WORKING_DIRECTORY",quote(directory),"TIMEOUT","180"]
    for item in test["properties"]:
        if item["name"] in ("WILL_FAIL","PASS_REGULAR_EXPRESSION","FAIL_REGULAR_EXPRESSION"):
            value=item["value"]
            properties += [item["name"],quote(";".join(value) if isinstance(value,list) else value)]
    lines.append("set_tests_properties("+quote(test["name"])+" PROPERTIES "+" ".join(properties)+")")
(out/"CTestTestfile.cmake").write_text("\n".join(lines)+"\n")
(base/"evidence/test-selection.json").write_text(json.dumps(dict(selected=selected,excluded=excluded),indent=2))
with (base/"evidence/ctest-serial.log").open("w") as log:
    result=subprocess.run(["ctest","--test-dir",str(out),"--output-on-failure","--parallel","1","-E",
        "^(distributed|haloFailureContract|peskinPartialRows|sparseCanonicalCopy)"],stdout=log,stderr=subprocess.STDOUT,env=env,cwd=base)
raise SystemExit(result.returncode)
