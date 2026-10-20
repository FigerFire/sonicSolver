#!/usr/bin/env python3
"""Render boundary mappings with inline first conditions for readability."""
import json
import subprocess
from pathlib import Path

def load(path):
    out=subprocess.run(["ruby","-rjson","-ryaml","-e",
        "puts JSON.generate(YAML.load_file(ARGV[0]))",str(path)],check=True,capture_output=True,text=True)
    return json.loads(out.stdout)

def value(v):
    if isinstance(v,list): return "["+", ".join(value(x) for x in v)+"]"
    if isinstance(v,str): return v if v not in {"true","false","null"} else json.dumps(v)
    if v is True:return "true"
    if v is False:return "false"
    if v is None:return "null"
    return str(v)

def dump(v, level=0):
    indent="  "*level; lines=[]
    for k, child in v.items():
        if isinstance(child,dict) and len(child)==1:
            ck,cv=next(iter(child.items()))
            lines.append(f"{indent}{k}: {{{ck}: {value(cv)}}}")
        elif isinstance(child,dict):
            lines.append(f"{indent}{k}:");lines.extend(dump(child,level+1).splitlines())
        else: lines.append(f"{indent}{k}: {value(child)}")
    return "\n".join(lines)

for path in Path("test").rglob("fields/boundaries.yaml"):
    if "ioRegistry" in path.parts: continue
    data=load(path)
    path.write_text("# Boundary conditions: patch: {fixedValue: value} or patch: zeroGradient.\n"+dump(data)+"\n")
