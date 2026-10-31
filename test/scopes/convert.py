"""Explicit copy-only migration using the recorded old Stage scopes, never a reader fallback."""
from pathlib import Path
import json
STAGES={'how-synchronous','how-synchronousReverse','how-stageInnerLoop','how-stageInnerConvergence','how-stageInitialGuess'}
def convert(case):
 case=Path(case)
 if case.name not in STAGES:return
 file=case/'solvers/how.yaml';v=json.loads(file.read_text())
 def mark(node):
  if isinstance(node,list):
   for c in node:mark(c)
  elif node.get('kind')!='Commit':node['frequency']='Stage'
 mark(v['root']);file.write_text(json.dumps(v,indent=2)+'\n')
 (case/'scope-migration.json').write_text(json.dumps({'source':'frozen lowering Stage-scoped cases','change':'explicit Stage frequency; same mathematical timing'},indent=2)+'\n')
