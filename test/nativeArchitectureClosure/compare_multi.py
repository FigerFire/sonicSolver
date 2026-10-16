from pathlib import Path
import csv,json,hashlib,re
B=Path(__file__).resolve().parent
out=[]
for a in sorted((B/'d2MultiBefore').glob('*/*')):
 if not (a/'run.log').exists() or not (a/'result').exists() or 'serial composite solve is not implemented' in (a/'run.log').read_text():continue
 b=B/'d2MultiAfter'/a.relative_to(B/'d2MultiBefore')
 def results(p):return {str(f.relative_to(p/'result')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'result').rglob('*') if f.is_file() and not f.name.startswith('._')}
 assert results(a)==results(b),('output first difference',a)
 def checkpoints(p):
  result={};pidRank={}
  for f in (p/'checkpoint').glob('stage-*.csv'):
   if f.name.startswith('._'):continue
   with f.open() as stream:row=next(csv.DictReader(stream))
   pidRank[f.stem.split('-')[1]]=row['rank']
  for f in (p/'checkpoint').glob('*.csv'):
   if f.name.startswith('._'):continue
   rows=list(csv.DictReader(f.open()))
   if f.name.startswith('levelSet-'):
    key=(rows[0]['rank'],rows[0]['sequence'],rows[0]['operation'],tuple((rows[0][v] for v in ('x','y','z'))))
   else:
    if f.name.startswith('stage-'):
     for row in rows:row.pop('patch',None)
     key=('stage',rows[0]['rank'],f.stem.split('-')[-1])
    else:key=('runtime',pidRank[f.stem.split('-')[1]])
   assert key not in result,(p,key)
   result[key]=rows
  return result
 x,y=checkpoints(a),checkpoints(b)
 assert x==y,('checkpoint first difference',a,set(x)^set(y))
 out.append(dict(case=str(a.relative_to(B/'d2MultiBefore')),result_files=len(results(a)),checkpoints=len(x),exact=True))
(B/'evidence/d2-multi-exact-comparison.json').write_text(json.dumps(out,indent=2))
print(len(out),'pairs: all result bytes and all checkpoint rows identical')
