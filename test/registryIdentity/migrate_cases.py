"""One-time explicit fixture migration; never used by the production reader."""
from pathlib import Path
import re, subprocess, os, json, shutil, sys
ROOT=Path(__file__).resolve().parents[2]
ROLE_NAMES={'stateRegistry':'state','equationRegistry':'equationRegistry','algorithmRegistry':'execution'}
def entries(text):
    starts=list(re.finditer(r'^([A-Za-z_][\w.-]*):([^\n]*)$',text,re.M))
    for i,m in enumerate(starts):
        if m.group(1)=='SonicFile':continue
        yield m.group(1),m.group(2).strip(),text[m.end():starts[i+1].start() if i+1<len(starts) else len(text)]
def title(text,obj,name):
    if text.lstrip().startswith("{"):
        data=json.loads(text);data["SonicFile"]={"object":obj,"type":name}
        return json.dumps(data,indent=2)+"\n"
    text=re.sub(r'^SonicFile:\s*\{object:\s*[^,}]+,\s*type:\s*[^}]+\}\s*$', 'SonicFile:\n  object: '+obj+'\n  type: '+name,text,flags=re.M)
    text,n=re.subn(r'^(  object:).*$',r'\1 '+obj,text,count=1,flags=re.M)
    if n!=1:raise ValueError('Missing explicit SonicFile.object')
    text,n=re.subn(r'^(  type:).*$',r'\1 '+name,text,count=1,flags=re.M)
    if n!=1:raise ValueError('Missing explicit SonicFile.type')
    return text
def migrate(case):
    changes=[]
    model_entries=[]
    for obj,folder,registry in [('solver','solvers','solvers.yaml'),('models','models','models.yaml')]:
        p=case/folder/registry;text=p.read_text();result='# Explicit document registry.\nSonicFile:\n  object: '+obj+'\n  type: registry\n'
        rows=[]
        for name,compact,block in [*entries(text), *(model_entries if obj=="models" else [])]:
            role=re.search(r'^  type:\s*(.*?)\s*$',block,re.M);file=re.search(r'^  file:\s*(.*?)\s*$',block,re.M)
            role=role.group(1).strip('"\'') if role else name
            source=file.group(1).strip('"\'') if file else folder+'/'+name+'.yaml'
            new_name=ROLE_NAMES.get(role, 'equations' if obj=='solver' and name=='flow' and role in ('equationSystem','BuiltinEquationSystem') else name)
            if obj=='solver' and role=='IBM':
                if file is None:source='models/IBM.yaml'
                model_entries.append((name,compact,'\n  type: IBM\n  file: '+source+'\n'))
                continue
            if obj=='solver' and role=='ILW' and file is None:source='models/ILW.yaml'
            old=case/source
            if not old.is_file():raise ValueError(str(old)+' is a registered missing file')
            dest=case/folder/(new_name+'.yaml')
            child=title(old.read_text(),obj,new_name)
            rows.append((new_name,role,dest.relative_to(case).as_posix(),compact if not file else ''))
            if dest!=old and dest.exists() and title(dest.read_text(),obj,new_name)!=child:
                raise ValueError('Conflicting registered output '+str(dest))
            dest.parent.mkdir(exist_ok=True);dest.write_text(child)
            if dest!=old:old.unlink();changes.append({'from':source,'to':dest.relative_to(case).as_posix()})
        # Explicitly migrate formerly auto-loaded legacy root selections once.
        if obj=='solver':
            for legacy,role in [('state/state.yaml','stateRegistry'),('equations/equations.yaml','equationRegistry'),('algorithms/algorithms.yaml','algorithmRegistry')]:
                q=case/legacy
                if q.exists():
                    name=ROLE_NAMES[role]
                    if any(row[0]==name for row in rows):raise ValueError('Duplicate declared selection '+str(case)+'/'+name)
                    dest=case/folder/(name+'.yaml');dest.write_text(title(q.read_text(),obj,name));q.unlink();rows.append((name,role,dest.relative_to(case).as_posix(),''));changes.append({'from':legacy,'to':dest.relative_to(case).as_posix()})
                    if not list(q.parent.iterdir()):q.parent.rmdir()
        names=[r[0] for r in rows]
        if len(names)!=len(set(names)):raise ValueError('Duplicate registry names '+str(case))
        for name,role,file,selection in rows:
            result+=f'{name}:\n  type: {role}\n  file: {file}\n'
            if selection:result+='  selection: '+selection+'\n'
        p.write_text(result)
    fields=case/'fields/fields.yaml';s=fields.read_text()
    for name,role in [('internalField','internalField'),('boundaries','boundary')]:
        child=case/'fields'/(name+'.yaml')
        if child.exists():
            child.write_text(title(child.read_text(),'fields',name))
            if not re.search(r'^'+name+r':',s,re.M):s+=f'\n{name}:\n  type: {role}\n  file: fields/{name}.yaml\n'
    fields.write_text(s)
    return {'case':str(case.relative_to(ROOT)),'moved':changes}
if __name__=='__main__':
    if sys.argv[1:] == ['--tracked']:
        paths=subprocess.check_output(['git','ls-files','*case.yaml'],cwd=ROOT,env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools'},text=True).splitlines()
        cases=[ROOT/p for p in paths if p.startswith('test/') and (ROOT/p).is_file()]
    else:cases=[ROOT/p/'case.yaml' for p in sys.argv[1:]]
    results=[]
    for case in cases:
        print(str(case.parent),file=sys.stderr)
        results.append(migrate(case.parent))
    print(json.dumps(results,indent=2))
