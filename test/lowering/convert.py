"""Explicit conversion of this regression's frozen A9 copies; never a production reader."""
from pathlib import Path
import json,copy

def convert(case):
    case=Path(case);read=lambda p:json.loads(p.read_text());save=lambda p,v:p.write_text(json.dumps(v,indent=2)+'\n')
    if not (case/'solvers/how.yaml').exists():return
    meta=read(case/'case.yaml')
    if meta['formatVersion']==2:return
    changed=[]
    def update(path,value):save(case/path,value);changed.append(path)
    howPath=case/'solvers/how.yaml'
    if howPath.exists():
        how=read(howPath);equations={};occurrences=set()
        def address(node,parent='',index=0):
            if isinstance(node,str):
                left,right=node.split('=>' if '=>' in node else '->')
                node={'kind':'BlockSolve','equations':[v.strip() for v in left.strip()[1:-1].split(',')],'unknowns':[v.strip() for v in right.strip()[1:-1].split(',')]} if '=>' in node else {'kind':'EquationCall','equation':left.strip(),'target':right.strip()}
            node=copy.deepcopy(node);node.setdefault('id',parent+'/'+str(index))
            if node.get('kind')=='BlockSolve':occurrences.add(node['id'])
            if node.get('kind')=='EquationCall':
                node.setdefault('occurrence',node['id']);node['id']=node['occurrence'];occurrences.add(node['occurrence']);equations.setdefault(node['equation'],[]).append(node['occurrence'])
            node['children']=[address(n,node['id'],i) for i,n in enumerate(node['children'])] if 'children' in node else node.get('children',None)
            if node.get('children') is None:node.pop('children',None)
            return node
        root=address(how['root']);assert root.get('repetitions',1)==1
        state=read(case/'solvers/state.yaml')
        for binding in state.get('calls',[]):
            if 'call' in binding and binding['call'] not in occurrences:binding['equation']=binding.pop('call')
        # Old x' correction becomes a distinct declared unknown and Working result.
        if case.name=='how-correctionInput':
            fields=read(case/'fields/fields.yaml');fields['deltaX']={'new':True,'type':'scalar','output':False,'dimensions':[0]*7};update('fields/fields.yaml',fields)
            initial=read(case/'fields/internalField.yaml');initial['deltaX']={'default':0.};update('fields/internalField.yaml',initial)
            state['use'].append('deltaX');what=read(case/'solvers/equationRegistry.yaml');what['use'].append('PredictX')
            def rename(e):
                if e.get('symbol')=='x':e['symbol']='deltaX'
                for a in e.get('args',[]):rename(a)
            what['add']['PredictX']=copy.deepcopy(what['add']['X']);rename(what['add']['PredictX']['lhs']);rename(what['add']['PredictX']['rhs']);rename(what['add']['Y']['lhs']);rename(what['add']['Y']['rhs']);update('solvers/equationRegistry.yaml',what)
            def correction(node):
                if node.get('occurrence')=='predictX':node.update(equation='PredictX',target='deltaX*')
                for c in node.get('children',[]):correction(c)
            correction(root)
            for binding in state.get('calls',[]):
                for item in binding['reads']:
                    if item['symbol']=="x'":item['symbol']='deltaX*'
            which=read(case/'solvers/which.yaml');which['bindings'].append({'equation':'PredictX','method':'LocalAffine'});update('solvers/which.yaml',which)
        def flatten(node):
            children=[]
            for child in node.get('children',[]):
                if child.get('kind')=='StageLoop':children.extend(flatten(c) for c in child['children'])
                else:children.append(flatten(child))
            if 'children' in node:node['children']=children
            return node
        root=flatten(root);how['root']=root['children'];update('solvers/how.yaml',how);update('solvers/state.yaml',state)
    meta['formatVersion']=2;update('case.yaml',meta)
    save(case/'migration.json',{'sourceContract':'A9 frozen input','targetContract':'formatVersion 2; recipe-owned HOW; typed STATE refs','changed':changed,'primeCorrectionConverted':case.name=='how-correctionInput'})
