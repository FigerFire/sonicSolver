#include "solver/run/SF_planExecutor.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
using namespace SF;using namespace System;using namespace Run;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f,const std::string& text){try{f();}catch(const std::exception& e){if(std::string(e.what()).find(text)!=std::string::npos)return;throw;}throw std::runtime_error("Expected failure: "+text);}
SolvePlanNode leaf(std::string id){SolvePlanNode n;n.id=n.operation=id;n.kind=PlanNodeKind::Update;n.provider="test.method";return n;}
ResidualDefinition metric(std::string id,std::string op){ResidualDefinition r;r.id=id;r.occurrence=op;r.provider="test.method";r.producerOperation=op;r.absoluteTolerance=1e-8;r.relativeTolerance=0;r.normalization=ResidualNormalization::None;return r;}
CompiledSolvePlan plan(SolvePlanNode root){CompiledSolvePlan p;p.root=std::move(root);return p;}
int main(){try{
    OpRegistry ops;double x=0;int calls=0;auto r=metric("dx","relax");
    ops.bind("relax","test.method",[&]{x=1+.5*x;++calls;});
    ops.bindResidual("dx","test.method",[&](const auto&){return std::vector<double>{x};});
    SolvePlanNode inner;inner.id="inner";inner.kind=PlanNodeKind::Loop;inner.minimumIterations=2;inner.repetitions=100;inner.residuals={r};inner.children={leaf("relax")};
    SolvePlanNode outer;outer.id="outer";outer.kind=PlanNodeKind::Loop;outer.repetitions=3;outer.children={inner};
    PlanExecutor::execute(plan(outer),ops);require(calls==32,"inner exit/reset affected outer or min iterations");
    x=0;calls=0;inner.repetitions=2;rejects([&]{PlanExecutor::execute(plan(inner),ops);},"maximum iterations");require(calls==2,"maximum was not honored");
    inner.failOnMaximum=false;PlanExecutor::execute(plan(inner),ops);require(calls==4,"explicit continue policy failed");
    auto stale=inner;stale.children={leaf("noop")};ops.bind("noop","test.method",[]{});
    rejects([&]{PlanExecutor::execute(plan(stale),ops);},"was not updated");
    auto missing=inner;missing.residuals.front().id="unbound";rejects([&]{PlanExecutor::execute(plan(missing),ops);},"no sampler");
    OpRegistry signals;int outerCalls=0,innerCalls=0;
    signals.bind("innerDone","test.method",[&](const ExecutionContext& c){++innerCalls;c.signals->publish("inside",true);});
    signals.bind("outerDone","test.method",[&](const ExecutionContext& c){c.signals->publish("outside",++outerCalls==3);});
    SolvePlanNode in;in.kind=PlanNodeKind::Loop;in.id="insideLoop";in.repetitions=8;in.terminationSignal="inside";in.children={leaf("innerDone")};
    auto out=in;out.id="outsideLoop";out.terminationSignal="outside";out.children={in,leaf("outerDone")};
    PlanExecutor::execute(plan(out),signals);require(outerCalls==3 && innerCalls==3,"nested signal exited the wrong scope");
    OpRegistry old;int iterations=0;old.bind("once","test.method",[&](const ExecutionContext& c){if(++iterations==1)c.signals->publish("stale",false);});
    auto staleSignal=in;staleSignal.terminationSignal="stale";staleSignal.children={leaf("once")};
    rejects([&]{PlanExecutor::execute(plan(staleSignal),old);},"was not published");
    // Criteria stay independent; norm, scale normalization and AND/OR are algorithm controls.
    for(bool all:{false,true}) {
        OpRegistry criteria;double value=1;int n=0;
        criteria.bind("halve","test.method",[&]{value*=.5;++n;});
        for(const auto* id:{"loose","tight"})criteria.bindResidual(id,"test.method",[&](const auto&){return std::vector<double>{value,value};});
        auto loose=metric("loose","halve"),tight=metric("tight","halve");
        loose.kind=tight.kind=ResidualKind::DiscreteEquation;loose.absoluteTolerance=.2;
        tight.absoluteTolerance=0;tight.relativeTolerance=.02;tight.normalization=ResidualNormalization::Scale;tight.scale=1;tight.norm=ResidualNorm::L2;
        auto loop=inner;loop.minimumIterations=1;loop.repetitions=20;loop.failOnMaximum=true;loop.residualAll=all;loop.residuals={loose,tight};loop.children={leaf("halve")};
        PlanExecutor::execute(plan(loop),criteria);require(n==(all?7:3),"AND/OR or L2 scale criterion changed meaning");
    }
    OpRegistry references;double v=1.;int total=0;
    references.bind("halve","test.method",[&]{v*=.5;++total;});
    references.bindResidual("relative","test.method",[&](const auto&){return std::vector<double>{v};});
    auto relative=metric("relative","halve");relative.absoluteTolerance=0;relative.relativeTolerance=.2;relative.normalization=ResidualNormalization::Initial;
    auto resetInner=inner;resetInner.minimumIterations=1;resetInner.repetitions=20;resetInner.residuals={relative};resetInner.children={leaf("halve")};
    auto resetOuter=outer;resetOuter.repetitions=2;resetOuter.children={resetInner};PlanExecutor::execute(plan(resetOuter),references);
    require(total==8,"relative reference did not reset on inner-loop entry");
    OpRegistry invalid;int mutated=0;invalid.bind("relax","test.method",[&]{++mutated;});
    invalid.bindResidual("dx","test.method",[](const auto&){return std::vector<double>{NAN};});
    auto invalidReference=inner;invalidReference.failOnMaximum=true;
    rejects([&]{PlanExecutor::execute(plan(invalidReference),invalid);},"Invalid initial residual reference");require(mutated==0,"invalid convergence reference executed numerical body");
    // One recipe stage encloses an ordinary iteration; no callback owns a stage loop.
    OpRegistry stages;std::vector<int> observed;stages.bind("sample","test.method",[&](const ExecutionContext& c){observed.push_back(c.stageIndex);require(c.stageCount==4,"nested loop changed recipe stages");});
    SolvePlanNode fixed;fixed.kind=PlanNodeKind::Loop;fixed.id="fixed";fixed.repetitions=2;fixed.children={leaf("sample")};
    SolvePlanNode stage;stage.kind=PlanNodeKind::StageLoop;stage.id="stages";stage.repetitions=4;stage.children={fixed};PlanExecutor::execute(plan(stage),stages);
    require(observed==std::vector<int>({0,0,1,1,2,2,3,3}),"stage/loop topology was flattened");
    auto sub=stage;sub.kind=PlanNodeKind::Subcycle;rejects([&]{PlanExecutor::execute(plan(sub),stages);},"does not implement");
    std::cout<<"A9 residual freshness, nested loops, maximum policy, stage context and unsupported subcycle passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
