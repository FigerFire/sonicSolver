#include "SF_declarations.h"
#include "core/config/types/SF_timeRecipe.h"
#include <cmath>
#include <sstream>
namespace SF::System::Declarations {
using P=Model::Parameters;
FormulaExpr formula(const P& node,const std::string& source,const FormulaSymbols& symbols,int depth) {
    using E=FormulaExpr;
    if(depth>64)throw std::runtime_error(source+": formula nesting exceeds 64.");
    if(!node.is_object())throw std::runtime_error(source+": typed formula node must be a mapping.");
    if(node.contains("symbol")) {
        Model::Schema{{{"symbol","string",true}},false}.validate(node,source);
        const auto name=node.at("symbol").get<std::string>();
        if(name.empty())throw std::runtime_error(source+": empty formula symbol.");
        return symbols.at(name).atOccurrence(source);
    }
    if(node.contains("constant")) {
        Model::Schema{{{"constant","number",true},{"dimensions","array"}},false}.validate(node,source);
        const auto value=node.at("constant").get<double>();
        if(!std::isfinite(value))throw std::runtime_error(source+": formula constant must be finite.");
        if(node.contains("dimensions")) {
            if(node.at("dimensions").size()!=7)throw std::runtime_error(source+": constant requires 7 dimension exponents");
            return E::constantValue(value,Math::signature(0,node.at("dimensions").get<std::array<double,7>>())).atOccurrence(source);
        }
        return E::constantValue(value).atOccurrence(source);
    }
    Model::Schema{{{"op","string",true},{"args","array",true},{"occurrence","string"}},false}.validate(node,source);
    const auto name=node.at("op").get<std::string>();
    std::vector<E> args;for(std::size_t i=0;i<node.at("args").size();++i)args.push_back(formula(node.at("args")[i],source+"/args/"+std::to_string(i),symbols,depth+1));
    return E::op(Math::resolveOperator(name),std::move(args)).atOccurrence(node.value("occurrence",source));
}

ExecutionProgram execution(const P& parameters,int formatVersion,Provenance origin) {
    Model::Schema{{{"root","any",true},{"residuals","array"}},false}.validate(parameters,origin.source);
    const auto parse=[&](const auto& self,const P& value)->ExecutionScope {
        if(value.is_array()) {
            ExecutionScope sequence;
            for(const auto& child:value)sequence.children.push_back(self(self,child));
            return sequence;
        }
        if(value.is_string()) {
            const auto text=value.get<std::string>();
            const auto trim=[](std::string s){auto first=s.find_first_not_of(" \t\r\n");return first==std::string::npos?std::string():s.substr(first,s.find_last_not_of(" \t\r\n")-first+1);};
            const bool block=text.find("=>")!=std::string::npos;const std::string arrow=block?"=>":"->";const auto pos=text.find(arrow);
            if(pos==std::string::npos || text.find(arrow,pos+2)!=std::string::npos)throw std::runtime_error("HOW binding requires one equation -> target or [equations] => [unknowns]");
            if(!block)return self(self,P{{"kind","EquationCall"},{"equation",trim(text.substr(0,pos))},{"target",trim(text.substr(pos+2))}});
            const auto list=[&](std::string part){part=trim(part);if(part.size()<2 || part.front()!='[' || part.back()!=']')throw std::runtime_error("Implicit block binding requires two explicit lists");auto body=trim(part.substr(1,part.size()-2));if(body.empty() || body.front()==',' || body.back()==',')throw std::runtime_error("Empty implicit block member");std::vector<std::string> out;std::istringstream input(body);std::string id;while(std::getline(input,id,',')){id=trim(id);if(id.empty())throw std::runtime_error("Empty implicit block member");out.push_back(id);}if(out.empty())throw std::runtime_error("Empty implicit block binding");return out;};
            return self(self,P{{"kind","BlockSolve"},{"equations",list(text.substr(0,pos))},{"unknowns",list(text.substr(pos+2))}});
        }
        Model::Schema{{{"kind","string"},{"id","string"},{"order","integer"},
            {"equation","string"},{"target","string"},{"targetKind","string"},{"occurrence","string"},
            {"children","array"},{"repetitions","integer"},{"minimumIterations","integer"},
            {"terminationSignal","string"},{"equations","array"},{"unknowns","array"},
            {"residuals","array"},{"residualAll","boolean"},{"failOnMaximum","boolean"},
            {"sources","array"},{"targets","array"},{"frequency","string"}},false}.validate(value,"HOW node");
        ExecutionScope node;
        const auto kind=value.value("kind",std::string(value.contains("equations")?"BlockSolve":value.contains("equation")?"EquationCall":"Sequence"));
        if (kind=="Sequence") node.kind=ExecutionKind::Sequence;
        else if (kind=="Loop") node.kind=ExecutionKind::Loop;
        else if (kind=="StageLoop") {
            if(formatVersion>=2)throw std::runtime_error("Public HOW StageLoop is retired in formatVersion: 2; use an ordered body and select the time recipe in WHICH/numerics.");
            node.kind=ExecutionKind::StageLoop;
        }
        else if (kind=="Commit") node.kind=ExecutionKind::Commit;
        else if (kind=="EquationCall") node.kind=ExecutionKind::EquationCall;
        else if (kind=="BlockSolve") node.kind=ExecutionKind::BlockSolve;
        else throw std::runtime_error("Unknown HOW node kind: "+kind);
        node.id=value.value("id",std::string());node.order=value.value("order",0);
        if(value.contains("frequency")) {
            if(formatVersion<2)throw std::runtime_error("HOW frequency requires formatVersion: 2; historical version 1 declares StageLoop explicitly");
            const auto frequency=value.at("frequency").get<std::string>();
            if(node.kind!=ExecutionKind::Sequence && node.kind!=ExecutionKind::Loop && node.kind!=ExecutionKind::EquationCall && node.kind!=ExecutionKind::BlockSolve)
                throw std::runtime_error("HOW frequency applies to a mathematical call/group, not a technical/Commit node");
            if(frequency=="Step")node.frequency=ExecutionFrequency::Step;
            else if(frequency=="Stage")node.frequency=ExecutionFrequency::Stage;
            else throw std::runtime_error("Unknown HOW frequency: "+frequency+"; select Step or Stage");
        }
        node.origin=origin;
        node.repetitions=value.value("repetitions",1);node.minimumIterations=value.value("minimumIterations",1);
        if(value.contains("terminationSignal"))throw std::runtime_error("Native HOW terminationSignal is an internal provider contract; use typed residuals for convergence.");
        if(value.contains("residuals") || value.contains("residualAll") || value.contains("failOnMaximum")) {
            if(node.kind!=ExecutionKind::Loop)throw std::runtime_error("HOW residual controls require a Loop scope.");
            node.residuals=value.value("residuals",std::vector<std::string>{});
            node.residualAll=value.value("residualAll",true);node.failOnMaximum=value.value("failOnMaximum",true);
        }
        if(value.contains("sources") || value.contains("targets")) {
            if(node.kind!=ExecutionKind::Commit || !value.contains("sources") || !value.contains("targets")
                || value.contains("children"))throw std::runtime_error("HOW variable Commit requires explicit sources/targets and no child scope.");
            for(const auto& source:value.at("sources"))node.commitSources.push_back(targetFromSyntax(source.get<std::string>()));
            for(const auto& target:value.at("targets"))node.commitTargets.push_back(targetFromSyntax(target.get<std::string>()));
            if(node.commitSources.empty() || node.commitSources.size()!=node.commitTargets.size())
                throw std::runtime_error("HOW variable Commit requires equally sized nonempty source/target lists.");
        }
        if (node.kind==ExecutionKind::BlockSolve) {
            if(!value.contains("equations") || !value.contains("unknowns") || value.contains("equation") || value.contains("target") || value.contains("children"))throw std::runtime_error("HOW BlockSolve requires equation/unknown sets, not sequential children");
            node.equations=value.at("equations").get<std::vector<std::string>>();node.unknowns=value.at("unknowns").get<std::vector<std::string>>();
        } else if(value.contains("equations") || value.contains("unknowns"))throw std::runtime_error("Block sets require a BlockSolve scope");
        if (node.kind==ExecutionKind::EquationCall) {
            if (!value.contains("equation") || !value.contains("target") || value.contains("children"))
                throw std::runtime_error("HOW EquationCall requires equation/target and no child scope.");
            node.step={value.at("equation").get<std::string>(),
                targetFromSyntax(value.at("target").get<std::string>()),value.value("occurrence",std::string())};
            if(!node.id.empty() && !node.step.occurrence.empty() && node.id!=node.step.occurrence)
                throw std::runtime_error("HOW has two different Step IDs: id and occurrence.");
            if(node.step.occurrence.empty())node.step.occurrence=node.id;
            if (value.contains("targetKind")) {
                if (node.step.target.kind!=TargetKind::Physical)
                    throw std::runtime_error("HOW targetKind cannot override a qualified target.");
                const auto targetKind=value.at("targetKind").get<std::string>();
                if (targetKind=="Physical") node.step.target.kind=TargetKind::Physical;
                else if (targetKind=="Working") node.step.target.kind=TargetKind::Working;
                else if (targetKind=="Correction") node.step.target.kind=TargetKind::Correction;
                else if (targetKind=="Workspace") node.step.target.kind=TargetKind::Workspace;
                else throw std::runtime_error("Unknown HOW target kind: "+targetKind);
            }
        } else {
            if (value.contains("equation") || value.contains("target") || value.contains("targetKind") || value.contains("occurrence"))
                throw std::runtime_error("HOW control scope cannot carry an equation target.");
            if (value.contains("children")) for (const auto& child:value.at("children")) {
                node.children.push_back(self(self,child));
            }
        }
        return node;
    };
    ExecutionProgram program;program.root=parse(parse,parameters.at("root"));
    program.recipeScopes=formatVersion>=2;
    if(parameters.contains("residuals"))for(const auto& entry:parameters.at("residuals")) {
        Model::Schema{{{"id","string",true},{"occurrence","string",true},{"kind","string"},{"norm","string"},
            {"normalization","string"},{"absoluteTolerance","number"},{"relativeTolerance","number"},{"scale","number"}},false}.validate(entry,"HOW residual definition");
        ResidualDefinition residual;residual.id=entry.at("id").get<std::string>();
        residual.occurrence=entry.at("occurrence").get<std::string>();
        const auto kind=entry.value("kind",std::string("StateChange"));
        if(kind=="StateChange")residual.kind=ResidualKind::StateChange;
        else if(kind=="DiscreteEquation")residual.kind=ResidualKind::DiscreteEquation;
        else if(kind=="LinearSolver")residual.kind=ResidualKind::LinearSolver;
        else throw std::runtime_error("Unknown residual metric kind: "+kind);
        const auto norm=entry.value("norm",std::string("Maximum"));
        if(norm=="Maximum")residual.norm=ResidualNorm::Maximum;
        else if(norm=="L2")residual.norm=ResidualNorm::L2;
        else throw std::runtime_error("Unknown residual norm: "+norm);
        const auto normalization=entry.value("normalization",std::string("Initial"));
        if(normalization=="Initial")residual.normalization=ResidualNormalization::Initial;
        else if(normalization=="Scale")residual.normalization=ResidualNormalization::Scale;
        else if(normalization=="None")residual.normalization=ResidualNormalization::None;
        else throw std::runtime_error("Unknown residual normalization: "+normalization);
        residual.absoluteTolerance=entry.value("absoluteTolerance",1e-10);residual.relativeTolerance=entry.value("relativeTolerance",1e-8);
        residual.scale=entry.value("scale",1.0);program.residuals.push_back(std::move(residual));
    }
    if (program.root.kind!=ExecutionKind::Sequence || program.root.children.empty())
        throw std::runtime_error("Authored HOW root requires a nonempty Sequence.");
    return program;
}

std::vector<NumericalBinding> numerics(const P& parameters,Provenance origin) {
    Model::Schema{{{"bindings","array",true}},false}.validate(parameters,origin.source);
    std::vector<NumericalBinding> result;
    for (const auto& value:parameters.at("bindings")) {
        if(value.contains("reads"))throw std::runtime_error("WHICH input versions were retired; declare STATE calls/reads for the execution occurrence.");
        Model::Schema{{{"equation","string"},{"step","string"},{"method","string"},{"inputs","array"},
            {"occurrence","string"},{"parameters","object"},{"operators","array"},{"time","string"}},false}.validate(value,"WHICH binding");
        if(!value.contains("equation") && !value.contains("step"))throw std::runtime_error("WHICH requires an explicit Equation or Step reference.");
        if(value.contains("step") && value.contains("occurrence"))throw std::runtime_error("WHICH duplicate Step reference: step and occurrence.");
        NumericalBinding binding;
        binding.origin=origin;
        binding.equation=value.value("equation",std::string("*"));binding.method=value.value("method",std::string());
        binding.occurrence=value.value(value.contains("step")?"step":"occurrence",std::string());
        binding.inputsDeclared=value.contains("inputs");
        if(value.contains("operators"))for(const auto& op:value.at("operators")) {
            Model::Schema{{{"occurrence","string",true},{"provider","string",true},{"equation","string"}},false}.validate(op,"WHICH operator binding");
            binding.operators.emplace_back(op.value("equation",binding.equation),op.at("occurrence").get<std::string>(),op.at("provider").get<std::string>());
        }
        if(value.contains("time")) {
            bool known=false;
            for(auto id:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::SSPRK3,FDM::TimeRecipeId::ClassicalRK4})
                if(value.at("time").get<std::string>()==FDM::toString(id)){binding.time=id;known=true;break;}
            if(!known)throw std::runtime_error("Unknown WHICH time recipe: "+value.at("time").get<std::string>());
        }
        if (value.contains("inputs")) binding.inputs=value.at("inputs").get<std::vector<std::string>>();
        if (value.contains("parameters")) for (const auto& entry:value.at("parameters").items()) {
            if (!entry.value().is_number()) throw std::runtime_error("WHICH parameters must be numeric: "+entry.key());
            binding.parameters.emplace(entry.key(),entry.value().get<double>());
        }
        result.push_back(std::move(binding));
    }
    return result;
}

}
