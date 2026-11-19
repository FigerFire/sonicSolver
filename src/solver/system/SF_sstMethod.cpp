#include "SF_sstMethod.h"
#include "core/interfaces/SF_scalarEquationData.h"
#include <algorithm>
#include <cmath>
namespace SF::System {
CompiledEquationCall compileSstClosure(const ExecutableEquationSystem& system,const EquationCall& call,
        const NumericalBinding& binding) {
    if((!call.stageScoped && !call.stepScoped) || call.target.kind!=TargetKind::Physical || !binding.operators.empty())
        throw std::runtime_error("SSTCentral2 requires an explicit nonintegrated Step/Stage call with Central2 inputs: "+call.occurrence);
    const std::array<std::string,5> operations={"sst.F1","sst.eddyViscosity","sst.kProduction","sst.omegaProduction","sst.crossDiffusion"};
    CompiledSstClosure contract;
    CompiledEquationCall result;result.source=call;result.target.symbol=call.target.symbol;result.stageLocal=call.stageScoped;
    const auto equations=binding.inputs.empty()?std::vector<std::string>{call.equation}:binding.inputs;
    const auto& first=system.registry.at(equations.front());
    contract.eddyViscosity=first.rhs.name()==operations[1];
    const auto& args=first.rhs.arguments();
    if(args.size()!=(contract.eddyViscosity?7:11))throw std::runtime_error("SST operands require explicit density/velocity/scalars/constants and, for production, published viscosity: "+call.occurrence);
    for(std::size_t i=0;i<args.size();++i)if(i!=1 && args[i].kind()!=FormulaExpr::Kind::Symbol)
        throw std::runtime_error("SST operands require registered field/constant bindings: "+call.occurrence);
    const std::array<Math::Signature,6> inputSignatures={Math::signature(0,{1,-3,0,0,0,0,0}),Math::signature(1,{0,1,-1,0,0,0,0}),Math::signature(0,{0,2,-2,0,0,0,0}),Math::signature(0,{0,0,-1,0,0,0,0}),Math::signature(0,{1,-1,-1,0,0,0,0}),Math::signature(0,{0,0,-1,0,0,0,0})};
    for(std::size_t i=0;i<(contract.eddyViscosity?7:10);++i)Math::requireCompatible(args[i].signature(),i<inputSignatures.size() && !(contract.eddyViscosity && i>=5)?inputSignatures[i]:Math::signature(0,{}),"SST operand "+std::to_string(i));
    const std::array<Math::Signature,5> outputSignatures={Math::signature(0,{}),Math::signature(0,{1,-1,-1,0,0,0,0}),Math::signature(0,{1,-1,-3,0,0,0,0}),Math::signature(0,{1,-3,-2,0,0,0,0}),Math::signature(0,{1,-3,-2,0,0,0,0})};
    for(const auto& name:equations) {
        const auto& equation=system.registry.at(name);
        const auto found=std::find(operations.begin(),operations.end(),equation.rhs.name());
        if(found==operations.end())throw std::runtime_error("SST implementation does not implement operator: "+equation.rhs.name());
        const int role=found-operations.begin();
        if((role==1)!=contract.eddyViscosity || equation.lhs.kind()!=FormulaExpr::Kind::Symbol || equation.rhs.kind()!=FormulaExpr::Kind::Operator
            || equation.rhs.arguments().size()!=args.size() || std::find(contract.outputRoles.begin(),contract.outputRoles.end(),role)!=contract.outputRoles.end())
            throw std::runtime_error("SST output group mixes independent viscosity with production or repeats mathematics: "+name);
        for(std::size_t i=0;i<args.size();++i)if(!sameFormulaExpression(args[i],equation.rhs.arguments()[i]))throw std::runtime_error("Shared SST operands differ: "+name);
        const auto id=equation.lhs.name();const auto& output=system.state.at(id);
        Math::requireCompatible(output.signature,outputSignatures[role],"SST output "+id);
        if(output.readOnly() || output.components!=1 || output.componentOffset || output.location!=VariableLocation::EulerianCell
            || (output.storageBinding!=StorageBinding::TransientWorkspace && output.storageBinding!=StorageBinding::ProviderDistributed && output.storageBinding!=StorageBinding::NamedDistributed))
            throw std::runtime_error("SST output needs writable scalar backing: "+id);
        for(const auto& previous:contract.outputs)if(system.state.at(previous).storageKey==output.storageKey)throw std::runtime_error("Aliased SST output: "+id);
        contract.outputs.push_back(id);contract.outputRoles.push_back(role);result.calls.push_back({name,id});result.writes.push_back(id);
        result.operatorBindings.push_back(operations[role]);result.stateEffects.push_back({id,true,false,true});
    }
    if(call.equation!=first.id || call.target.symbol!=contract.outputs.front())throw std::runtime_error("SST call must target its declared first output: "+call.occurrence);
    if(!contract.eddyViscosity) {
        contract.viscosity=args[10].name();Math::requireCompatible(args[10].signature(),outputSignatures[1],"published SST dynamic viscosity");
    }
    const auto& velocity=args[1];
    std::string velocityInput;
    if(velocity.kind()==FormulaExpr::Kind::Divide && velocity.arguments()[0].kind()==FormulaExpr::Kind::Symbol
        && sameFormulaExpression(velocity.arguments()[1],args[0])) {
        velocityInput=velocity.arguments()[0].name();contract.packedVelocity=true;
    } else if(velocity.kind()==FormulaExpr::Kind::Symbol)velocityInput=velocity.name();
    else throw std::runtime_error("SST velocity requires a registered vector or explicit momentum/density relation");
    for(const auto& id:{args[0].name(),velocityInput,args[2].name(),args[3].name()}) {
        auto read=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& use){return use.symbol==id;});
        if(read==call.inputs.end() || (read->kind!=TargetKind::Physical && !(call.stepScoped && read->kind==TargetKind::Working)) || read->halo || !read->snapshot.empty() || read->perIteration
            || (read->version!=StateVersion::Stage && read->version!=StateVersion::Current) || (call.stepScoped && read->version!=StateVersion::Current))
            throw std::runtime_error("SST requires explicit Current/Stage input: "+id);
        result.stateUses.push_back(*read);
    }
    if(call.inputs.size()!=(contract.eddyViscosity?4:5))throw std::runtime_error("Duplicate/unconsumed SST inputs: "+call.occurrence);
    if(!contract.eddyViscosity) {
        const auto read=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& use){return use.symbol==contract.viscosity;});
        if(read==call.inputs.end() || read->version!=(call.stageScoped?StateVersion::Stage:StateVersion::Current) || read->kind!=TargetKind::Physical || read->halo || read->perIteration || !read->snapshot.empty())
            throw std::runtime_error("SST production requires its explicitly published viscosity at this evaluation version: "+call.occurrence);
        result.stateUses.push_back(*read);
    }
    const auto& density=system.state.at(args[0].name());const auto& momentum=system.state.at(velocityInput);
    for(const auto& [name,value]:binding.parameters) {
        if(name!="constantDensity" || value!=1 || !density.constantValue || !std::isfinite(*density.constantValue) || *density.constantValue<=0)
            throw std::runtime_error("Selected closure constantDensity requirement needs a declared positive constant density: "+call.occurrence);
    }
    if(contract.packedVelocity && (result.stateUses[1].kind!=TargetKind::Physical || density.storageBinding!=StorageBinding::PackedDistributed || momentum.storageBinding!=StorageBinding::PackedDistributed
        || density.storageKey!="conservative" || momentum.storageKey!=density.storageKey || density.componentOffset!=0 || momentum.componentOffset!=1))
        throw std::runtime_error("SST momentum/density velocity currently requires the original packed backing; named inputs must declare velocity directly");
    if(density.components!=1 || momentum.components!=3)throw std::runtime_error("SST density/vector shape mismatch");
    std::array<double,7> parameters{};
    const std::vector<int> parameterSlots=contract.eddyViscosity?std::vector<int>{0,2,6}:std::vector<int>{0,1,2,3,4,5};
    for(std::size_t n=0;n<parameterSlots.size();++n) {
        const int i=parameterSlots[n];
        const auto& state=system.state.at(args[n+4].name());
        if(!state.constantValue || !std::isfinite(*state.constantValue) || *state.constantValue<0 || (i && *state.constantValue==0))
            throw std::runtime_error("SST requires finite explicit viscosity/nonzero closure parameter: "+state.id);
        parameters[i]=*state.constantValue;result.reads.push_back(state.id);
    }
    for(int i=2;i<4;++i) {
        const auto& scalar=system.state.at(args[i].name());
        if(scalar.components!=1 || scalar.componentOffset || scalar.producer || scalar.constantValue
            || scalar.location!=VariableLocation::EulerianCell)
            throw std::runtime_error("SST scalar input needs actual unsliced boundary-prepared backing: "+scalar.id);
    }
    contract.parameters={parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6]};
    contract.density=density.id;contract.momentum=momentum.id;contract.targets={args[2].name(),args[3].name()};
    result.equationMethod="SSTCentral2";result.backendProvider="closure.sst-central2";
    result.backendOperation=instanceOperation("closure.evaluate",call.occurrence);result.temporalPublish=instanceOperation("closure.publish",call.occurrence);
    result.coefficientPrepare=result.backendOperation;result.coefficientPublish=result.temporalPublish;
    if(call.stepScoped) {
        result.temporalPublish.clear();result.coefficientPublish.clear();
        result.coefficientVersion=StateVersion::Current;
    }
    result.coefficientProvides={instanceOperation("closure.outputs",call.occurrence)};result.coefficientStates=contract.outputs;
    for(const auto& read:result.stateUses)result.reads.push_back(read.symbol);
    FieldDataRequirement fields;
    fields.target=contract.outputs.front();
    fields.states=contract.outputs;
    for(const auto& read:result.stateUses)fields.states.push_back(read.symbol);
    if(!contract.packedVelocity)fields.boundaries={contract.momentum,contract.targets[0],contract.targets[1]};
    result.physicalResources.push_back({call.occurrence+".inputs","native.scalar",typeid(FDM::ScalarEquationData),false,std::move(fields)});
    const bool packedGradient=contract.packedVelocity,eddy=contract.eddyViscosity;
    result.requiredHaloWidth=1;result.providerContract=std::move(contract);
    for(const auto& op:{result.backendOperation,result.temporalPublish})if(!op.empty())result.operations.push_back({op,"evaluate only declared SST outputs",OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,"SSTCentral2"}});
    result.operatorBindings.push_back(packedGradient?"velocity-gradient.Central2.packed":"velocity-gradient.Central2.sampled");if(!eddy)result.operatorBindings.push_back("scalar-gradient.Central2.sampled");
    result.requirements={"serial single patch; stationary non-IBM; actual density/velocity and k/omega Current/Working or Stage input backing", "geometry: original metrics/wallDistance; absent physical wall distance uses declared original far-field 1e20 convention", "explicit WHAT constants; only explicitly called outputs are published; production reads the independent viscosity publication; Stage final Current refresh; no temporal integration"};
    return result;
}

}
