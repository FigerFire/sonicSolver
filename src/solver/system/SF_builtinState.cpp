#include "SF_builtinState.h"
#include "core/state/SF_distributedField.h"

namespace SF::System {
namespace {
// All projections borrow the Field-owned EOS cache. Binding never closes EOS;
// the original Q invalidation/version contract governs cache hits and misses.
std::shared_ptr<const StateProducer> cacheProjection(const StateSymbol& symbol,
        Math::Signature (*result)(const std::vector<Math::Signature>&),
        std::function<double(const Physics::FluidStateModel::ThermodynamicState&)> project) {
    auto producer=std::make_shared<StateProducer>();
    producer->id="closure.thermodynamic."+symbol.id;
    std::vector<FormulaExpr> arguments;
    for (const auto& dependency:symbol.dependencies) {
        if (dependency=="rho") arguments.push_back(builtinState(BuiltinState::Density).expression());
        else if (dependency=="rhoU") arguments.push_back(builtinState(BuiltinState::Momentum).expression());
        else if (dependency=="rhoE") arguments.push_back(builtinState(BuiltinState::TotalEnergy).expression());
        else throw std::runtime_error("Unsupported thermodynamic cache operand: "+dependency);
    }
    const auto operation="Field.thermodynamicState."+symbol.id;
    const Math::OperatorRule rule{operation.c_str(),result};
    producer->expression=FormulaExpr::op(rule,std::move(arguments));
    producer->versions={StateViewKind::Physical,StateViewKind::Stage};
    producer->bind=[project=std::move(project)](Field& field,const auto&,State::DistributedFieldView& view)->std::function<void()> {
        view.read=[&field,project](int cell,int) {
            int i=0,j=0,k=0;field.getIJK(cell,i,j,k);
            return project(field.thermodynamicState(i,j,k));
        };
        return {};
    };
    return producer;
}
void requireConservativeCacheArguments(const std::vector<Math::Signature>& args) {
    if(args.size()!=3)throw std::runtime_error("Thermodynamic projection requires density, momentum and energy");
    const BuiltinState variables[]={BuiltinState::Density,BuiltinState::Momentum,BuiltinState::TotalEnergy};
    for(unsigned i=0;i<3;++i)Math::requireCompatible(args[i],builtinState(variables[i]).signature,"thermodynamic cache input");
}
void requireOwnerCacheArguments(const std::vector<Math::Signature>& args) {
    if(!args.empty())throw std::runtime_error("Field cache projection has no independent STATE operands");
}
}

StateSymbol builtinState(BuiltinState variable) {
    StateSymbol symbol;
    symbol.origin={OriginKind::BuiltinPreset,"builtin state catalog"};
    symbol.nameSpace="fluid";
    symbol.initializationRequired=false;
    symbol.boundaryRequired=false;
    symbol.runtimeStorageRequired=false;
    symbol.role=StateRole::Derived;
    switch (variable) {
    case BuiltinState::Density:
        symbol.id="rho";symbol.name="density";symbol.storageKey="conservative";break;
    case BuiltinState::Momentum:
        symbol.id="rhoU";symbol.name="momentum density";
        symbol.components=3;symbol.componentOffset=1;symbol.storageKey="conservative";break;
    case BuiltinState::TotalEnergy:
        symbol.id="rhoE";symbol.name="total energy density";
        symbol.componentOffset=4;symbol.storageKey="conservative";break;
    case BuiltinState::Velocity:
        symbol.id="U";symbol.name="velocity";symbol.components=3;
        break;
    case BuiltinState::Pressure:
        symbol.id="p";symbol.name="pressure";break;
    case BuiltinState::Temperature:
        symbol.id="T";symbol.name="temperature";break;
    case BuiltinState::Enthalpy:
        symbol.id="h";symbol.name="specific enthalpy";break;
    case BuiltinState::FaceFlux:
        symbol.id="phi";symbol.name="canonical volume face flux";
        symbol.location=VariableLocation::EulerianFace;
        symbol.ownership=OwnershipKind::CanonicalFace;symbol.storageKey="pressureFaceFlux";break;
    case BuiltinState::VolumeFraction:
        symbol.id="alpha";symbol.name="volume fraction";symbol.storageKey="alpha";break;
    case BuiltinState::TurbulentEnergy:
        symbol.id="k";symbol.name="turbulent kinetic energy";symbol.storageKey="k";break;
    case BuiltinState::SpecificDissipation:
        symbol.id="omega";symbol.name="specific dissipation";symbol.storageKey="omega";break;
    case BuiltinState::DynamicViscosity:
        symbol.id="mu";symbol.name="dynamic viscosity";break;
    case BuiltinState::KinematicViscosity:
        symbol.id="nu";symbol.name="kinematic viscosity";break;
    case BuiltinState::ThermalConductivity:
        symbol.id="conductivity";symbol.name="thermal conductivity";break;
    }
    switch(variable) {
    case BuiltinState::Density:symbol.signature=Math::signature(0,{1,-3,0,0,0,0,0});break;
    case BuiltinState::Momentum:symbol.signature=Math::signature(1,{1,-2,-1,0,0,0,0});break;
    case BuiltinState::TotalEnergy:case BuiltinState::Pressure:symbol.signature=Math::signature(0,{1,-1,-2,0,0,0,0});break;
    case BuiltinState::Velocity:symbol.signature=Math::signature(1,{0,1,-1,0,0,0,0});break;
    case BuiltinState::Temperature:symbol.signature=Math::signature(0,{0,0,0,1,0,0,0});break;
    case BuiltinState::Enthalpy:case BuiltinState::TurbulentEnergy:symbol.signature=Math::signature(0,{0,2,-2,0,0,0,0});break;
    case BuiltinState::SpecificDissipation:symbol.signature=Math::signature(0,{0,0,-1,0,0,0,0});break;
    case BuiltinState::DynamicViscosity:symbol.signature=Math::signature(0,{1,-1,-1,0,0,0,0});break;
    case BuiltinState::KinematicViscosity:symbol.signature=Math::signature(0,{0,2,-1,0,0,0,0});break;
    case BuiltinState::ThermalConductivity:symbol.signature=Math::signature(0,{1,1,-3,-1,0,0,0});break;
    case BuiltinState::FaceFlux:symbol.signature=Math::signature(0,{0,3,-1,0,0,0,0});break;
    case BuiltinState::VolumeFraction:symbol.signature=Math::signature(0,{});break;
    }
    if (variable==BuiltinState::Velocity) {
        symbol.evaluation=StateEvaluation::Lazy;symbol.dependencies={"rho","rhoU"};
        auto producer=std::make_shared<StateProducer>();producer->id=velocityProducer;
        producer->expression=FormulaExpr::divide(FormulaExpr::symbol("rhoU",Math::signature(1,{1,-2,-1,0,0,0,0})),
            FormulaExpr::symbol("rho",Math::signature(0,{1,-3,0,0,0,0,0})));
        producer->bind=[](Field& field,const auto&,State::DistributedFieldView& view)->std::function<void()> {
            // Preserve the original EOS cache, including its validation and miss timing.
            view.read=[&field](int cell,int component) {
                int i=0,j=0,k=0;field.getIJK(cell,i,j,k);
                return field.thermodynamicState(i,j,k).velocity.at((std::size_t)component);
            };
            return {};
        };
        symbol.producer=std::move(producer);
    }
    if (variable==BuiltinState::Pressure || variable==BuiltinState::Temperature || variable==BuiltinState::Enthalpy) {
        symbol.evaluation=StateEvaluation::Lazy;symbol.dependencies={"rho","rhoU","rhoE"};
    }
    // Transport properties also borrow this cache, including for primitive-state
    // providers that have no independent conservative STATE registration.
    // Their zero-argument projection names the bound Field owner, not fake Q symbols.
    switch(variable) {
    case BuiltinState::Pressure:
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireConservativeCacheArguments(a);return Math::signature(0,{1,-1,-2,0,0,0,0});},
            [](const auto& v){return v.pressure;});break;
    case BuiltinState::Temperature:
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireConservativeCacheArguments(a);return Math::signature(0,{0,0,0,1,0,0,0});},
            [](const auto& v){return v.temperature;});break;
    case BuiltinState::Enthalpy:
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireConservativeCacheArguments(a);return Math::signature(0,{0,2,-2,0,0,0,0});},
            [](const auto& v){return (v.internalEnergyDensity+v.pressure)/v.density;});break;
    case BuiltinState::DynamicViscosity:
        symbol.evaluation=StateEvaluation::Lazy;
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireOwnerCacheArguments(a);return Math::signature(0,{1,-1,-1,0,0,0,0});},
            [](const auto& v){return v.dynamicViscosity;});break;
    case BuiltinState::KinematicViscosity:
        symbol.evaluation=StateEvaluation::Lazy;
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireOwnerCacheArguments(a);return Math::signature(0,{0,2,-1,0,0,0,0});},
            [](const auto& v){return v.dynamicViscosity/v.density;});break;
    case BuiltinState::ThermalConductivity:
        symbol.evaluation=StateEvaluation::Lazy;
        symbol.producer=cacheProjection(symbol,[](const auto& a){requireOwnerCacheArguments(a);return Math::signature(0,{1,1,-3,-1,0,0,0});},
            [](const auto& v){return v.thermalConductivity;});break;
    default:break;
    }
    symbol.shape=symbol.components==1 ? ValueShape::Scalar : ValueShape::Vector;
    if (!symbol.storageKey.empty()) symbol.storageBinding=StorageBinding::NamedDistributed;
    return symbol;
}

BuiltinStateCatalog::BuiltinStateCatalog() {
    for (const auto variable:{BuiltinState::Density,BuiltinState::Momentum,
            BuiltinState::TotalEnergy,BuiltinState::Velocity,BuiltinState::Pressure,
            BuiltinState::Temperature,BuiltinState::Enthalpy,BuiltinState::FaceFlux,
            BuiltinState::VolumeFraction,BuiltinState::TurbulentEnergy,
            BuiltinState::SpecificDissipation,BuiltinState::DynamicViscosity,
            BuiltinState::KinematicViscosity,BuiltinState::ThermalConductivity}) {
        symbols_.push_back(builtinState(variable));
    }
    for(const auto& item:std::vector<std::pair<std::string,Math::Signature>>{
            {"massFlux",Math::signature(1,{1,-2,-1,0,0,0,0})},
            {"momentumFlux",Math::signature(2,{1,-1,-2,0,0,0,0})},
            {"energyFlux",Math::signature(1,{1,0,-3,0,0,0,0})}}) {
        StateSymbol symbol;symbol.id=symbol.name=item.first;symbol.signature=item.second;
        symbol.components=item.second.rank==2?9:3;symbol.shape=item.second.rank==2?ValueShape::Tensor:ValueShape::Vector;
        symbol.role=StateRole::Derived;symbol.runtimeStorageRequired=false;
        symbol.initializationRequired=symbol.boundaryRequired=symbol.outputEligible=symbol.restartEligible=false;
        symbols_.push_back(std::move(symbol));
    }
}

bool BuiltinStateCatalog::contains(std::string_view id) const {
    return std::any_of(symbols_.begin(),symbols_.end(),
        [&](const auto& symbol) { return symbol.id==id; });
}

const StateSymbol& BuiltinStateCatalog::at(std::string_view id) const {
    const auto found=std::find_if(symbols_.begin(),symbols_.end(),
        [&](const auto& symbol) { return symbol.id==id; });
    if (found==symbols_.end())
        throw std::runtime_error("Unknown builtin STATE symbol '"+std::string(id)
            +"'; declare custom metadata with addState.");
    return *found;
}

void BuiltinStateCatalog::require(StateRegistry& active,std::string_view id) const {
    if (!active.contains(id)) {
        auto symbol=at(id);
        if (std::find(symbol.dependencies.begin(),symbol.dependencies.end(),"rho")!=symbol.dependencies.end()
            && !active.contains("rho")) {
            std::vector<std::string> densityComponents;
            for (const auto& item:active.symbols()) if (item.id.rfind("partialDensity.",0)==0) densityComponents.push_back(item.id);
            if (!densityComponents.empty()) {
                symbol.dependencies.erase(std::remove(symbol.dependencies.begin(),symbol.dependencies.end(),"rho"),symbol.dependencies.end());
                symbol.dependencies.insert(symbol.dependencies.end(),densityComponents.begin(),densityComponents.end());
                if (symbol.producer) {
                    auto density=active.expression(densityComponents.front());
                    for (std::size_t i=1;i<densityComponents.size();++i) density=density+active.expression(densityComponents[i]);
                    auto producer=std::make_shared<StateProducer>(*symbol.producer);
                    const auto substitute=[&](const auto& self,const FormulaExpr& expression)->FormulaExpr {
                        if (expression.kind()==FormulaExpr::Kind::Symbol && expression.name()=="rho") return density;
                        if (expression.arguments().empty()) return expression;
                        std::vector<FormulaExpr> args;
                        for (const auto& a:expression.arguments()) args.push_back(self(self,a));
                        return expression.withArguments(std::move(args));
                    };
                    producer->expression=substitute(substitute,producer->expression);symbol.producer=std::move(producer);
                }
            }
        }
        active.add(std::move(symbol));
    }
}

StateSymbol BuiltinStateCatalog::solution(std::string_view id) const {
    auto symbol=at(id);
    if(symbol.storageKey.empty() && !symbol.readOnly()
        && !symbol.outputEligible && !symbol.restartEligible)
        throw std::runtime_error("Logical mathematical quantity has no writable STATE storage: "+std::string(id));
    symbol.producer.reset();
    symbol.evaluation=StateEvaluation::Direct;symbol.dependencies.clear();
    symbol.initializationRequired=true;
    symbol.boundaryRequired=true;
    symbol.runtimeStorageRequired=true;
    symbol.role=StateRole::Primary;
    if (id=="rho" || id=="rhoU" || id=="rhoE")
        symbol.storageBinding=StorageBinding::PackedDistributed;
    else if (id=="U") {
        symbol.storageBinding=StorageBinding::PackedDistributed;
        symbol.storageKey="velocity";
    } else if (id=="p") {
        symbol.role=StateRole::Algebraic;
        symbol.storageBinding=StorageBinding::NamedDistributed;
        symbol.storageKey="pressure";
    } else {
        symbol.storageBinding=StorageBinding::NamedDistributed;
        symbol.storageKey=std::string(id);
    }
    return symbol;
}

void requireTargetStates(StateRegistry& active,const ExecutionScope& scope) {
    const BuiltinStateCatalog catalog;
    const auto request=[&](const auto& self,const ExecutionScope& node)->void {
        if (node.kind==ExecutionKind::EquationCall && node.step.target.kind!=TargetKind::Workspace)
            catalog.require(active,node.step.target.symbol);
        for (const auto& child:node.children) self(self,child);
    };
    request(request,scope);
}
} // namespace SF::System
