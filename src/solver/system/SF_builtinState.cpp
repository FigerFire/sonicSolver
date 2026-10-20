#include "SF_builtinState.h"

namespace SF::System {
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
        symbol.derivation=StateDerivation::Velocity;break;
    case BuiltinState::Pressure:
        symbol.id="p";symbol.name="pressure";symbol.derivation=StateDerivation::Pressure;break;
    case BuiltinState::Temperature:
        symbol.id="T";symbol.name="temperature";symbol.derivation=StateDerivation::Temperature;break;
    case BuiltinState::Enthalpy:
        symbol.id="h";symbol.name="specific enthalpy";symbol.derivation=StateDerivation::Enthalpy;break;
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
        symbol.id="mu";symbol.name="dynamic viscosity";symbol.derivation=StateDerivation::DynamicViscosity;break;
    case BuiltinState::KinematicViscosity:
        symbol.id="nu";symbol.name="kinematic viscosity";symbol.derivation=StateDerivation::KinematicViscosity;break;
    case BuiltinState::ThermalConductivity:
        symbol.id="conductivity";symbol.name="thermal conductivity";symbol.derivation=StateDerivation::ThermalConductivity;break;
    }
    if (symbol.derivation==StateDerivation::Velocity) {
        symbol.evaluation=StateEvaluation::Lazy;symbol.dependencies={"rho","rhoU"};
    }
    if (symbol.derivation==StateDerivation::Pressure || symbol.derivation==StateDerivation::Temperature
        || symbol.derivation==StateDerivation::Enthalpy) {
        symbol.evaluation=StateEvaluation::Lazy;symbol.dependencies={"rho","rhoU","rhoE"};
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
        if (symbol.evaluation==StateEvaluation::Lazy && !active.contains("rho")) {
            std::vector<std::string> densityComponents;
            for (const auto& item:active.symbols()) if (item.id.rfind("partialDensity.",0)==0) densityComponents.push_back(item.id);
            if (!densityComponents.empty()) {
                symbol.dependencies.erase(std::remove(symbol.dependencies.begin(),symbol.dependencies.end(),"rho"),symbol.dependencies.end());
                symbol.dependencies.insert(symbol.dependencies.end(),densityComponents.begin(),densityComponents.end());
            }
        }
        active.add(std::move(symbol));
    }
}

StateSymbol BuiltinStateCatalog::solution(std::string_view id) const {
    auto symbol=at(id);
    symbol.derivation=StateDerivation::None;
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
