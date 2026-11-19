#include "SF_model.h"
#include "compatibility/SF_compatibility.h"
#include "infrastructure/io/case/SF_case.h"
#include "solver/system/SF_singleFluidPreset.h"
#include "solver/system/SF_presets.h"
#include "models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.h"
#include "models/physics/multiphase/SF_multiphaseSystemContribution.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "models/ibm/SF_ibmSystemContribution.h"
#include "models/ibm/descriptor/SF_algorithmDescriptor.h"
#include <algorithm>
#include <filesystem>
#include <cctype>
namespace SF::Application::ModelLoader {
namespace {
using P=Model::Parameters;
std::string turbulenceName(FDM::TurbulenceModelKind model) {
    switch (model) {
        case FDM::TurbulenceModelKind::None: return "none";
        case FDM::TurbulenceModelKind::kEpsilon: return "kEpsilon";
        case FDM::TurbulenceModelKind::kOmegaSST: return "kOmegaSST";
        case FDM::TurbulenceModelKind::Smagorinsky: return "Smagorinsky";
        case FDM::TurbulenceModelKind::DNS: return "DNS";
    }
    throw std::runtime_error("Unknown turbulence model.");
}

Model::Description expandPresets(const Model::Description& input) {
    auto result=input;
    std::map<std::string,Model::ObjectDescriptor> defaults;
    P numerical=P::object();
    std::map<std::string,std::string> origins;
    const auto role=[](std::string name) {
        for (auto& c:name)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if(name=="builtinequationsystem")name="equationsystem";
        return name;
    };
    const auto explicitlySet=[&](const std::string& path) {
        const std::string numericalPrefix="solver/numerics";
        if(path.rfind(numericalPrefix+"/",0)==0)
            return input.numerics.contains(P::json_pointer(path.substr(numericalPrefix.size())));
        for(const auto& object:input.objects) {
            const auto prefix="object/"+role(object.type);
            if(object.type!="compositionPreset" && path.rfind(prefix+"/",0)==0
                && object.parameters.contains(P::json_pointer(path.substr(prefix.size()))))return true;
        }
        return false;
    };
    const auto record=[&](const std::string& path,const std::string& source,bool explicitValue) {
        const auto old=origins.find(path);
        const auto previous=old==origins.end()?std::string():old->second;
        origins[path]=(explicitValue?"explicit ":"default ")+source;
        if(explicitValue && !previous.empty())origins[path]+="; overrides "+previous;
        else if(!explicitValue && !previous.empty() && previous!=origins[path])origins[path]=previous+"; also default "+source;
    };
    const auto merge=[&](const auto& self,P& target,const P& value,const std::string& path,
                         const std::string& source,bool explicitValue)->void {
        if (value.is_object()) {
            if (!target.is_object()) {
                if (!target.is_null() && !target.empty()) throw std::runtime_error("Conflicting setting shape at "+path+": "+origins[path]+" and "+source);
                target=P::object();
            }
            if(value.empty())record(path,source,explicitValue);
            for (const auto& item:value.items()) self(self,target[item.key()],item.value(),path+"/"+item.key(),source,explicitValue);
            return;
        }
        if (!target.is_null() && target!=value && !explicitValue && !explicitlySet(path))
            throw std::runtime_error("Ambiguous preset default at "+path+": "+origins[path]+" and "+source);
        target=value;record(path,source,explicitValue);
    };
    std::vector<Model::ObjectDescriptor> explicitObjects;
    for (const auto& object:input.objects) {
        if (object.type!="compositionPreset") {explicitObjects.push_back(object);continue;}
        Model::Schema{{{"preset","string",true}},false}.validate(object.parameters,object.source);
        const auto name=object.parameters.at("preset").get<std::string>();
        if(name!="singleFluidLaminar" && name!="singleFluidSST")
            throw std::runtime_error(object.source+": unknown built-in composition preset "+name);
        const bool sst=name=="singleFluidSST";
        const auto source=object.source+" [builtin "+name+"]";
        merge(merge,numerical,P{{"time",{{"default",sst?"classicalRK4":"forwardEuler"}}},
            {"terms",{{"convection","weno5Rusanov"},{"diffusion","central2Explicit"}}}},"solver/numerics",source,false);
        const auto add=[&](const char* id,const char* category,const char* type,P parameters) {
            auto& value=defaults[role(type)];
            if(value.type.empty())value={id,category,type};
            merge(merge,value.parameters,parameters,"object/"+role(type),source,false);
            value.source=source;
        };
        const auto base=System::Preset::singleFluidDefaults();
        add("equations","solver","equationSystem",{{"preset","conservativeExplicit"}});
        add("state","solver","stateRegistry",{{"use",base.solutionVariables}});
        add("equationRegistry","solver","equationRegistry",{{"use",base.equations}});
        add("execution","solver","algorithmRegistry",{{base.algorithm,P::object()}});
        add("turbulence","models","turbulence",sst
            ? P{{"simulationType","RAS"},{"RAS",{{"model","kOmegaSST"}}},{"transportMethod","stageConservative"}}
            : P{{"simulationType","laminar"}});

    }
    const auto numericsSource=input.selectionOrigins.find("solver/numerics");
    merge(merge,numerical,input.numerics,"solver/numerics",numericsSource==input.selectionOrigins.end()?"registered numerics":numericsSource->second,true);
    result.numerics=std::move(numerical);
    result.objects.clear();
    std::map<std::string,std::string> explicitRoles;
    for(auto object:explicitObjects) {
        const auto key=role(object.type);
        if(defaults.count(key)) {
            const auto [prior,inserted]=explicitRoles.emplace(key,object.source);
            if(!inserted)throw std::runtime_error("Conflicting explicit target "+key+": "+prior->second+" and "+object.source);
            auto parameters=defaults.at(key).parameters;
            // Explicit algorithm/model choices replace alternative branches, not merge two algorithms.
            if(key=="algorithmregistry" || key=="equationsystem" || (key=="turbulence" && object.parameters.contains("simulationType"))) {
                for(const auto& item:parameters.items())if(!object.parameters.contains(item.key())) {
                    const auto prefix="object/"+key+"/"+item.key();
                    for(auto it=origins.begin();it!=origins.end();) {
                        if(it->first==prefix || it->first.rfind(prefix+"/",0)==0)it=origins.erase(it);else ++it;
                    }
                }
                parameters=P::object();
            }
            merge(merge,parameters,object.parameters,"object/"+key,object.source,true);
            object.parameters=std::move(parameters);
        } else {
            P empty;merge(merge,empty,object.parameters,"object/"+key,object.source,true);
        }
        result.objects.push_back(std::move(object));
    }
    for(auto& [key,object]:defaults)if(!explicitRoles.count(key))result.objects.push_back(std::move(object));
    origins["solver/numerics"]=numericsSource==input.selectionOrigins.end()?"registered numerics":numericsSource->second;
    if(!result.numerics.contains("time"))origins["solver/numerics/time/default"]="builtin NumericsConfig default; no explicit time selection";
    result.selectionOrigins=std::move(origins);
    return result;
}
} // namespace
std::vector<System::SystemContribution> modelDeclarations(const CaseConfig& config) {
    std::vector<System::SystemContribution> contributions;
    const bool active=config.multiPhaseEnabled && config.multiPhase.enabled;
    const bool levelSet=active && Physics::Multiphase::isLevelSetType(config.multiPhase.type);
    const bool homogeneous=active && Physics::Multiphase::isHomogeneousType(config.multiPhase.type);
    const bool mixture=active && Physics::Multiphase::isMixtureType(config.multiPhase.type);
    const bool eulerian=active && Physics::Multiphase::isEulerianEulerianType(config.multiPhase.type);
    if (levelSet) {
        System::SystemContribution contribution;
        Physics::InterfaceModels::LevelSetContribution::contribute(
            contribution,{Physics::Multiphase::normalizeModelType(
                config.multiPhase.levelSet.surfaceTensionModel)=="ghostfluid",
                config.multiPhase.levelSet.advectionOrder,
                config.multiPhase.levelSet.reinitializationOrder,
                config.multiPhase.levelSet.reinitializationSteps,
                config.multiPhase.levelSet.pseudoTimeStep,
                config.multiPhase.levelSet.wenoEpsilon,
                config.multiPhase.levelSet.wenoPower,
                config.multiPhase.levelSet.signSmoothingFactor,
                Physics::Multiphase::normalizeModelType(config.multiPhase.levelSet.surfaceTensionModel)=="csf",
                config.multiPhase.levelSet.surfaceTension,config.multiPhase.levelSet.interfaceThickness});
        contributions.push_back(std::move(contribution));
    }
    if (config.multiPhaseEnabled && config.multiPhase.enabled && Physics::Multiphase::isThermalType(config.multiPhase.type)) {
        System::SystemContribution contribution;
        Physics::Multiphase::contributeThermal(contribution,config.multiPhase,
            config.solver.boundaries.ilwEnabled,config.solver.boundaries.ilwOrder);
        contributions.push_back(std::move(contribution));
    }
    if (homogeneous || mixture) {
        System::SystemContribution contribution;
        if (homogeneous)
            Physics::Multiphase::contributeHomogeneous(contribution,config.multiPhase);
        else Physics::Multiphase::contributeMixture(contribution,config.multiPhase);
        contributions.push_back(std::move(contribution));
    }
    if (config.solver.turbulence.enabled) {
        if(config.solver.turbulence.fixedVolumeTransport) {
            if(!config.solver.thermophysical || config.solver.thermophysical->selection.equationOfState!="rhoConst")
                throw std::runtime_error("fixedVolume SST requires an explicitly declared rhoConst closure; variable capacity/mass transport is not implemented");
            if(config.solver.turbulence.coupleEnergy)
                throw std::runtime_error("fixedVolume SST has no energy transport feedback; explicitly disable coupleEnergy");
            if(!config.solver.turbulence.closureDensity.empty() && config.solver.turbulence.closureDensity!="rho")
                throw std::runtime_error("fixedVolume SST currently consumes the declared rhoConst symbol rho; alternate capacity mapping is not implemented");
        }
        System::SystemContribution contribution;
        const auto& turbulence=config.solver.turbulence;
        Turbulence::SystemContributionSpec spec{
            turbulenceName(turbulence.model),turbulence.phaseNames,eulerian,
            config.multiPhase.eulerianEulerian.phaseNames,turbulence.stageConservativeTransport,
            turbulence.laminarDynamicViscosity,turbulence.coefficients.omegaFloor,
            turbulence.closureDensity,turbulence.closureVelocity,turbulence.fixedVolumeTransport,turbulence.sharedVolumeFlux};
        spec.viscosityOutput=turbulence.viscosityOutput;
        spec.viscosityEquation=turbulence.viscosityEquation;
        const bool explicitStress=!turbulence.stressEquation.empty() || !turbulence.stressOccurrence.empty();
        if(turbulence.fixedVolumeTransport && turbulence.coupleMomentum) {
            if(explicitStress) {
                if(turbulence.stressEquation.empty() || turbulence.stressOccurrence.empty())
                    throw std::runtime_error("momentumStress requires both equation and occurrence");
                spec.stressTarget=std::pair{turbulence.stressEquation,turbulence.stressOccurrence};
            } else {
                if(!config.composition.predictorEquation.empty())
                    throw std::runtime_error("Explicit pressure predictor="+config.composition.predictorEquation+
                        " requires turbulence/momentumStress with its actual equation and occurrence");
                spec.stressTarget=std::pair{System::Compose::pressureMomentumEquation,System::Compose::pressureMomentumDiffusion};
            }
        } else if(explicitStress) {
            throw std::runtime_error("momentumStress requires active fixedVolume momentum coupling");
        }
        Turbulence::contribute(contribution,spec);
        contributions.push_back(std::move(contribution));
    }
    if (config.solver.ibm.enabled) {
        auto immersedDescriptor = config.solver.ibm.method == FDM::IBMMethod::Ghost
            ? IBM::Descriptor::ghostCell()
            : IBM::Descriptor::variational(config.solver.ibm.forcing);
        immersedDescriptor.wallClosure=config.solver.ibm.wallClosure;
        System::SystemContribution contribution;
        IBM::SystemContribution::contribute(contribution,immersedDescriptor);
        contributions.push_back(std::move(contribution));
    }
    return contributions;
}

std::vector<EquationSystemInstanceConfig> equationInstances(
        const Model::Description& model) {
    std::vector<EquationSystemInstanceConfig> result;
    for (const auto& object : model.objects) {
        if (object.type != "equationSystem"
            && object.type != "BuiltinEquationSystem") continue;
        if (object.name.empty() || !object.parameters.contains("preset")) {
            throw std::runtime_error(
                "An equation-system instance requires a name and a built-in preset.");
        }
        const std::string preset=object.parameters.at("preset").get<std::string>();
        // Instance metadata records provenance; execution remains provider-driven.
        const std::string type=preset=="compressible" || preset=="conservativeExplicit"
            ? "singleFluid" : preset=="eulerianEulerian" ? "eulerian" : std::string{};
        if (type.empty()) {
            throw std::runtime_error(
                "Equation-system instance '"+object.name
                +"' has unsupported type '"+preset+"'.");
        }
        for (const auto& existing : result) {
            if (existing.name == object.name) {
                throw std::runtime_error(
                    "Duplicate equation-system instance name '"
                    +object.name+"'.");
            }
        }
        result.push_back({object.name,type});
    }
    return result;
}

Model::FactoryRegistry<CaseConfig>& objectRegistry() {
    static Model::FactoryRegistry<CaseConfig> registry;
    return registry;
}
CaseConfig build(const Model::Description& model) {
    auto builtins=expandPresets(model);
    const auto instances = equationInstances(builtins);
    const auto objects=builtins.objects;
    std::vector<Model::ObjectDescriptor> extensions;
    builtins.objects.clear();
    for(const auto& object:objects) {
        if(objectRegistry().contains(object.type))extensions.push_back(object);
        else builtins.objects.push_back(object);
    }
    auto result=CaseAdapter(model.caseDir).build(builtins);
    for(const auto& object:extensions)objectRegistry().create(object,result);
    auto authored=model;authored.selectionOrigins=builtins.selectionOrigins;
    result.modelDescription=std::make_shared<const Model::Description>(std::move(authored));
    result.composition.instances=instances;
    return result;
}
CaseConfig read(const std::string& path) {
    const std::filesystem::path p(path);
    const auto dir=std::filesystem::is_directory(p)?p:p.parent_path();
    if(std::filesystem::exists(dir/"case.yaml"))return build(CaseIO::read(path));
    throw std::runtime_error(
        "Legacy OpenFOAM case loading was removed; convert the case to case.yaml "
        "(see SF_nativeDecode.cpp / decodeNativeCase) before loading: "+path);
}
}
