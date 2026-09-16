/// @file SF_systemBuilder.cpp
/// @brief 显式声明单流体、多相、湍流和 IBM 对数学系统的贡献。

#include "SF_systemBuilder.h"
#include "SF_equationContribution.h"

#include "SF_config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::System {
namespace {

std::string timeName(FDM::TimeScheme scheme) {
    switch (scheme) {
        case FDM::TimeScheme::Euler: return "Euler";
        case FDM::TimeScheme::SSPRK3: return "SSP-RK3";
        case FDM::TimeScheme::RK4: return "RK4";
    }
    throw std::runtime_error("Unknown time integration scheme.");
}

void addUnknown(
        ResolvedSimulationSystem& system,
        std::string id,
        std::string name,
        int components = 1,
        UnknownRole role = UnknownRole::Primary,
        StorageBinding binding = StorageBinding::SpecializedExecutor,
        std::string storageKey = {},
        int componentOffset = 0,
        std::string nameSpace = {}) {
    UnknownDescriptor unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = role;
    unknown.storageBinding = binding;
    unknown.storageKey = std::move(storageKey);
    unknown.componentOffset = componentOffset;
    unknown.nameSpace = std::move(nameSpace);
    EquationSystemBuilder(system).addUnknown(std::move(unknown));
}

void addEquation(
        ResolvedSimulationSystem& system,
        EquationDescriptor descriptor,
        Equation::Definition definition) {
    EquationSystemBuilder(system).addEquation(
        std::move(descriptor),std::move(definition));
}

void addDensityBasedFluid(ResolvedSimulationSystem& system) {
    addUnknown(system,"rho","density",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",0,"fluid");
    addUnknown(system,"rhoU","momentum",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"rhoE","total energy",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",4,"fluid");
    addEquation(system,{"E_MASS","mass","conservation",{"rho"}},
        Equation::named("E_MASS",
            Equation::ddt({"rho"}) + Equation::div({"massFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_MOMENTUM","momentum","conservation",{"rhoU"}},
        Equation::named("E_MOMENTUM",
            Equation::ddt({"rhoU"}) + Equation::div({"momentumFlux"})
                + Equation::diffusion({"mu"},{"U"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_ENERGY","total energy","conservation",{"rhoE"}},
        Equation::named("E_ENERGY",
            Equation::ddt({"rhoE"}) + Equation::div({"energyFlux"})
                + Equation::diffusion({"conductivity"},{"T"})
                == Equation::Symbol{"zero"}));
    system.solveBlocks.push_back({
        "S_FLUID","compressible conservative predictor",
        system.timeIntegrator,
        {"E_MASS","E_MOMENTUM","E_ENERGY"},{},
        {"rho","rhoU","rhoE"},
        FDM::SolveStrategyKind::ExplicitTimeIntegration});
}

void addPressureBasedFluid(
        ResolvedSimulationSystem& system,
        const std::string& pressureStrategy) {
    addUnknown(system,"rho","density",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",0,"fluid");
    addUnknown(system,"U","velocity",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"rhoE","total energy",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",4,"fluid");
    addUnknown(system,"pPrime","pressure correction",1,
               UnknownRole::Algebraic,StorageBinding::SpecializedExecutor,
               {},0,"pressure");
    addEquation(system,{"E_MASS","mass","conservation",{"rho"}},
        Equation::named("E_MASS",
            Equation::ddt({"rho"}) + Equation::div({"massFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_MOMENTUM","momentum predictor","conservation",{"U"}},
        Equation::named("E_MOMENTUM",
            Equation::ddt({"U"}) + Equation::div({"momentumFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_ENERGY","total energy","conservation",{"rhoE"}},
        Equation::named("E_ENERGY",
            Equation::ddt({"rhoE"}) + Equation::div({"energyFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_PRESSURE","pressure correction","constraint",{"pPrime","U"}},
        Equation::named("E_PRESSURE",
            Equation::constraint({"pressureVelocityConsistency"})
                == Equation::Symbol{"zero"}));
    system.solveBlocks.push_back({
        "S_PREDICTOR","compressible predictor",system.timeIntegrator,
        {"E_MASS","E_MOMENTUM","E_ENERGY"},{},
        {"rho","U","rhoE"},FDM::SolveStrategyKind::SegregatedPredictor});
    system.solveBlocks.push_back({
        "S_PRESSURE","pressure-velocity correction",pressureStrategy,
        {"E_PRESSURE"},{},{"pPrime","U"},
        FDM::SolveStrategyKind::PressureCorrection});
}

void addPhaseEquationPack(
        ResolvedSimulationSystem& system,
        const std::string& phase,
        std::size_t phaseIndex) {
    const std::string suffix = "."+phase;
    const std::string storage = "phase"+std::to_string(phaseIndex)+".";
    addUnknown(system,"phaseMass"+suffix,"phase mass "+phase,1,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"mass",0,phase);
    addUnknown(system,"momentum"+suffix,"phase momentum "+phase,3,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"momentum",0,phase);
    addUnknown(system,"enthalpy"+suffix,"phase total enthalpy "+phase,1,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"enthalpy",0,phase);
    addEquation(system,{
        "E_CONTINUITY"+suffix,"phase continuity "+phase,
        "conservation",{"phaseMass"+suffix}},
        Equation::named("E_CONTINUITY"+suffix,
            Equation::ddt({"phaseMass"+suffix})
                + Equation::div({"phaseMassFlux"+suffix})
                == Equation::Symbol{"phaseMassSources"+suffix}));
    addEquation(system,{
        "E_MOMENTUM"+suffix,"phase momentum "+phase,
        "conservation",{"momentum"+suffix}},
        Equation::named("E_MOMENTUM"+suffix,
            Equation::ddt({"momentum"+suffix})
                + Equation::div({"phaseMomentumFlux"+suffix})
                + Equation::gradient({"alphaPressure"+suffix})
                + Equation::diffusion(
                    {"effectivePhaseViscosity"+suffix},{"U"+suffix})
                == Equation::Symbol{"phaseMomentumSources"+suffix}));
    addEquation(system,{
        "E_ENTHALPY"+suffix,"phase enthalpy "+phase,
        "conservation",{"enthalpy"+suffix}},
        Equation::named("E_ENTHALPY"+suffix,
            Equation::ddt({"enthalpy"+suffix})
                + Equation::div({"phaseEnthalpyFlux"+suffix})
                + Equation::diffusion(
                    {"effectivePhaseConductivity"+suffix},{"h"+suffix})
                == Equation::Symbol{"phaseEnthalpySources"+suffix}));
}

void addSharedPressureConstraint(ResolvedSimulationSystem& system) {
    addEquation(system,{
        "E_SHARED_PRESSURE","shared pressure correction","constraint",{"p"}},
        Equation::named("E_SHARED_PRESSURE",
            Equation::constraint({"sharedPressure"})
                == Equation::Symbol{"zero"}));
    system.constraints.push_back({
        "C_SHARED_PRESSURE","shared pressure relationship",
        "p.phase = p",""});
    system.constraints.push_back({
        "C_VOLUME_FRACTION","volume-fraction closure",
        "sum(alpha.phase) = 1",""});
}

void addEulerianEulerianSolveBlock(ResolvedSimulationSystem& system) {
    system.solveBlocks.push_back({
        "S_EE_PIMPLE","Eulerian-Eulerian pressure coupling",
        "PIMPLE phase predictor/corrector",{},
        {"C_SHARED_PRESSURE","C_VOLUME_FRACTION"},{"p"},
        FDM::SolveStrategyKind::PressureVelocityCoupling});
    for (const auto& equation : system.equations) {
        system.solveBlocks.back().equations.push_back(equation.id);
        for (const auto& unknown : equation.solvedUnknowns) {
            if (unknown != "p") {
                system.solveBlocks.back().unknowns.push_back(unknown);
            }
        }
    }
}

void addEulerianEulerianTemplate(
        ResolvedSimulationSystem& system,
        const std::vector<std::string>& names) {
    if (names.size() < 2) {
        throw std::runtime_error(
            "Resolved Eulerian-Eulerian system requires at least two phases.");
    }
    addUnknown(system,"p","shared pressure",1,UnknownRole::Algebraic,
               StorageBinding::NamedDistributed,"pressure",0,"pressure");
    for (std::size_t phase = 0; phase < names.size(); ++phase) {
        addPhaseEquationPack(system,names[phase],phase);
        const std::string prefix = "phase"+std::to_string(phase)+".";
        system.workspaceRequirements.push_back({
            prefix+"momentumDiagonal",1,VariableLocation::EulerianCell,
            OwnershipKind::EulerianGlobalDof});
        system.workspaceRequirements.push_back({
            prefix+"volumeFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
        system.workspaceRequirements.push_back({
            prefix+"massFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
    }
    system.workspaceRequirements.push_back({
        "pressureCorrection",1,VariableLocation::EulerianCell,
        OwnershipKind::EulerianGlobalDof});
    addSharedPressureConstraint(system);
    addEulerianEulerianSolveBlock(system);
}

void attachEquationToPrimarySolveBlock(
        ResolvedSimulationSystem& system,
        const std::string& equation,
        const std::string& unknown) {
    const auto found = std::find_if(
        system.solveBlocks.begin(),system.solveBlocks.end(),
        [](const SolveBlock& block) {
            return block.id == "S_FLUID" || block.id == "S_PREDICTOR";
        });
    if (found == system.solveBlocks.end()) {
        throw std::runtime_error(
            "Auxiliary equation '"+equation
            +"' has no primary fluid solve block.");
    }
    found->equations.push_back(equation);
    found->unknowns.push_back(unknown);
}

void addTurbulenceEquations(
        ResolvedSimulationSystem& system,
        const BuildRequest& request) {
    if (!request.turbulence) return;
    const bool kEpsilon = request.turbulenceModel == "kEpsilon";
    const bool kOmega = request.turbulenceModel == "kOmegaSST";
    if (!kEpsilon && !kOmega) {
        system.closures.push_back(
            request.turbulenceModel.empty()
                ? "turbulence closure"
                : "turbulence: "+request.turbulenceModel);
        return;
    }

    SolveBlock block{
        "S_TURBULENCE","transported turbulence equations",
        system.timeIntegrator,{},{},{},
        FDM::SolveStrategyKind::ExplicitTimeIntegration};
    std::vector<std::string> phases = request.turbulencePhaseNames;
    if (request.templateOrigin != PhysicsTemplateKind::EulerianEulerian) {
        phases = {""};
    } else if (phases.empty()) {
        throw std::runtime_error(
            "Eulerian transported turbulence equations require explicit phases.");
    }
    for (std::size_t phaseIndex = 0; phaseIndex < phases.size(); ++phaseIndex) {
        const std::string& phase = phases[phaseIndex];
        const std::string suffix = phase.empty() ? "" : "."+phase;
        const std::string k = "k"+suffix;
        const std::string second = (kEpsilon ? "epsilon" : "omega")+suffix;
        const std::string kEquation = "E_TURB_K"+suffix;
        const std::string secondEquation =
            (kEpsilon ? "E_TURB_EPSILON" : "E_TURB_OMEGA")+suffix;
        const std::string storagePrefix = phase.empty()
            ? "" : "phase"+std::to_string(phaseIndex)+".";
        addUnknown(system,k,"turbulent kinetic energy"+(phase.empty()
            ? "" : " "+phase),1,UnknownRole::Transported,
            StorageBinding::NamedDistributed,storagePrefix+"k",0,
            phase.empty() ? "turbulence" : phase);
        addUnknown(system,second,
            std::string(kEpsilon ? "turbulence dissipation "
                                 : "specific dissipation ")+phase,
            1,UnknownRole::Transported,StorageBinding::NamedDistributed,
            storagePrefix+(kEpsilon ? "epsilon" : "omega"),0,
            phase.empty() ? "turbulence" : phase);
        addEquation(system,{kEquation,"turbulent kinetic energy transport",
                            "conservation",{k}},
            Equation::named(kEquation,
                Equation::ddt({k}) + Equation::div({"flux."+k})
                    + Equation::diffusion({"diffusivity."+k},{k})
                    == Equation::Symbol{"source."+k}));
        addEquation(system,{secondEquation,
                            kEpsilon ? "epsilon transport" : "omega transport",
                            "conservation",{second}},
            Equation::named(secondEquation,
                Equation::ddt({second}) + Equation::div({"flux."+second})
                    + Equation::diffusion(
                        {"diffusivity."+second},{second})
                    == Equation::Symbol{"source."+second}));
        block.equations.insert(
            block.equations.end(),{kEquation,secondEquation});
        block.unknowns.insert(block.unknowns.end(),{k,second});
    }
    system.closures.push_back("mu_t from "+request.turbulenceModel);
    system.solveBlocks.push_back(std::move(block));
}

struct SourceEquationContribution {
    FDM::SourceKind kind;
    const char* symbol;
    bool energy;
};

void addSourceContributions(
        ResolvedSimulationSystem& system,
        const FDM::SourceConfig& config) {
    static const std::array<SourceEquationContribution,3> registry{{
        {FDM::SourceKind::Gravity,"gravity",false},
        {FDM::SourceKind::MRF,"MRF",false},
        {FDM::SourceKind::WallHeat,"wallHeat",true}
    }};
    for (FDM::SourceKind kind : config.enabled) {
        const auto contribution = std::find_if(
            registry.begin(),registry.end(),
            [kind](const SourceEquationContribution& item) {
                return item.kind == kind;
            });
        if (contribution == registry.end()) {
            throw std::runtime_error(
                "No resolved equation contribution is registered for SourceKind.");
        }
        for (const auto& equation : system.equations) {
            const bool target = contribution->energy
                ? (equation.id == "E_ENERGY"
                   || equation.id.rfind("E_ENTHALPY.",0) == 0)
                : (equation.id == "E_MOMENTUM"
                   || equation.id.rfind("E_MOMENTUM.",0) == 0);
            if (target) {
                EquationSystemBuilder(system).addTerm(
                    equation.id,Equation::source({contribution->symbol}));
            }
        }
    }
}

VariableLocation location(FDM::ImmersedVariableLocation value) {
    switch (value) {
        case FDM::ImmersedVariableLocation::EulerianGlobalDof:
            return VariableLocation::EulerianCell;
        case FDM::ImmersedVariableLocation::BodyConstraintDof:
            return VariableLocation::BodyConstraint;
        case FDM::ImmersedVariableLocation::SurfaceConstraintDof:
            return VariableLocation::SurfaceConstraint;
        case FDM::ImmersedVariableLocation::SolidGlobalDof:
            return VariableLocation::SolidGlobal;
    }
    throw std::runtime_error("Unknown immersed variable location.");
}

OwnershipKind ownership(FDM::ImmersedOwnershipKind value) {
    switch (value) {
        case FDM::ImmersedOwnershipKind::EulerianOwner:
            return OwnershipKind::EulerianGlobalDof;
        case FDM::ImmersedOwnershipKind::ConstraintOwner:
            return OwnershipKind::ConstraintGlobalDof;
        case FDM::ImmersedOwnershipKind::SolidOwner:
            return OwnershipKind::SolidGlobalDof;
    }
    throw std::runtime_error("Unknown immersed ownership kind.");
}

void addImmersed(
        ResolvedSimulationSystem& system,
        const FDM::ImmersedAlgorithmDescriptor& immersed) {
    system.immersedAlgorithm = immersed.id;
    system.immersedReference = immersed.referenceName;
    system.immersedSupport = FDM::toString(immersed.support);
    system.immersedRepresentation = FDM::toString(immersed.representation);
    system.immersedEnforcement = FDM::toString(immersed.enforcement);
    system.immersedSolid = FDM::toString(immersed.solid);
    system.immersedFunctional = immersed.variational.stationaryFunctional;
    for (const auto& source : immersed.unknowns) {
        UnknownDescriptor unknown;
        unknown.id = source.id;
        unknown.name = source.name;
        unknown.location = location(source.location);
        unknown.components = source.components;
        unknown.ownership = ownership(source.ownership);
        unknown.shape = source.components == 1
            ? ValueShape::Scalar : ValueShape::Vector;
        unknown.role = unknown.location == VariableLocation::BodyConstraint
                || unknown.location == VariableLocation::SurfaceConstraint
            ? UnknownRole::Multiplier : UnknownRole::Algebraic;
        unknown.storageBinding = StorageBinding::SpecializedExecutor;
        unknown.nameSpace = "immersed";
        system.unknowns.push_back(std::move(unknown));
    }
    for (const auto& source : immersed.constraints) {
        system.constraints.push_back({
            source.id,source.name,source.equation,source.multiplierUnknown});
    }
    for (const auto& source : immersed.equations) {
        const std::string unknown = source.solvedUnknowns.empty()
            ? source.id : source.solvedUnknowns.front();
        addEquation(system,
            {source.id,source.name,source.form,source.solvedUnknowns},
            Equation::named(source.id,
                Equation::constraint({source.form.empty()
                    ? source.id : source.form})
                    == Equation::Symbol{unknown}));
    }
    for (const auto& source : immersed.solveBlocks) {
        SolveBlock block;
        block.id = source.id;
        block.name = source.name;
        block.strategy = source.strategy;
        if (source.strategy == "boundaryStencilClosure") {
            block.strategyKind = FDM::SolveStrategyKind::BoundaryClosure;
        } else if (source.strategy == "monolithicKKT"
                   || source.strategy == "augmentedLagrangianKKT") {
            block.strategyKind = FDM::SolveStrategyKind::MonolithicKKT;
        } else if (source.strategy == "explicitLaggedMultiplier"
                   || source.strategy == "fractionalVariationalProjection"
                   || source.strategy == "surfaceMassProjection"
                   || source.strategy == "dissipativePenalty") {
            block.strategyKind = FDM::SolveStrategyKind::ConstraintSolve;
        } else {
            throw std::runtime_error(
                "Unknown immersed solve strategy '"+source.strategy+"'.");
        }
        block.equations = source.equations;
        block.unknowns = source.unknowns;
        block.constraints = source.constraints;
        if (immersed.monolithic) {
            block.unknowns.insert(block.unknowns.begin(),{"pPrime","U"});
            block.equations.insert(
                block.equations.begin(),{"E_MOMENTUM","E_PRESSURE"});
        }
        system.solveBlocks.push_back(std::move(block));
    }
}

} // namespace

ResolvedSimulationSystem build(
        const FDM::SolverConfig& config,
        const BuildRequest& request) {
    ResolvedSimulationSystem result;
    result.formulation = config.numerics.solver;
    result.templateOrigin = request.templateOrigin;
    result.timeIntegrator = timeName(config.numerics.time);

    if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
        addEulerianEulerianTemplate(result,request.phaseNames);
    } else if (config.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        addDensityBasedFluid(result);
    } else {
        addPressureBasedFluid(
            result, FDM::toString(config.pressure.workflow.algorithm));
    }

    addSourceContributions(result,config.sources);

    if (request.homogeneousThermodynamics) {
        addUnknown(result,"phaseMassAux","homogeneous phase mass auxiliary",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"homogeneous");
        addEquation(result,{
            "E_PHASE_MASS","homogeneous phase-mass transport",
            "conservation",{"phaseMassAux"}},
            Equation::named("E_PHASE_MASS",
                Equation::ddt({"phaseMassAux"})
                    + Equation::div({"phaseMassFlux"})
                    == Equation::Symbol{"phaseMassSources"}));
        attachEquationToPrimarySolveBlock(
            result,"E_PHASE_MASS","phaseMassAux");
    } else if (request.legacyMixture) {
        addUnknown(result,"alphaAux","legacy transported volume fraction",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"legacyMultiphase");
        addEquation(result,{
            "E_LEGACY_ALPHA","legacy volume-fraction transport",
            "conservation",{"alphaAux"}},
            Equation::named("E_LEGACY_ALPHA",
                Equation::ddt({"alphaAux"})
                    + Equation::div({"alphaFlux"})
                    == Equation::Symbol{"zero"}));
        attachEquationToPrimarySolveBlock(
            result,"E_LEGACY_ALPHA","alphaAux");
    }
    if (request.levelSet) {
        addUnknown(result,"phi","level-set geometry",1,
                   UnknownRole::Transported,StorageBinding::NamedDistributed,
                   "phi",0,"interface");
        addEquation(result,{
            "E_LEVEL_SET","level-set advection/reinitialization",
            "geometry transport",{"phi"}},
            Equation::named("E_LEVEL_SET",
                Equation::ddt({"phi"}) + Equation::div({"levelSetFlux"})
                    == Equation::Symbol{"reinitialization"}));
        attachEquationToPrimarySolveBlock(result,"E_LEVEL_SET","phi");
        result.closures.push_back("surface normal and curvature from phi");
        if (request.interfaceGhostFluid) {
            result.closures.push_back("ghost-fluid interface closure");
        }
    }
    addTurbulenceEquations(result,request);
    if (request.immersed) {
        addImmersed(result,*request.immersed);
        if (request.immersed->monolithic) {
            result.solveBlocks.erase(
                std::remove_if(
                    result.solveBlocks.begin(), result.solveBlocks.end(),
                    [](const SolveBlock& block) {
                        return block.id == "S_PRESSURE";
                    }),
                result.solveBlocks.end());
        }
    }

    result.requirements.push_back({
        "EulerianGlobalDof",true,true,"canonical cell ownership"});
    result.requirements.push_back({
        "CanonicalFace",request.parallel,true,
        "one owner computes each shared face flux"});
    const bool constraintDof = request.immersed
        && request.immersed->introducesMultiplier && request.parallel;
    result.requirements.push_back({
        "ConstraintGlobalDof",constraintDof,
        !constraintDof || request.constraintGlobalDofAvailable,
        "unique owner for Lambda/lambda and sparse J/S edges"});
    const bool distributedSolve = request.immersed
        && request.immersed->monolithic;
    result.requirements.push_back({
        "DistributedLinearSystem",distributedSolve,
        request.distributedLinearSystemAvailable,
        "GlobalDofId is mapped to backend rows outside equation assembly"});
    const bool homogeneousUnsupported = request.homogeneousThermodynamics
        && (config.numerics.solver != FDM::SolverAlgorithm::DensityBased
            || (config.turbulence.enabled
                && config.turbulence.family != FDM::TurbulenceFamily::DNS));
    result.requirements.push_back({
        "HomogeneousEquationExecution",request.homogeneousThermodynamics,
        !homogeneousUnsupported,
        "current homogeneous EquationSet execution requires density formulation "
        "without transported turbulence"});
    const bool eulerianExecution =
        request.templateOrigin == PhysicsTemplateKind::EulerianEulerian;
    result.requirements.push_back({
        "EulerianEquationExecution",eulerianExecution,
        !eulerianExecution
            || (config.numerics.solver == FDM::SolverAlgorithm::PressureBased
                && std::isfinite(config.numerics.maxDeltaT)
                && config.numerics.maxDeltaT > 0.0),
        "current Eulerian equation executor requires pressure formulation and "
        "a finite positive maxDeltaT"});
    const bool unsupportedPhaseChange = request.phaseChange
        && (request.levelSet || request.legacyMixture);
    result.requirements.push_back({
        "PhaseChangeExecution",request.phaseChange,
        !unsupportedPhaseChange,
        "phase change is implemented for homogeneous or Eulerian equation systems"});
    result.requirements.push_back({
        "CanonicalScalarInterfaceFlux",
        request.transportedLegacyAlpha && request.parallel,
        true,
        "legacy transported alpha on coupled patches is not implemented"});
    return result;
}

bool hasUnknown(const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(system.unknowns.begin(), system.unknowns.end(),
        [&](const UnknownDescriptor& value) { return value.id == id; });
}

bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(system.equations.begin(), system.equations.end(),
        [&](const EquationDescriptor& value) { return value.id == id; });
}

const Equation::Definition& equationDefinition(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return system.equationDefinitions.at(std::string(id));
}

bool hasEquationPrefix(
        const ResolvedSimulationSystem& system, std::string_view prefix) {
    return std::any_of(system.equations.begin(), system.equations.end(),
        [&](const EquationDescriptor& value) {
            return value.id.compare(0, prefix.size(), prefix) == 0;
        });
}

bool hasConstraint(const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(system.constraints.begin(), system.constraints.end(),
        [&](const ConstraintDescriptor& value) { return value.id == id; });
}

bool hasSolveBlock(const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(system.solveBlocks.begin(), system.solveBlocks.end(),
        [&](const SolveBlock& value) { return value.id == id; });
}

bool hasSolveStrategy(
        const ResolvedSimulationSystem& system, FDM::SolveStrategyKind kind) {
    return std::any_of(system.solveBlocks.begin(),system.solveBlocks.end(),
        [kind](const SolveBlock& block) { return block.strategyKind == kind; });
}

bool hasRequirement(
        const ResolvedSimulationSystem& system, std::string_view name) {
    return std::any_of(system.requirements.begin(), system.requirements.end(),
        [&](const ExecutionRequirement& value) { return value.name == name; });
}

bool requiresCapability(
        const ResolvedSimulationSystem& system, std::string_view name) {
    return std::any_of(system.requirements.begin(), system.requirements.end(),
        [&](const ExecutionRequirement& value) {
            return value.name == name && value.required;
        });
}

const char* toString(PhysicsTemplateKind kind) {
    switch (kind) {
        case PhysicsTemplateKind::SingleFluid: return "singleFluid";
        case PhysicsTemplateKind::HomogeneousMixture:
            return "homogeneousMixture";
        case PhysicsTemplateKind::OneFluidInterface: return "oneFluidInterface";
        case PhysicsTemplateKind::EulerianEulerian: return "eulerianEulerian";
    }
    return "unknown";
}

const char* toString(VariableLocation value) {
    switch (value) {
        case VariableLocation::EulerianCell: return "EulerianCell";
        case VariableLocation::EulerianFace: return "EulerianFace";
        case VariableLocation::BodyConstraint: return "BodyConstraint";
        case VariableLocation::SurfaceConstraint: return "SurfaceConstraint";
        case VariableLocation::SolidGlobal: return "SolidGlobal";
    }
    return "UnknownLocation";
}

const char* toString(OwnershipKind value) {
    switch (value) {
        case OwnershipKind::EulerianGlobalDof: return "EulerianGlobalDof";
        case OwnershipKind::CanonicalFace: return "CanonicalFace";
        case OwnershipKind::ConstraintGlobalDof: return "ConstraintGlobalDof";
        case OwnershipKind::SolidGlobalDof: return "SolidGlobalDof";
    }
    return "UnknownOwnership";
}

} // namespace SF::System
