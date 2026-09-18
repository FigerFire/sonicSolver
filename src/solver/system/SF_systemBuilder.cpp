/// @file SF_systemBuilder.cpp
/// @brief 显式声明单流体、多相、湍流和 IBM 对数学系统的贡献。

#include "SF_systemBuilder.h"
#include "SF_equationContribution.h"
#include "SF_solvePlan.h"
#include "SF_transformation.h"

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

bool hasRawEquationPrefix(
        const RawEquationSystem& system, std::string_view prefix) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& value) {
            return value.id.compare(0,prefix.size(),prefix) == 0;
        });
}

bool usesEquation(
        const EquationCompositionConfig& composition,
        std::string_view name) {
    return std::find(composition.equations.begin(),composition.equations.end(),name)
        != composition.equations.end();
}

std::string pressureAlgorithmName(const std::string& kind) {
    if (kind == "SIMPLE" || kind == "PISO" || kind == "PIMPLE") return kind;
    throw std::runtime_error(
        "A pressure-constraint composition requires SIMPLE, PISO, or PIMPLE.");
}

void validateSemanticComposition(const EquationCompositionConfig& composition) {
    if (!composition.declared) return;
    if (composition.equations.empty()) {
        throw std::runtime_error(
            "equations.yaml is required when semantic composition is declared.");
    }
    const bool continuity = usesEquation(composition,"Continuity");
    const bool momentum = usesEquation(composition,"Momentum");
    const bool energy = usesEquation(composition,"Energy");
    if (!continuity || !momentum) {
        throw std::runtime_error(
            "Current single-fluid composition requires Continuity and Momentum.");
    }
    const auto providers = builtinCompositionProviders();
    if (composition.thermoDynamics.equationOfState.empty()) {
        throw std::runtime_error(
            "Active fluid equations require models/thermoDynamics.yaml with equationOfState.");
    }
    if (!providers.hasEquationOfState(composition.thermoDynamics.equationOfState)) {
        throw std::runtime_error("No registered equation-of-state provider for '"
            +composition.thermoDynamics.equationOfState+"'.");
    }
    if (energy && composition.thermoDynamics.thermo.empty()) {
        throw std::runtime_error(
            "Energy equation requires caloric thermodynamics, but no thermo model is registered.");
    }
    if (energy && !providers.hasCaloricThermo(composition.thermoDynamics.thermo)) {
        throw std::runtime_error("No registered caloric-thermo provider for '"
            +composition.thermoDynamics.thermo+"'.");
    }
    if (composition.thermoDynamics.equationOfState == "rhoConst") {
        if (energy) throw std::runtime_error(
            "rhoConst composition currently supports Momentum + Continuity; "
            "select an explicit Energy/Enthalpy equation before adding caloric closure.");
        if (!providers.hasTransport(composition.thermoDynamics.transport)) {
            throw std::runtime_error(
                "rhoConst Momentum requires transport; no transport model is registered.");
        }
        (void)pressureAlgorithmName(composition.algorithm);
    } else if (composition.thermoDynamics.equationOfState == "perfectGas") {
        if (!energy) throw std::runtime_error(
            "perfectGas single-fluid execution currently requires Energy; "
            "a non-caloric variable-density equation pack is unsupported.");
        if (!composition.algorithm.empty() && composition.algorithm != "Explicit") {
            throw std::runtime_error(
                "The selected pressure algorithm requires an incompressibility constraint; "
                "the perfectGas conservative system has none.");
        }
    }
}

void addUnknown(
        SystemCompositionBuilder& system,
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
    system.addUnknown(std::move(unknown));
}

void addEquation(
        SystemCompositionBuilder& system,
        EquationDescriptor descriptor,
        Equation::Definition definition) {
    system.addEquation(std::move(descriptor),std::move(definition));
}

void addDensityBasedFluid(
        SystemCompositionBuilder& system,
        const std::string& timeIntegrator) {
    system.recordContribution("builtin.navierStokes","core Navier-Stokes preset");
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
    ExecutionPolicy policy;
    policy.id = "S_FLUID";
    policy.name = "conservative explicit equations";
    policy.kind = ExecutionPolicyKind::ExplicitStages;
    policy.strategyName = timeIntegrator;
    policy.strategyKind = FDM::SolveStrategyKind::ExplicitTimeIntegration;
    policy.equations = {"E_MASS","E_MOMENTUM","E_ENERGY"};
    policy.unknowns = {"rho","rhoU","rhoE"};
    system.addExecutionPolicy(std::move(policy));
}

void addPressureBasedFluid(
        SystemCompositionBuilder& system,
        const std::string& timeIntegrator,
        const FDM::SolverPropertiesConfig& pressure) {
    const std::string pressureStrategy = FDM::toString(pressure.algorithm);
    system.recordContribution("builtin.navierStokes","core Navier-Stokes preset");
    system.recordContribution(
        "preset.coupling."+pressureStrategy,
        pressureStrategy+" pressure-constraint preset");
    addUnknown(system,"rho","density",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",0,"fluid");
    addUnknown(system,"U","velocity",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"rhoE","total energy",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",4,"fluid");
    addUnknown(system,"p","thermodynamic pressure",1,
               UnknownRole::Derived,StorageBinding::SpecializedExecutor,
               "thermodynamicClosure",0,"pressure");
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
    ConstraintDescriptor pressureConstraint{
        "C_INCOMPRESSIBILITY","pressure/continuity constraint",
        "pressureVelocityConsistency = 0","p"};
    system.addConstraint(std::move(pressureConstraint));
    system.requestTransformation({
        "pressureConstraint",pressureStrategy+" pressure transformation",
        100,true,{}});

    ExecutionPolicy predictor;
    predictor.id = "S_PREDICTOR";
    predictor.name = "conservative predictor";
    predictor.kind = ExecutionPolicyKind::ExplicitStages;
    predictor.strategyName = timeIntegrator;
    predictor.strategyKind = FDM::SolveStrategyKind::SegregatedPredictor;
    // The current production predictor is a fused conservative update: the
    // pressure transformation names its momentum role explicitly, while the
    // same operation also advances mass and energy with the existing kernel.
    predictor.equations = {
        "E_MASS","E_MOMENTUM","E_ENERGY","E_MOMENTUM_PREDICTOR"};
    predictor.unknowns = {"rho","U","rhoE"};
    predictor.priority = 10;
    system.addExecutionPolicy(std::move(predictor));

    ExecutionPolicy correction;
    correction.id = "S_PRESSURE";
    correction.name = "pressure-velocity correction";
    correction.kind = pressureStrategy == "PIMPLE"
        ? ExecutionPolicyKind::PressureVelocityFixedPoint
        : ExecutionPolicyKind::SegregatedPressureCorrection;
    correction.strategyName = pressureStrategy;
    correction.strategyKind = FDM::SolveStrategyKind::PressureCorrection;
    correction.equations = {"E_PRESSURE"};
    correction.constraints = {"C_INCOMPRESSIBILITY"};
    correction.unknowns = {"pPrime","U"};
    correction.priority = 20;
    correction.repeatCount = pressure.pressureCorrectors;
    correction.nestedRepeatCount = pressure.nonOrthogonalCorrectors+1;
    system.addExecutionPolicy(std::move(correction));
}

/// @brief Raw U/p composition for rho=rho0.  It intentionally does not
/// fabricate rhoE or a PerfectGas closure: execution is enabled only after a
/// dedicated constant-density predictor/corrector operator is bound.
void addConstantDensityFluid(
        SystemCompositionBuilder& system,
        const EquationCompositionConfig& composition,
        const std::string& timeIntegrator) {
    const std::string algorithm = pressureAlgorithmName(composition.algorithm);
    system.recordContribution("builtin.continuity","builtin Continuity equation");
    system.recordContribution("builtin.momentum","builtin Momentum equation");
    system.recordContribution("model.rhoConst","constant-density EOS closure");
    system.recordContribution("preset.coupling."+algorithm,
                              algorithm+" pressure-constraint preset");
    addUnknown(system,"rho","constant density closure",1,UnknownRole::Derived,
               StorageBinding::SpecializedExecutor,"rhoConst",0,"fluid");
    addUnknown(system,"U","velocity",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"p","pressure multiplier",1,UnknownRole::Multiplier,
               StorageBinding::NamedDistributed,"pressure",0,"pressure");
    addEquation(system,{"E_CONTINUITY","continuity","constraint",{"U"}},
        Equation::named("E_CONTINUITY",
            Equation::div({"U"}) == Equation::Symbol{"zero"}));
    addEquation(system,{"E_MOMENTUM","momentum","conservation",{"U"}},
        Equation::named("E_MOMENTUM",
            Equation::ddt({"U"}) + Equation::div({"momentumFlux"})
                + Equation::diffusion({"nu"},{"U"})
                == Equation::Symbol{"zero"}));
    system.addConstraint({"C_INCOMPRESSIBILITY","constant-density continuity",
                          "div(U) = 0","p"});
    system.requestTransformation({"pressureConstraint",
                                  algorithm+" pressure transformation",100,true,{}});

    ExecutionPolicy predictor;
    predictor.id = "S_PREDICTOR";
    predictor.name = "constant-density momentum predictor";
    predictor.kind = ExecutionPolicyKind::ExplicitStages;
    predictor.strategyName = timeIntegrator;
    predictor.strategyKind = FDM::SolveStrategyKind::SegregatedPredictor;
    predictor.equations = {
        "E_CONTINUITY","E_MOMENTUM","E_MOMENTUM_PREDICTOR"};
    predictor.unknowns = {"U"};
    predictor.priority = 10;
    system.addExecutionPolicy(std::move(predictor));

    ExecutionPolicy correction;
    correction.id = "S_PRESSURE";
    correction.name = "constant-density pressure correction";
    correction.kind = algorithm == "PIMPLE"
        ? ExecutionPolicyKind::PressureVelocityFixedPoint
        : ExecutionPolicyKind::SegregatedPressureCorrection;
    correction.strategyName = algorithm;
    correction.strategyKind = FDM::SolveStrategyKind::PressureCorrection;
    correction.equations = {"E_PRESSURE"};
    correction.constraints = {"C_INCOMPRESSIBILITY"};
    correction.unknowns = {"pPrime","U"};
    correction.priority = 20;
    correction.repeatCount = composition.pressureCorrectors;
    correction.nestedRepeatCount = composition.nonOrthogonalCorrectors + 1;
    system.addExecutionPolicy(std::move(correction));
}

void addPhaseEquationPack(
        SystemCompositionBuilder& system,
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

void addSharedPressureConstraint(SystemCompositionBuilder& system) {
    system.addConstraint({
        "C_SHARED_PRESSURE","shared pressure relationship",
        "p.phase = p",""});
    system.addConstraint({
        "C_VOLUME_FRACTION","volume-fraction closure",
        "sum(alpha.phase) = 1",""});
    system.requestTransformation({
        "sharedPressureConstraint","PIMPLE shared-pressure transformation",
        110,true,{}});
}

void addEulerianEulerianSolvePolicy(
        SystemCompositionBuilder& system, const FDM::SolverConfig& config) {
    ExecutionPolicy policy;
    policy.id = "S_EE_PIMPLE";
    policy.name = "Eulerian-Eulerian pressure coupling";
    policy.kind = ExecutionPolicyKind::PressureVelocityFixedPoint;
    policy.strategyName = "PIMPLE phase predictor/corrector";
    policy.strategyKind = FDM::SolveStrategyKind::PressureVelocityCoupling;
    policy.constraints = {"C_SHARED_PRESSURE","C_VOLUME_FRACTION"};
    policy.unknowns = {"p"};
    for (const auto& equation : system.rawSystem().equations) {
        policy.equations.push_back(equation.id);
        for (const auto& unknown : equation.solvedUnknowns) {
            if (unknown != "p") {
                policy.unknowns.push_back(unknown);
            }
        }
    }
    // Algorithmic pressure equation is generated after raw composition.
    policy.equations.push_back("E_SHARED_PRESSURE");
    policy.repeatCount = config.pressure.workflow.outerCorrectors;
    policy.nestedRepeatCount = config.pressure.workflow.pressureCorrectors;
    // The existing stepper executes zero plus N non-orthogonal passes.
    policy.innerRepeatCount = config.pressure.workflow.nonOrthogonalCorrectors + 1;
    system.addExecutionPolicy(std::move(policy));
}

void addEulerianEulerianTemplate(
        SystemCompositionBuilder& system,
        ResolvedSimulationSystem& resolved,
        const std::vector<std::string>& names, const FDM::SolverConfig& config) {
    if (names.size() < 2) {
        throw std::runtime_error(
            "Resolved Eulerian-Eulerian system requires at least two phases.");
    }
    system.recordContribution(
        "preset.eulerianEulerian","Eulerian-Eulerian equation preset");
    system.recordContribution(
        "preset.coupling.PIMPLE","PIMPLE shared-pressure preset");
    addUnknown(system,"p","shared pressure",1,UnknownRole::Algebraic,
               StorageBinding::NamedDistributed,"pressure",0,"pressure");
    for (std::size_t phase = 0; phase < names.size(); ++phase) {
        addPhaseEquationPack(system,names[phase],phase);
        const std::string prefix = "phase"+std::to_string(phase)+".";
        resolved.workspaceRequirements.push_back({
            prefix+"momentumDiagonal",1,VariableLocation::EulerianCell,
            OwnershipKind::EulerianGlobalDof});
        resolved.workspaceRequirements.push_back({
            prefix+"volumeFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
        resolved.workspaceRequirements.push_back({
            prefix+"massFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
    }
    resolved.workspaceRequirements.push_back({
        "pressureCorrection",1,VariableLocation::EulerianCell,
        OwnershipKind::EulerianGlobalDof});
    addSharedPressureConstraint(system);
    addEulerianEulerianSolvePolicy(system,config);
}

void attachEquationToPrimarySolveBlock(
        SystemCompositionBuilder& system,
        const std::string& equation,
        const std::string& unknown) {
    const char* policy = hasRawEquationPrefix(system.rawSystem(),"E_CONTINUITY.")
        ? "S_EE_PIMPLE"
        : (hasConstraint(system.rawSystem(),"C_INCOMPRESSIBILITY")
            ? "S_PREDICTOR" : "S_FLUID");
    try {
        system.extendExecutionPolicy(policy,equation,unknown);
    } catch (const std::runtime_error&) {
        throw std::runtime_error(
            "Auxiliary equation '"+equation
            +"' has no primary fluid solve block.");
    }
}

void addTurbulenceEquations(
        SystemCompositionBuilder& system,
        const BuildRequest& request) {
    if (!request.turbulence) return;
    system.recordContribution(
        "model.turbulence","turbulence equation/closure contribution");
    const bool kEpsilon = request.turbulenceModel == "kEpsilon";
    const bool kOmega = request.turbulenceModel == "kOmegaSST";
    if (!kEpsilon && !kOmega) {
        system.addClosure(
            request.turbulenceModel.empty()
                ? "turbulence closure"
                : "turbulence: "+request.turbulenceModel);
        return;
    }

    ExecutionPolicy policy;
    policy.id = "S_TURBULENCE";
    policy.name = "transported turbulence equations";
    policy.kind = ExecutionPolicyKind::ExplicitStages;
    policy.strategyName = "transported turbulence legacy schedule";
    policy.strategyKind = FDM::SolveStrategyKind::ExplicitTimeIntegration;
    policy.priority = 30;
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
        if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
            // PIMPLE owns the one global schedule.  Turbulence contributes
            // equations to that schedule instead of creating another stage loop.
            system.extendExecutionPolicy("S_EE_PIMPLE",kEquation,k);
            system.extendExecutionPolicy("S_EE_PIMPLE",secondEquation,second);
        } else {
            policy.equations.insert(
                policy.equations.end(),{kEquation,secondEquation});
            policy.unknowns.insert(policy.unknowns.end(),{k,second});
        }
    }
    system.addClosure("mu_t from "+request.turbulenceModel);
    if (request.templateOrigin != PhysicsTemplateKind::EulerianEulerian) {
        system.addExecutionPolicy(std::move(policy));
    }
}

struct SourceEquationContribution {
    FDM::SourceKind kind;
    const char* symbol;
    bool energy;
};

void addSourceContributions(
        SystemCompositionBuilder& system,
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
        for (const auto& equation : system.rawSystem().equations) {
            const bool target = contribution->energy
                ? (equation.id == "E_ENERGY"
                   || equation.id.rfind("E_ENTHALPY.",0) == 0)
                : (equation.id == "E_MOMENTUM"
                   || equation.id.rfind("E_MOMENTUM.",0) == 0);
            if (target) {
                system.extendEquation(
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
        SystemCompositionBuilder& system,
        ResolvedSimulationSystem& resolved,
        const FDM::ImmersedAlgorithmDescriptor& immersed) {
    system.recordContribution(
        "model.ibm."+immersed.id,"immersed-boundary contribution");
    resolved.immersedAlgorithm = immersed.id;
    resolved.immersedReference = immersed.referenceName;
    resolved.immersedSupport = FDM::toString(immersed.support);
    resolved.immersedRepresentation = FDM::toString(immersed.representation);
    resolved.immersedEnforcement = FDM::toString(immersed.enforcement);
    resolved.immersedSolid = FDM::toString(immersed.solid);
    resolved.immersedFunctional = immersed.variational.stationaryFunctional;
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
        system.addUnknown(std::move(unknown));
    }
    for (const auto& source : immersed.constraints) {
        system.addConstraint({
            source.id,source.name,source.equation,source.multiplierUnknown});
    }
    for (const auto& source : immersed.equations) {
        const std::string unknown = source.solvedUnknowns.empty()
            ? source.id : source.solvedUnknowns.front();
        EquationDescriptor descriptor{
            source.id,source.name,source.form,source.solvedUnknowns};
        descriptor.category = EquationCategory::ConstraintEquation;
        addEquation(system,std::move(descriptor),
            Equation::named(source.id,
                Equation::constraint({source.form.empty()
                    ? source.id : source.form})
                    == Equation::Symbol{unknown}));
    }
    bool hasConstraintTransformation = false;
    for (const auto& source : immersed.solveBlocks) {
        ExecutionPolicy policy;
        policy.id = source.id;
        policy.name = source.name;
        policy.strategyName = source.strategy;
        if (source.strategy == "boundaryStencilClosure") {
            policy.kind = ExecutionPolicyKind::BoundaryClosure;
            policy.strategyKind = FDM::SolveStrategyKind::BoundaryClosure;
        } else if (source.strategy == "monolithicKKT"
                   || source.strategy == "augmentedLagrangianKKT") {
            policy.kind = ExecutionPolicyKind::MonolithicKKT;
            policy.strategyKind = FDM::SolveStrategyKind::MonolithicKKT;
            hasConstraintTransformation = true;
        } else if (source.strategy == "explicitLaggedMultiplier"
                   || source.strategy == "fractionalVariationalProjection"
                   || source.strategy == "surfaceMassProjection"
                   || source.strategy == "dissipativePenalty") {
            policy.kind = ExecutionPolicyKind::ConstraintProjection;
            policy.strategyKind = FDM::SolveStrategyKind::ConstraintSolve;
            // Brinkman/BP is a direct momentum penalty contribution and has
            // no multiplier constraint for the transformer to rewrite.
            hasConstraintTransformation = hasConstraintTransformation
                || source.strategy != "dissipativePenalty";
        } else {
            throw std::runtime_error(
                "Unknown immersed solve strategy '"+source.strategy+"'.");
        }
        policy.equations = source.equations;
        policy.unknowns = source.unknowns;
        policy.constraints = source.constraints;
        policy.priority = 200;
        if (immersed.monolithic) {
            policy.unknowns.insert(policy.unknowns.begin(),{"pPrime","U"});
            policy.equations.insert(
                policy.equations.begin(),{"E_MOMENTUM","E_PRESSURE"});
        }
        system.addExecutionPolicy(std::move(policy));
    }
    if (hasConstraintTransformation) {
        system.requestTransformation({
            "immersedConstraint","immersed constraint transformation",
            200,true,{}});
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
    validateSemanticComposition(request.composition);

    std::vector<TransformationDescriptor> transformationRequests;
    SystemCompositionBuilder builtin(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {OriginKind::BuiltinPreset,"default system"});

    if (request.composition.declared) {
        if (request.templateOrigin != PhysicsTemplateKind::SingleFluid) {
            throw std::runtime_error(
                "Semantic single-fluid equation composition cannot be combined "
                "with the current multiphase template.");
        }
        // Equations choose the raw mathematical system. EOS contributes a
        // closure and may simplify that system; it never selects a solver
        // family or swaps in a different equation pack.
        if (usesEquation(request.composition,"Energy")) {
            addDensityBasedFluid(builtin,result.timeIntegrator);
        } else {
            addConstantDensityFluid(builtin,request.composition,result.timeIntegrator);
        }
        const auto& eos = request.composition.thermoDynamics.equationOfState;
        if (eos == "rhoConst") {
            result.densityBehavior = "constant";
            result.thermodynamicCompressibility = "zero";
            result.rawSystem.closures.push_back("rho = rho0 (rhoConst)");
        } else if (eos == "perfectGas") {
            result.densityBehavior = "variable";
            result.thermodynamicCompressibility = "available";
            result.rawSystem.closures.push_back("p = p(rho,T) (perfectGas)");
        } else {
            throw std::runtime_error("Semantic composition has no EOS provider.");
        }
    } else if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
        addEulerianEulerianTemplate(builtin,result,request.phaseNames,config);
    } else if (config.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        addDensityBasedFluid(builtin,result.timeIntegrator);
    } else {
        addPressureBasedFluid(
            builtin,result.timeIntegrator,config.pressure.workflow);
    }

    if (result.densityBehavior.empty()) {
        result.densityBehavior = "legacy-configured";
        result.thermodynamicCompressibility = "legacy-configured";
    }

    SystemCompositionBuilder models(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {OriginKind::Model,"configured models"});
    addSourceContributions(models,config.sources);

    if (request.homogeneousThermodynamics) {
        models.recordContribution(
            "model.homogeneousMultiphase","homogeneous multiphase equations");
        addUnknown(models,"phaseMassAux","homogeneous phase mass auxiliary",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"homogeneous");
        addEquation(models,{
            "E_PHASE_MASS","homogeneous phase-mass transport",
            "conservation",{"phaseMassAux"}},
            Equation::named("E_PHASE_MASS",
                Equation::ddt({"phaseMassAux"})
                    + Equation::div({"phaseMassFlux"})
                    == Equation::Symbol{"phaseMassSources"}));
        attachEquationToPrimarySolveBlock(
            models,"E_PHASE_MASS","phaseMassAux");
    } else if (request.legacyMixture) {
        models.recordContribution(
            "model.legacyMixture","legacy mixture equations");
        addUnknown(models,"alphaAux","legacy transported volume fraction",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"legacyMultiphase");
        addEquation(models,{
            "E_LEGACY_ALPHA","legacy volume-fraction transport",
            "conservation",{"alphaAux"}},
            Equation::named("E_LEGACY_ALPHA",
                Equation::ddt({"alphaAux"})
                    + Equation::div({"alphaFlux"})
                    == Equation::Symbol{"zero"}));
        attachEquationToPrimarySolveBlock(
            models,"E_LEGACY_ALPHA","alphaAux");
    }
    if (request.levelSet) {
        models.recordContribution("model.levelSet","level-set equations");
        addUnknown(models,"phi","level-set geometry",1,
                   UnknownRole::Transported,StorageBinding::NamedDistributed,
                   "phi",0,"interface");
        addEquation(models,{
            "E_LEVEL_SET","level-set advection/reinitialization",
            "geometry transport",{"phi"}},
            Equation::named("E_LEVEL_SET",
                Equation::ddt({"phi"}) + Equation::div({"levelSetFlux"})
                    == Equation::Symbol{"reinitialization"}));
        attachEquationToPrimarySolveBlock(models,"E_LEVEL_SET","phi");
        models.addClosure("surface normal and curvature from phi");
        if (request.interfaceGhostFluid) {
            models.addClosure("ghost-fluid interface closure");
        }
    }
    addTurbulenceEquations(models,request);
    if (request.immersed) {
        if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
            if (request.immersed->enforcement
                == FDM::IBMEnforcement::GhostCell) {
                throw std::runtime_error(
                    "Eulerian multiphase + Ghost IBM is unsupported: "
                    "phase-wise ghost-state boundary closure is unavailable.");
            }
            throw std::runtime_error(
                "Eulerian multiphase IBM constraint exists, but required "
                "phase-wise IBM fluid-port assembly is unavailable.");
        }
        addImmersed(models,result,*request.immersed);
        if (request.immersed->monolithic) {
            result.executionPolicies.erase(
                std::remove_if(
                    result.executionPolicies.begin(),result.executionPolicies.end(),
                    [](const ExecutionPolicy& policy) {
                        return policy.id == "S_PRESSURE";
                    }),
                result.executionPolicies.end());
        }
    }

    TransformerRegistry transformers;
    transformers.registerTransformer(makePressureConstraintTransformer());
    transformers.registerTransformer(makeSharedPressureTransformer());
    transformers.registerTransformer(makeImmersedConstraintTransformer());
    result.executableSystem = TransformationPipeline::apply(
        result.rawSystem,transformationRequests,transformers,
        result.executionPolicies,result.transformations);
    result.solvePlan = SolvePlanner::compile(
        result.executableSystem,result.executionPolicies,result.timeIntegrator);
    result.executionCapabilities = SolvePlanner::capabilities(
        result.executableSystem,result.executionPolicies);
    // Runtime availability is a compilation report, never a property of the
    // plan IR.  The plan remains a complete mathematical control-flow object
    // even when its numerical provider is intentionally unavailable.
    result.runtime.requiredOperations =
        SolvePlanner::requiredOperations(result.solvePlan);
    const bool hasOperations = !result.runtime.requiredOperations.empty();
    const auto hasLegacyControl = [](const SolvePlanNode& root) {
        const auto visit = [&](const auto& self,
                               const SolvePlanNode& node) -> bool {
            if (node.kind == PlanNodeKind::Subcycle) return true;
            return std::any_of(
                node.children.begin(),node.children.end(),
                [&](const SolvePlanNode& child) { return self(self,child); });
        };
        return visit(visit,root);
    };
    if (result.executionCapabilities.constantDensity
        && result.executionCapabilities.pressureConstraint) {
        result.runtime.status = RuntimeStatus::Unsupported;
        result.runtime.requiredAdapters.push_back(
            "UnsupportedConstantDensityPressureExecution");
        result.runtime.reason =
            "Constant-density pressure execution has no numerical operator provider.";
    } else if (hasOperations && !hasLegacyControl(result.solvePlan.root)) {
        result.runtime.status = RuntimeStatus::Runnable;
    } else if (hasEquation(result.executableSystem,"E_PRESSURE")) {
        result.runtime.requiredAdapters.push_back("LegacyPressureExecutionAdapter");
    } else {
        result.runtime.requiredAdapters.push_back("LegacyDensityExecutionAdapter");
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
    return hasUnknown(system.executableSystem,id);
}

bool hasUnknown(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.unknowns.begin(),system.unknowns.end(),
        [&](const UnknownDescriptor& value) { return value.id == id; });
}

bool hasUnknown(const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.unknowns.begin(),system.unknowns.end(),
        [&](const UnknownDescriptor& value) { return value.id == id; });
}

bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasEquation(system.executableSystem,id);
}

bool hasEquation(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& value) { return value.id == id; });
}

bool hasEquation(const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& value) { return value.id == id; });
}

const Equation::Definition& equationDefinition(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return system.executableSystem.equationDefinitions.at(std::string(id));
}

bool hasEquationPrefix(
        const ResolvedSimulationSystem& system, std::string_view prefix) {
    return std::any_of(
        system.executableSystem.equations.begin(),
        system.executableSystem.equations.end(),
        [&](const EquationDescriptor& value) {
            return value.id.compare(0, prefix.size(), prefix) == 0;
        });
}

bool hasConstraint(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasConstraint(system.executableSystem,id);
}

bool hasConstraint(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.constraints.begin(),system.constraints.end(),
        [&](const ConstraintDescriptor& value) { return value.id == id; });
}

bool hasConstraint(
        const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.constraints.begin(),system.constraints.end(),
        [&](const ConstraintDescriptor& value) { return value.id == id; });
}

bool hasSolveBlock(const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(
        system.solvePlan.blocks.begin(),system.solvePlan.blocks.end(),
        [&](const SolveBlock& value) { return value.id == id; });
}

bool hasSolveStrategy(
        const ResolvedSimulationSystem& system, FDM::SolveStrategyKind kind) {
    return std::any_of(
        system.solvePlan.blocks.begin(),system.solvePlan.blocks.end(),
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

const char* toString(OriginKind value) {
    switch (value) {
        case OriginKind::BuiltinDefault: return "builtin-default";
        case OriginKind::BuiltinPreset: return "builtin-preset";
        case OriginKind::Model: return "model";
        case OriginKind::User: return "user";
        case OriginKind::Generated: return "generated";
    }
    return "unknown-origin";
}

const char* toString(ModificationKind value) {
    switch (value) {
        case ModificationKind::Add: return "add";
        case ModificationKind::Extend: return "extend";
        case ModificationKind::Replace: return "replace";
        case ModificationKind::Disable: return "disable";
    }
    return "unknown-modification";
}

const char* toString(EquationCategory value) {
    switch (value) {
        case EquationCategory::PhysicalEquation: return "physical";
        case EquationCategory::ConstraintEquation: return "constraint";
        case EquationCategory::AlgorithmicDerivedEquation:
            return "algorithmic-derived";
        case EquationCategory::AlgebraicRelation: return "algebraic-relation";
    }
    return "unknown-equation-category";
}

const char* toString(TransformationState value) {
    switch (value) {
        case TransformationState::NotRegistered: return "not-registered";
        case TransformationState::RegisteredAndActive:
            return "registered-active";
        case TransformationState::RegisteredButNotApplicable:
            return "registered-not-applicable";
        case TransformationState::RegisteredButInvalid:
            return "registered-invalid";
        case TransformationState::Applied: return "applied";
    }
    return "unknown-transformation-state";
}

const char* toString(ExecutionPolicyKind value) {
    switch (value) {
        case ExecutionPolicyKind::ExplicitStages: return "explicit-stages";
        case ExecutionPolicyKind::SegregatedPressureCorrection:
            return "segregated-pressure-correction";
        case ExecutionPolicyKind::PressureVelocityFixedPoint:
            return "pressure-velocity-fixed-point";
        case ExecutionPolicyKind::ConstraintProjection:
            return "constraint-projection";
        case ExecutionPolicyKind::MonolithicKKT: return "monolithic-kkt";
        case ExecutionPolicyKind::BoundaryClosure: return "boundary-closure";
    }
    return "unknown-execution-policy";
}

const char* toString(PlanNodeKind value) {
    switch (value) {
        case PlanNodeKind::Sequence: return "Sequence";
        case PlanNodeKind::Loop: return "Loop";
        case PlanNodeKind::StageLoop: return "StageLoop";
        case PlanNodeKind::Subcycle: return "Subcycle";
        case PlanNodeKind::BlockSolve: return "BlockSolve";
        case PlanNodeKind::Assemble: return "Assemble";
        case PlanNodeKind::Solve: return "Solve";
        case PlanNodeKind::Correct: return "Correct";
        case PlanNodeKind::Update: return "Update";
        case PlanNodeKind::Synchronize: return "Synchronize";
        case PlanNodeKind::Reduction: return "Reduction";
        case PlanNodeKind::Commit: return "Commit";
        case PlanNodeKind::ConvergenceCheck: return "ConvergenceCheck";
    }
    return "UnknownPlanNode";
}

const char* toString(ResourceAccessMode value) {
    switch (value) {
        case ResourceAccessMode::Read: return "read";
        case ResourceAccessMode::Write: return "write";
        case ResourceAccessMode::ReadWrite: return "read-write";
    }
    return "unknown-access";
}

const char* toString(SynchronizationRequirement value) {
    switch (value) {
        case SynchronizationRequirement::None: return "none";
        case SynchronizationRequirement::ReadHalo: return "read-halo";
        case SynchronizationRequirement::WriteOwned: return "write-owned";
    }
    return "unknown-synchronization";
}

} // namespace SF::System
