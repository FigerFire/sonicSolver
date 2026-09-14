/// @file SF_systemBuilder.cpp
/// @brief 显式声明单流体、多相、湍流和 IBM 对数学系统的贡献。

#include "SF_systemBuilder.h"

#include "SF_config.h"

#include <algorithm>
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
        int components = 1) {
    system.unknowns.push_back({
        std::move(id),std::move(name),VariableLocation::EulerianCell,
        components,OwnershipKind::EulerianGlobalDof});
}

void addDensityBasedFluid(ResolvedSimulationSystem& system) {
    addUnknown(system,"rho","density");
    addUnknown(system,"rhoU","momentum",3);
    addUnknown(system,"rhoE","total energy");
    system.equations.push_back({"E_MASS","mass","conservation",{"rho"}});
    system.equations.push_back(
        {"E_MOMENTUM","momentum","conservation",{"rhoU"}});
    system.equations.push_back(
        {"E_ENERGY","total energy","conservation",{"rhoE"}});
    system.solveBlocks.push_back({
        "S_FLUID","compressible conservative predictor",
        system.timeIntegrator,
        {"E_MASS","E_MOMENTUM","E_ENERGY"},{},
        {"rho","rhoU","rhoE"}});
}

void addPressureBasedFluid(
        ResolvedSimulationSystem& system,
        const std::string& pressureStrategy) {
    addUnknown(system,"rho","density");
    addUnknown(system,"U","velocity",3);
    addUnknown(system,"rhoE","total energy");
    addUnknown(system,"pPrime","pressure correction");
    system.equations.push_back({"E_MASS","mass","conservation",{"rho"}});
    system.equations.push_back(
        {"E_MOMENTUM","momentum predictor","conservation",{"U"}});
    system.equations.push_back(
        {"E_ENERGY","total energy","conservation",{"rhoE"}});
    system.equations.push_back(
        {"E_PRESSURE","pressure correction","constraint",{"pPrime","U"}});
    system.solveBlocks.push_back({
        "S_PREDICTOR","compressible predictor",system.timeIntegrator,
        {"E_MASS","E_MOMENTUM","E_ENERGY"},{},
        {"rho","U","rhoE"}});
    system.solveBlocks.push_back({
        "S_PRESSURE","pressure-velocity correction",pressureStrategy,
        {"E_PRESSURE"},{},{"pPrime","U"}});
}

void addEulerianEulerian(
        ResolvedSimulationSystem& system,
        const std::vector<std::string>& names) {
    if (names.size() < 2) {
        throw std::runtime_error(
            "Resolved Eulerian-Eulerian system requires at least two phases.");
    }
    addUnknown(system,"p","shared pressure");
    for (const std::string& phase : names) {
        const std::string suffix = "."+phase;
        addUnknown(system,"phaseMass"+suffix,"phase mass "+phase);
        addUnknown(system,"momentum"+suffix,"phase momentum "+phase,3);
        addUnknown(system,"enthalpy"+suffix,"phase total enthalpy "+phase);
        system.equations.push_back({
            "E_CONTINUITY"+suffix,"phase continuity "+phase,
            "conservation",{"phaseMass"+suffix}});
        system.equations.push_back({
            "E_MOMENTUM"+suffix,"phase momentum "+phase,
            "conservation",{"momentum"+suffix}});
        system.equations.push_back({
            "E_ENTHALPY"+suffix,"phase enthalpy "+phase,
            "conservation",{"enthalpy"+suffix}});
    }
    system.equations.push_back({
        "E_SHARED_PRESSURE","shared pressure correction","constraint",{"p"}});
    system.solveBlocks.push_back({
        "S_EE_PIMPLE","Eulerian-Eulerian pressure coupling",
        "PIMPLE phase predictor/corrector",{}, {}, {"p"}});
    for (const auto& equation : system.equations) {
        system.solveBlocks.back().equations.push_back(equation.id);
        for (const auto& unknown : equation.solvedUnknowns) {
            if (unknown != "p") {
                system.solveBlocks.back().unknowns.push_back(unknown);
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
        system.unknowns.push_back({
            source.id,source.name,location(source.location),
            source.components,ownership(source.ownership)});
    }
    for (const auto& source : immersed.constraints) {
        system.constraints.push_back({
            source.id,source.name,source.equation,source.multiplierUnknown});
    }
    for (const auto& source : immersed.equations) {
        system.equations.push_back({
            source.id,source.name,source.form,source.solvedUnknowns});
    }
    for (const auto& source : immersed.solveBlocks) {
        SolveBlock block;
        block.id = source.id;
        block.name = source.name;
        block.strategy = source.strategy;
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
    result.flow = config.numerics.solver;
    result.physics = request.physics;
    result.timeIntegrator = timeName(config.numerics.time);

    if (request.physics == PhysicsStateKind::EulerianEulerian) {
        addEulerianEulerian(result,request.phaseNames);
    } else if (config.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        addDensityBasedFluid(result);
    } else {
        addPressureBasedFluid(
            result, FDM::toString(config.pressure.workflow.algorithm));
    }

    if (request.physics == PhysicsStateKind::HomogeneousMixture) {
        addUnknown(result,"phaseMassAux","homogeneous phase mass auxiliary");
        result.equations.push_back({
            "E_PHASE_MASS","homogeneous phase-mass transport",
            "conservation",{"phaseMassAux"}});
    }
    if (request.levelSet) {
        addUnknown(result,"phi","level-set geometry");
        result.equations.push_back({
            "E_LEVEL_SET","level-set advection/reinitialization",
            "geometry transport",{"phi"}});
        result.closures.push_back("surface normal and curvature from phi");
    }
    if (request.turbulence) {
        result.closures.push_back(
            request.turbulenceModel.empty()
                ? "transported turbulence equations"
                : "turbulence: "+request.turbulenceModel);
    }
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
    return result;
}

const char* toString(PhysicsStateKind kind) {
    switch (kind) {
        case PhysicsStateKind::SingleFluid: return "singleFluid";
        case PhysicsStateKind::HomogeneousMixture:
            return "homogeneousMixture";
        case PhysicsStateKind::OneFluidInterface: return "oneFluidInterface";
        case PhysicsStateKind::EulerianEulerian: return "eulerianEulerian";
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
