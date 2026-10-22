#include "SF_levelSetSystemContribution.h"
#include <utility>

namespace SF::Physics::InterfaceModels::LevelSetContribution {
void contribute(System::SystemContribution& system,const Spec& spec) {
    using namespace System;
    using E=FormulaExpr;
    system.recordContribution("model.levelSet","Hamilton-Jacobi transport and interface geometry");
    system.requireProvider("equation.level-set","bind phi, geometry and stage assembly ports");
    system.requireState("U");
    for (const auto& item:std::vector<std::pair<std::string,std::string>>{
            {"phi","phi"},{"interfaceCurvature","levelSetCurvature"},{"interfaceNormal","levelSetNormal"}}) {
        StateSymbol state;state.id=state.name=item.first;state.storageKey=item.second;
        state.nameSpace="interface";state.storageBinding=StorageBinding::NamedDistributed;
        state.role=item.first=="phi" ? StateRole::Transported : StateRole::Derived;
        if (item.first!="phi") {state.evaluation=StateEvaluation::Materialized;state.dependencies={"phi"};}
        if (item.first=="interfaceNormal") {state.components=3;state.shape=ValueShape::Vector;}
        state.initializationRequired=state.boundaryRequired=item.first=="phi";
        state.restartEligible=false;
        state.signature=Physics::InterfaceModels::LevelSetContribution::Mathematics::signature(state.id);system.addState(std::move(state));
    }
    system.addEquation(System::Equation{"levelSetAdvection",E::add(E::op(Math::Op::Ddt,{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")}),
        Physics::InterfaceModels::LevelSetContribution::Mathematics::op("advect",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("U"),Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")})),E::constantValue(0.0),{}});
    ExecutionScope advection;advection.kind=ExecutionKind::EquationCall;
    advection.step={"levelSetAdvection",{"phi"}};system.addExecution(std::move(advection));
    NumericalBinding transport{"levelSetAdvection","LevelSetAdvection"};
    transport.parameters={{"order",spec.advectionOrder},{"epsilon",spec.wenoEpsilon},{"power",spec.wenoPower}};
    transport.parameters.emplace("csf",spec.continuousSurfaceForce ? 1.0 : 0.0);
    transport.parameters.emplace("ghostFluid",spec.ghostFluid ? 1.0 : 0.0);
    transport.parameters.emplace("sigma",spec.surfaceTension);
    transport.parameters.emplace("width",spec.interfaceThickness);
    if (spec.continuousSurfaceForce) {
        auto force=E::multiply(E::multiply(E::constantValue(-spec.surfaceTension,Math::signature(0,{1,0,-2,0,0,0,0})),Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceCurvature")),
            E::multiply(Physics::InterfaceModels::LevelSetContribution::Mathematics::op("regularizedDelta",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi"),E::constantValue(spec.interfaceThickness,Math::signature(0,{0,1,0,0,0,0,0}))}),Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceNormal")));
        system.addEquation(System::Equation{"interfaceMomentumSource",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceForce"),force,{}});
        system.addEquation(System::Equation{"interfaceEnergySource",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfacePower"),E::op(Math::Op::Dot,{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("U"),force}),{}});
    }
    if (spec.ghostFluid)
        system.addEquation(System::Equation{"interfacePressureJump",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("pLiquidMinusGas"),
            E::multiply(E::constantValue(-spec.surfaceTension,Math::signature(0,{1,0,-2,0,0,0,0})),Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceCurvature")),{}});
    system.bindNumerics(std::move(transport));
    ExecutionScope tail;tail.kind=ExecutionKind::Sequence;tail.id="interface.finalize";
    if (spec.reinitializationSteps>0) {
        system.addEquation(System::Equation{"levelSetReference",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi0"),Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi"),{}});
        system.addEquation(System::Equation{"levelSetReinitialization",Physics::InterfaceModels::LevelSetContribution::Mathematics::op("ddtau",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")}),
            E::negate(E::multiply(Physics::InterfaceModels::LevelSetContribution::Mathematics::op("smoothedSign",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi0")}),
                E::subtract(Physics::InterfaceModels::LevelSetContribution::Mathematics::op("gradientNorm",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")}),E::constantValue(1.0)))),{}});
        ExecutionScope snapshot;snapshot.kind=ExecutionKind::EquationCall;
        snapshot.step={"levelSetReference",{"phi0",TargetKind::Workspace}};tail.children.push_back(snapshot);
        NumericalBinding reference{"levelSetReference","LevelSetReference"};
        reference.parameters={{"order",spec.reinitializationOrder},{"epsilon",spec.wenoEpsilon},
            {"power",spec.wenoPower},{"pseudoDt",spec.pseudoTimeStep},{"signFactor",spec.signSmoothingFactor}};
        system.bindNumerics(std::move(reference));
        ExecutionScope loop;loop.kind=ExecutionKind::Loop;loop.id="interface.pseudoTime";
        loop.repetitions=spec.reinitializationSteps;
        ExecutionScope stage;stage.kind=ExecutionKind::EquationCall;
        stage.step={"levelSetReinitialization",{"phi"}};loop.children.push_back(stage);tail.children.push_back(loop);
        NumericalBinding reinit{"levelSetReinitialization","LevelSetReinitialization"};
        reinit.inputs={"levelSet.phi0"};
        system.bindNumerics(std::move(reinit));
    }
    system.addEquation(System::Equation{"interfaceNormal",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceNormal"),
        Physics::InterfaceModels::LevelSetContribution::Mathematics::op("normalizedGradient",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")}),{}});
    system.addEquation(System::Equation{"interfaceGeometry",Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("interfaceCurvature"),
        Physics::InterfaceModels::LevelSetContribution::Mathematics::op("meanCurvature",{Physics::InterfaceModels::LevelSetContribution::Mathematics::symbol("phi")}),{}});
    ExecutionScope geometry;geometry.kind=ExecutionKind::EquationCall;
    geometry.step={"interfaceGeometry",{"interfaceCurvature"}};tail.children.push_back(geometry);
    system.placement.push_back({{},"equation:levelSetAdvection",{"equation:energy"},{"equation:interfaceGeometry"}});
    system.placement.push_back({{},"equation:interfaceGeometry",{"equation:levelSetAdvection"},{"kind:Commit"}});
    system.addExecution(std::move(tail));
    system.bindNumerics({"interfaceGeometry","LevelSetGeometry"});
    system.addClosure("normal = normalized gradient(phi); curvature and material properties from phi",{"interfaceNormal","interfaceGeometry"});
    if (spec.ghostFluid) system.addClosure("ghost-fluid interface pressure jump from phi/curvature",{"interfaceGeometry"});
}
}
