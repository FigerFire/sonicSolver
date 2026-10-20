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
        system.addState(std::move(state));
    }
    system.addEquation(System::Equation{"levelSetAdvection",E::add(E::op("ddt",{E::symbol("phi")}),
        E::op("advect",{E::symbol("U"),E::symbol("phi")})),E::constantValue(0.0),{},true});
    ExecutionScope advection;advection.kind=ExecutionKind::EquationCall;
    advection.step={"levelSetAdvection",{"phi"}};system.addExecution(std::move(advection));
    NumericalBinding transport{"levelSetAdvection","LevelSetAdvection"};
    transport.parameters={{"order",spec.advectionOrder},{"epsilon",spec.wenoEpsilon},{"power",spec.wenoPower}};
    transport.parameters.emplace("csf",spec.continuousSurfaceForce ? 1.0 : 0.0);
    transport.parameters.emplace("ghostFluid",spec.ghostFluid ? 1.0 : 0.0);
    transport.parameters.emplace("sigma",spec.surfaceTension);
    transport.parameters.emplace("width",spec.interfaceThickness);
    if (spec.continuousSurfaceForce) {
        auto force=E::multiply(E::multiply(E::constantValue(-spec.surfaceTension),E::symbol("interfaceCurvature")),
            E::multiply(E::op("regularizedDelta",{E::symbol("phi"),E::constantValue(spec.interfaceThickness)}),E::symbol("interfaceNormal")));
        system.addEquation(System::Equation{"interfaceMomentumSource",E::symbol("interfaceForce"),force,{},true});
        system.addEquation(System::Equation{"interfaceEnergySource",E::symbol("interfacePower"),E::op("dot",{E::symbol("U"),force}),{},true});
    }
    if (spec.ghostFluid)
        system.addEquation(System::Equation{"interfacePressureJump",E::symbol("pLiquidMinusGas"),
            E::multiply(E::constantValue(-spec.surfaceTension),E::symbol("interfaceCurvature")),{},true});
    system.bindNumerics(std::move(transport));
    ExecutionScope tail;tail.kind=ExecutionKind::Sequence;tail.id="interface.finalize";
    if (spec.reinitializationSteps>0) {
        system.addEquation(System::Equation{"levelSetReference",E::symbol("phi0"),E::symbol("phi"),{},true});
        system.addEquation(System::Equation{"levelSetReinitialization",E::op("ddtau",{E::symbol("phi")}),
            E::negate(E::multiply(E::op("smoothedSign",{E::symbol("phi0")}),
                E::subtract(E::op("gradientNorm",{E::symbol("phi")}),E::constantValue(1.0)))),{},true});
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
    system.addEquation(System::Equation{"interfaceNormal",E::symbol("interfaceNormal"),
        E::op("normalizedGradient",{E::symbol("phi")}),{},true});
    system.addEquation(System::Equation{"interfaceGeometry",E::symbol("interfaceCurvature"),
        E::op("meanCurvature",{E::symbol("phi")}),{},true});
    ExecutionScope geometry;geometry.kind=ExecutionKind::EquationCall;
    geometry.step={"interfaceGeometry",{"interfaceCurvature"}};tail.children.push_back(geometry);
    system.placement.push_back({{},"equation:levelSetAdvection",{"equation:energy"},{"equation:interfaceGeometry"}});
    system.placement.push_back({{},"equation:interfaceGeometry",{"equation:levelSetAdvection"},{"kind:Commit"}});
    system.addExecution(std::move(tail));
    system.bindNumerics({"interfaceGeometry","LevelSetGeometry"});
    system.addClosure("normal = normalized gradient(phi); curvature and material properties from phi");
    if (spec.ghostFluid) system.addClosure("ghost-fluid interface pressure jump from phi/curvature");
}
}
