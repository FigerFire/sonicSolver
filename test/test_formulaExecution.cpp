#include "core/system/SF_formula.h"
#include "solver/system/SF_formulaCompiler.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_formulaStorage.h"
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_equationContribution.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "solver/run/SF_planExecutor.h"
#include "infrastructure/mpi/SF_parallelContext.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace SF;
using System::FormulaExpr;

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

System::Equation diffusion(std::string id,std::string target) {
    auto gradient=FormulaExpr::op("grad",{FormulaExpr::symbol(target)});
    auto flux=FormulaExpr::multiply(FormulaExpr::symbol("nu"),
                                    std::move(gradient));
    return {std::move(id),
        FormulaExpr::negate(FormulaExpr::op("div",{std::move(flux)},"diffusion")),
        FormulaExpr::symbol("S"),{}};
}

System::FormulaValues valuesFor(std::string target,int components,
                                std::vector<double>& field) {
    System::FormulaValues values;
    values.ownedCells=9;
    values.components=components;
    values.read=[&field,target,components](const std::string& symbol,int cell,int c) {
        if (symbol=="nu") return 1.0;
        if (symbol=="phi") return 1.0;
        if (symbol=="S") return 0.0;
        if (symbol!=target) throw std::runtime_error("Unbound Equation symbol: "+symbol);
        return field.at(static_cast<std::size_t>(cell)*components+c);
    };
    values.write=[&field,target,components](const std::string& symbol,int cell,int c,double value) {
        if (symbol!=target) throw std::runtime_error("Wrong Equation write target.");
        field.at(static_cast<std::size_t>(cell)*components+c)=value;
    };
    values.dof=[target](const std::string& symbol,int cell,int c) {
        if (symbol!=target) throw std::runtime_error("Wrong Equation DOF target.");
        return LinearAlgebra::GlobalDofId::make(
            LinearAlgebra::GlobalDofSpace::ScalarTransport,cell,c);
    };
    values.boundaryValue=[](const std::string&,int,int,int,int c) {
        return static_cast<double>(c+1);
    };
    return values;
}

LinearAlgebra::SolveResult solve(const LinearAlgebra::GlobalDofSystem& system,
                                int components) {
    using namespace LinearAlgebra;
    const int count=9*components;
    std::vector<std::pair<GlobalDofId,std::int64_t>> entries;
    for (int cell=0;cell<9;++cell)
        for (int c=0;c<components;++c)
            entries.push_back({GlobalDofId::make(
                GlobalDofSpace::ScalarTransport,cell,c),cell*components+c});
    StaticDistributedNumbering numbering(0,count-1,count,entries);
    FDM::LinearSolverConfig config;
    config.absoluteTolerance=1e-12;
    config.relativeTolerance=1e-11;
    DistributedLinearSystem solver(config,numbering);
    return solver.solve(system);
}

} // namespace

int main(int argc,char** argv) {
    try {
        Parallel::ParallelContext parallel(argc,argv,true);
        System::EquationRegistry formulas;
        formulas.add(diffusion("anyName","Y_O2"));
        System::FormulaOperatorCatalog catalog;
        catalog.add(Discretization::central2FormulaDiffusion({3,3,1.0,1.0}));
        auto alternate=Discretization::central2FormulaDiffusion({3,3,1.0,1.0});
        alternate.id="central2.cartesian.alternate";
        catalog.add(std::move(alternate));
        const std::vector<System::FormulaOperatorBinding> bindings{
            {"","div","central2.cartesian.alternate"},
            {"anyName","","central2.cartesian.alternate"},
            {"anyName","diffusion","central2.cartesian.diffusion"}};
        System::ExecutableEquationSystem executable;
        System::StateSymbol oxygen;
        oxygen.id="Y_O2";
        oxygen.storageKey="scalarStorage";
        executable.state.add(oxygen);
        executable.registry=formulas;
        const auto lowered=System::compileExecutionProgram(executable,
            System::ExecutionProgram{System::EquationCall{"anyName",{"Y_O2"}}},
            {{"anyName","LinearEquation"}},System::builtinProviders());
        const auto compiled=System::FormulaCompiler::compile(formulas,
            {lowered.steps.front().calls.front().equation,
             lowered.steps.front().calls.front().target,
             System::Legacy::FormulaMode::Implicit},catalog,bindings);
        require(compiled.backend=="formula.linear-assembly",
                "Implicit Equation did not compile to linear assembly.");
        require(compiled.numericalProviders.size()==1
                    && compiled.numericalProviders.front()
                        =="central2.cartesian.diffusion",
                "Occurrence binding did not override Equation/global defaults.");
        bool badPathRejected=false;
        try {
            (void)System::FormulaCompiler::compile(formulas,
                {"anyName","Y_O2",System::Legacy::FormulaMode::Implicit},catalog,
                {{"anyName","difffusion","central2.cartesian.diffusion"}});
        } catch (const std::runtime_error& error) {
            badPathRejected=std::string(error.what()).find("does not address")
                !=std::string::npos;
        }
        require(badPathRejected,"Unknown numerical operator path was ignored.");
        bool duplicatePathRejected=false;
        try {
            System::EquationRegistry invalid;
            invalid.add({"duplicatePath",
                FormulaExpr::add(
                    FormulaExpr::op("grad",{FormulaExpr::symbol("C")},"same"),
                    FormulaExpr::op("div",{FormulaExpr::symbol("C")},"same")),
                FormulaExpr::constantValue(0.0),{}});
        } catch (const std::runtime_error& error) {
            duplicatePathRejected=std::string(error.what()).find(
                "duplicate operator occurrence")!=std::string::npos;
        }
        require(duplicatePathRejected,"Duplicate operator paths were accepted.");
        std::vector<double> scalar(9,0.0);
        auto scalarValues=valuesFor("Y_O2",1,scalar);
        auto matrix=compiled.assemble(scalarValues);
        require(matrix.rows.size()==9,"Scalar Equation row count is wrong.");
        require(matrix.rows.front().coefficients().size()==3,
                "Central2 corner stencil is wrong.");
        System::CompiledSolvePlan implicitPlan;
        implicitPlan.root.kind=System::PlanNodeKind::Sequence;
        implicitPlan.root.id="formula.implicit.test";
        implicitPlan.root.children.push_back({});
        implicitPlan.root.children.back().kind=System::PlanNodeKind::Solve;
        implicitPlan.root.children.back().id="formula.implicit";
        implicitPlan.root.children.back().operation="formula.linear.solve";
        implicitPlan.root.children.back().equationCalls.push_back(
            {compiled.call.equation,compiled.call.target});
        Run::OpRegistry implicitOperations;
        implicitOperations.bind("formula.linear.solve",[&] {
            compiled.writeSolution(scalarValues,solve(compiled.assemble(scalarValues),1));
        });
        Run::PlanExecutor::execute(implicitPlan,implicitOperations);
        for (double value:scalar)
            require(std::abs(value-1.0)<1e-9,
                    "Implicit arbitrary scalar Equation solution is wrong.");

        formulas.add(diffusion("vectorFormula","V"));
        std::vector<double> vector(18,0.0);
        auto vectorValues=valuesFor("V",2,vector);
        const auto vectorCall=System::FormulaCompiler::compile(formulas,
            {"vectorFormula","V",System::Legacy::FormulaMode::Implicit},catalog,
            {{"","div","central2.cartesian.alternate"},
             {"vectorFormula","","central2.cartesian.diffusion"}});
        require(vectorCall.numericalProviders.size()==1
                    && vectorCall.numericalProviders.front()
                        =="central2.cartesian.diffusion",
                "Equation binding did not override global default.");
        vectorCall.writeSolution(vectorValues,solve(vectorCall.assemble(vectorValues),2));
        for (int cell=0;cell<9;++cell)
            for (int c=0;c<2;++c)
                require(std::abs(vector[cell*2+c]-(c+1.0))<1e-9,
                        "Implicit vector Equation solution is wrong.");

        // Use actual registered Field storage and StateRealization, rather
        // than a test-owned vector, for the same scalar/vector method path.
        Field physical;
        physical.setup(3,3,1,1,3);
        State::StateBundle physicalState;
        physicalState.patches={&physical};
        physicalState.registerConservativeState();
        System::ExecutableEquationSystem storedSystem;
        System::StateSymbol storedScalar;
        storedScalar.id="Y_O2";
        storedScalar.storageBinding=System::StorageBinding::PackedDistributed;
        storedScalar.storageKey="conservative";
        storedScalar.componentOffset=0;
        storedSystem.state.add(storedScalar);
        System::StateSymbol storedVector=storedScalar;
        storedVector.id="V";
        storedVector.components=2;
        storedVector.componentOffset=1;
        storedVector.shape=System::ValueShape::Vector;
        storedSystem.state.add(storedVector);
        System::StateSymbol viscosity;
        viscosity.id="nu";
        viscosity.runtimeStorageRequired=false;
        viscosity.constantValue=1.0;
        storedSystem.state.add(viscosity);
        System::StateSymbol source=viscosity;
        source.id="S";
        source.constantValue=0.0;
        storedSystem.state.add(source);
        const auto realized=System::realizeState(
            storedSystem.state,System::RuntimeRequirements{},physicalState);
        const auto boundary=[](const std::string&,int,int,int,int component) {
            return static_cast<double>(component+1);
        };
        auto boundScalar=System::bindFormulaCellValues(
            realized,{"Y_O2"},boundary);
        compiled.writeSolution(boundScalar,solve(compiled.assemble(boundScalar),1));
        auto boundVector=System::bindFormulaCellValues(
            realized,{"V"},boundary);
        vectorCall.writeSolution(boundVector,solve(vectorCall.assemble(boundVector),2));
        for (int j=1;j<=3;++j)
            for (int i=1;i<=3;++i) {
                require(std::abs(physical(i,j,1,0)-1.0)<1e-9,
                        "Registered scalar Field binding did not receive linear solution.");
                require(std::abs(physical(i,j,1,1)-1.0)<1e-9
                            && std::abs(physical(i,j,1,2)-2.0)<1e-9,
                        "Registered vector Field binding did not receive linear solution.");
            }
        formulas.add({"storedAssignment",FormulaExpr::symbol("Y_O2"),
            FormulaExpr::add(FormulaExpr::symbol("Y_O2"),
                             FormulaExpr::constantValue(2.0)),{}});
        storedSystem.registry=formulas;
        const auto storedMethod=System::compileExecutionProgram(storedSystem,
            System::ExecutionProgram{System::EquationCall{"storedAssignment",{"Y_O2"}}},
            {{"storedAssignment","DirectEvaluation"}},
            System::builtinProviders());
        const auto storedAssign=System::FormulaCompiler::compile(formulas,
            {storedMethod.steps.front().calls.front().equation,
             storedMethod.steps.front().calls.front().target,
             System::Legacy::FormulaMode::Assign},catalog,{});
        storedAssign.assign(boundScalar);
        for (int j=1;j<=3;++j)
            for (int i=1;i<=3;++i)
                require(std::abs(physical(i,j,1,0)-3.0)<1e-9,
                        "Registered scalar Field binding did not receive direct assignment.");

        // The same equation writes lazy Working/Correction views. Neither
        // registration nor a view write may publish into physical state.
        System::ExecutionProgram temporaryProgram{
            System::EquationCall{"storedAssignment",System::targetFromSyntax("Y_O2*")},
            System::EquationCall{"storedAssignment",System::targetFromSyntax("Y_O2'")}};
        const auto temporary=System::compileExecutionProgram(storedSystem,temporaryProgram,
            {{"storedAssignment","DirectEvaluation"}},System::builtinProviders());
        require(temporary.stateViews.size()==2 && temporaryProgram.root.children.size()==2,
                "STATE realization changed HOW or created unrequested versions.");
        const auto physicalBeforeTemporary=physical(1,1,1,0);
        const auto temporaryState=System::realizeState(storedSystem.state,
            System::RuntimeRequirements{},physicalState,temporary.stateViews);
        for (const auto& method:temporary.steps) {
            auto values=System::bindFormulaCellValues(temporaryState,method.target,boundary);
            storedAssign.assign(values);
        }
        for (int j=1;j<=3;++j)
            for (int i=1;i<=3;++i) {
                const int cell=physical.getIdx(i,j,1);
                require(std::abs(physical(i,j,1,0)-3.0)<1e-9
                    && std::abs(temporaryState.view("Y_O2",System::StateViewKind::Working).fields.front()->read(cell,0)-5.0)<1e-9
                    && temporaryState.view("Y_O2",System::StateViewKind::Correction).fields.front()->read(cell,0)==2.0,
                    "Temporary assignment corrupted the physical STATE authority.");
            }

        require(physical(1,1,1,0)==physicalBeforeTemporary,
                "Working/Correction writes changed the physical state bitwise.");

        // A standard symbol is available without any user addState call.
        System::ExecutableEquationSystem common;
        const System::BuiltinStateCatalog builtinCatalog;
        require(builtinCatalog.contains("U") && builtinCatalog.contains("k")
                    && common.state.size()==0,
                "Available builtin metadata activated case STATE.");
        System::requireTargetStates(common.state,
            System::ExecutionProgram{System::EquationCall{"pressureDelta",System::targetFromSyntax("p'")}}.root);
        common.registry.add({"pressureDelta",FormulaExpr::symbol("p"),
            FormulaExpr::constantValue(7.0),{}});
        System::ExecutableEquationSystem mathematicsOnly;
        mathematicsOnly.registry=common.registry;
        const auto commonMethod=System::compileExecutionProgram(mathematicsOnly,common.state,
            System::ExecutionProgram{System::EquationCall{"pressureDelta",System::targetFromSyntax("p'")}},
            {{"pressureDelta","DirectEvaluation"}},System::builtinProviders());
        require(common.state.contains("p") && common.state.size()==1
                    && !common.state.contains("U") && !common.state.contains("k")
                    && !common.state.contains("omega") && !common.state.contains("rhoE")
                    && !common.state.contains("p'") && !common.state.contains("pPrime")
                    && commonMethod.stateViews.size()==1,
                "Builtin STATE catalog requires user registration or synthetic symbols.");
        const auto commonState=System::realizeState(common.state,System::RuntimeRequirements{},
            physicalState,commonMethod.stateViews);
        auto commonValues=System::bindFormulaCellValues(commonState,commonMethod.steps.front().target,boundary);
        const auto commonAssign=System::FormulaCompiler::compile(common.registry,
            {"pressureDelta","p",System::Legacy::FormulaMode::Assign},catalog,{});
        commonAssign.assign(commonValues);
        require(commonState.view("p",System::StateViewKind::Correction).fields.front()->read(
                    physical.getIdx(1,1,1),0)==7.0 && std::abs(physical(1,1,1,0)-3.0)<1e-9,
                "Builtin correction view was not realized independently of physical state.");
        bool missingBaseRejected=false;
        try {
            System::compileExecutionProgram(common,
                System::ExecutionProgram{System::EquationCall{"pressureDelta",{"missing",System::TargetKind::Working}}},
                {{"pressureDelta","DirectEvaluation"}},System::builtinProviders());
        } catch (const std::runtime_error& error) {
            missingBaseRejected=std::string(error.what()).find("base STATE symbol")!=std::string::npos;
        }
        require(missingBaseRejected,"Working target bypassed STATE base lookup.");

        formulas.add({"assignment",FormulaExpr::symbol("C"),
            FormulaExpr::add(FormulaExpr::symbol("C"),
                             FormulaExpr::constantValue(2.0)),{}});
        std::vector<double> assignment(9,1.0);
        auto assignmentValues=valuesFor("C",1,assignment);
        const auto assign=System::FormulaCompiler::compile(formulas,
            {"assignment",{},System::Legacy::FormulaMode::Assign},catalog,{});
        System::CompiledSolvePlan plan;
        plan.root.kind=System::PlanNodeKind::Sequence;
        plan.root.id="formula.test";
        plan.root.children.push_back({});
        plan.root.children.back().kind=System::PlanNodeKind::Update;
        plan.root.children.back().id="formula.assignment";
        plan.root.children.back().operation="formula.assignment";
        plan.root.children.back().equationCalls.push_back(
            {assign.call.equation,assign.call.target});
        Run::OpRegistry operations;
        operations.bind("formula.assignment",[&] { assign.assign(assignmentValues); });
        Run::PlanExecutor::execute(plan,operations);
        for (double value:assignment)
            require(value==3.0,"Equation assignment or PlanExecutor dispatch failed.");

        formulas.add({"explicitTransport",
            FormulaExpr::add(
                FormulaExpr::op("ddt",{FormulaExpr::symbol("C")}),
                diffusion("unused","C").lhs),
            FormulaExpr::constantValue(0.0),{}});
        const auto explicitCall=System::FormulaCompiler::compile(formulas,
            {"explicitTransport","C",System::Legacy::FormulaMode::Explicit},catalog,
            {{"","div","central2.cartesian.diffusion"}});
        require(explicitCall.backend=="formula.explicit-residual",
                "Explicit Equation did not compile a residual.");
        require(std::abs(explicitCall.evaluateRhs(assignmentValues,0,0)+4.0)<1e-12,
                "Explicit Central2 Equation residual is wrong at a boundary cell.");
        for (double value:assignment) require(value==3.0,"Explicit compile changed state.");

        formulas.add({"scalarTransport",
            FormulaExpr::add(
                FormulaExpr::add(
                    FormulaExpr::op("ddt",{FormulaExpr::symbol("C")}),
                    FormulaExpr::op("div",{FormulaExpr::symbol("phi"),
                                            FormulaExpr::symbol("C")},"convection")),
                diffusion("unused","C").lhs),
            FormulaExpr::symbol("S"),{}});
        System::FormulaOperatorCatalog transportCatalog;
        transportCatalog.add(
            Discretization::central2FormulaDiffusion({3,3,1.0,1.0}));
        transportCatalog.add({"test.upwind1.scalar","div",
            [](const FormulaExpr& expression) {
                return expression.arguments.size()==2
                    && expression.arguments[0].kind==FormulaExpr::Kind::Symbol
                    && expression.arguments[1].kind==FormulaExpr::Kind::Symbol;
            },
            [](const FormulaExpr& expression) {
                const auto flux=expression.arguments[0].name;
                const auto target=expression.arguments[1].name;
                return System::FormulaValueKernel{
                    [flux,target](const System::FormulaValues& values,int cell,int c) {
                        const double faceFlux=values.read(flux,cell,0);
                        if (faceFlux<0.0)
                            throw std::runtime_error("Test upwind provider needs positive flux.");
                        const double upstream=cell%3==0
                            ? values.boundaryValue(target,cell,0,-1,c)
                            : values.read(target,cell-1,c);
                        return faceFlux*(values.read(target,cell,c)-upstream);
                    }};
            },{}});
        const auto transport=System::FormulaCompiler::compile(formulas,
            {"scalarTransport","C",System::Legacy::FormulaMode::Explicit},
            transportCatalog,{});
        require(transport.numericalProviders.size()==2
                    && std::abs(transport.evaluateRhs(assignmentValues,0,0)+6.0)
                        <1e-12,
                "Generic scalar transport did not combine convection and diffusion.");

        System::FormulaOperatorCatalog weno;
        weno.add({"weno.evaluate","div",
            [](const FormulaExpr&) { return true; },
            [](const FormulaExpr&) { return System::FormulaValueKernel{
                [](const System::FormulaValues&,int,int) { return 0.0; }}; },{}});
        formulas.add({"unsupportedWeno",
            FormulaExpr::op("div",{FormulaExpr::symbol("C")},"convection"),
            FormulaExpr::constantValue(0.0),{}});
        bool rejected=false;
        try {
            (void)System::FormulaCompiler::compile(formulas,
                {"unsupportedWeno","C",System::Legacy::FormulaMode::Implicit},weno,{});
        } catch (const std::runtime_error& error) {
            rejected=std::string(error.what()).find("linear assembly capability")
                !=std::string::npos;
        }
        require(rejected,"Implicit WENO silently became an explicit RHS.");

        const auto old=SF::Equation::named("manual",
            SF::Equation::ddt({"C"})+SF::Equation::div({"flux"})
                == SF::Equation::Symbol{"zero"});
        const auto builtIn=System::formulaFromEquation(old,
            {System::OriginKind::BuiltinPreset,"builtin"});
        System::Equation manual{"myOwnFormula",
            FormulaExpr::add(
                FormulaExpr::op("ddt",{FormulaExpr::symbol("C")}),
                FormulaExpr::op("div",{FormulaExpr::symbol("flux")})),
            FormulaExpr::constantValue(0.0),
            {System::OriginKind::User,"manual"}};
        require(System::canonicalFormula(builtIn)
                    ==System::canonicalFormula(manual),
                "Built-in and manual Equation ASTs differ by provenance.");
        System::RawEquationSystem composed;
        std::vector<System::TransformationDescriptor> transformations;
        std::vector<System::LegacyExecutionPolicy> policies;
        System::SystemContribution builtinContribution;
        builtinContribution.addEquation(builtIn);
        System::SystemCompositionBuilder builtinBuilder(composed,transformations,
            policies,{System::OriginKind::BuiltinPreset,"builtin"});
        builtinBuilder.applyContribution(std::move(builtinContribution));
        System::SystemContribution userContribution;
        userContribution.addEquation(manual);
        System::SystemCompositionBuilder userBuilder(composed,transformations,
            policies,{System::OriginKind::User,"manual"});
        userBuilder.applyContribution(std::move(userContribution));
        System::FormulaOperatorCatalog convection;
        convection.add({"test.flux.div","div",
            [](const FormulaExpr& expression) {
                return expression.arguments.size()==1
                    && expression.arguments[0].kind==FormulaExpr::Kind::Symbol
                    && expression.arguments[0].name=="flux";
            },
            [](const FormulaExpr&) { return System::FormulaValueKernel{
                [](const System::FormulaValues&,int,int) { return 2.0; }}; },{}});
        const auto presetCall=System::FormulaCompiler::compile(composed.registry,
            {"manual","C",System::Legacy::FormulaMode::Explicit},convection,{});
        const auto userCall=System::FormulaCompiler::compile(composed.registry,
            {"myOwnFormula","C",System::Legacy::FormulaMode::Explicit},convection,{});
        require(presetCall.backend==userCall.backend
                    && presetCall.numericalProviders==userCall.numericalProviders
                    && presetCall.evaluateRhs(assignmentValues,0,0)
                        ==userCall.evaluateRhs(assignmentValues,0,0),
                "Built-in and user Equation did not compile through the same backend.");
        std::cout << "Equation execution contract passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
