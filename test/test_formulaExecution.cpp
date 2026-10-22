#include "test/typedSymbols.h"
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
    auto gradient=FormulaExpr::op("grad",{SF::Test::symbol(target)});
    auto flux=FormulaExpr::multiply(SF::Test::symbol("nu"),
                                    std::move(gradient));
    return {std::move(id),
        FormulaExpr::negate(FormulaExpr::op("div",{std::move(flux)},"diffusion")),
        SF::Test::symbol("S"),{}};
}

System::FormulaValues valuesFor(std::string target,int components,
                                std::vector<double>& field) {
    System::FormulaValues values;
    values.ownedCells=9;
    values.components=components;
    values.read=[&field,target,components](const std::string& symbol,int cell,int c) {
        if (symbol=="nu") return 1.0;
        if (symbol=="flux") return 1.0;
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
        auto equation=diffusion("arbitraryEquation","Y_O2");
        System::FormulaOperatorCatalog catalog;
        auto provider=Discretization::central2FormulaDiffusion({3,3,1.,1.});catalog.add(provider);
        const auto& occurrence=equation.lhs.arguments()[0];
        const auto& selected=catalog.resolve(equation,occurrence,{{equation.id,"diffusion",provider.id}});
        const auto linear=selected.compileLinear(occurrence,"Y_O2");
        std::vector<double> field(9,0);auto values=valuesFor("Y_O2",1,field);
        LinearAlgebra::GlobalDofSystem matrix;
        for(int cell=0;cell<9;++cell) {
            const auto term=linear(values,cell,0);
            LinearAlgebra::GlobalDofRow row(values.dof("Y_O2",cell,0));
            for(const auto& coefficient:term.coefficients)row.add(coefficient.column,-coefficient.value);
            row.setRightHandSide(term.constant);row.setInitialGuess(0);matrix.rows.push_back(row);
        }
        const auto result=solve(matrix,1);require(result.converged,"Existing Central2/HYPRE linear solve failed");
        for(double x:result.solution)require(std::abs(x-1)<1e-9,"Central2 fixed boundary solution changed");
        const auto explicitKernel=selected.compileValue(occurrence);
        require(explicitKernel(values,0,0)==2 && explicitKernel(values,4,0)==0,"Existing Central2 value kernel changed");
        bool rejected=false;
        try{catalog.resolve(equation,occurrence,{{equation.id,"diffusion","missing"}});}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Explicit incompatible provider silently replaced");
        rejected=false;try{catalog.resolve(equation,occurrence,{{equation.id,"missingOccurrence",provider.id}});}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Unmatched explicit occurrence was silently ignored");
        auto c=SF::Test::symbol("C");auto k=FormulaExpr::constantValue(.7,Math::signature(0,{0,0,-1,0,0,0,0}));
        auto rhs=FormulaExpr::op(Math::Op::Diffusion,{FormulaExpr::constantValue(.03,Math::signature(0,{0,2,-1,0,0,0,0})),c})-k*FormulaExpr::power(c,2);
        System::Equation nonlinear{"reaction",FormulaExpr::op(Math::Op::Ddt,{c}),rhs};
        require(System::inspectFormula(rhs,"C").targetDegree==2,"Nonlinear structure not extracted");
        require(System::formulaText({"",System::explicitFormulaRhs(nonlinear,"C"),rhs})==System::formulaText({"",rhs,rhs}),"Bare ddt reordered RHS");
        rejected=false;try{System::Equation wrong{"wrongK",nonlinear.lhs,FormulaExpr::op(Math::Op::Diffusion,{FormulaExpr::constantValue(.03,Math::signature(0,{0,2,-1,0,0,0,0})),c})-FormulaExpr::constantValue(.7)*FormulaExpr::power(c,2)};}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Wrong reaction coefficient dimensions bypassed construction");
        for(const auto& invalid:std::vector<System::Equation>{{"algebraic",c,c},{"massMatrix",FormulaExpr::constantValue(2)*nonlinear.lhs,rhs}}) {
            rejected=false;try{System::explicitFormulaRhs(invalid,"C");}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"Unsupported temporal structure silently inverted or converted");
        }
        System::ExecutableEquationSystem system;System::StateSymbol state;state.id=state.storageKey="C";state.signature=c.signature();state.storageBinding=System::StorageBinding::NamedDistributed;system.state.add(state);system.registry.add(nonlinear);
        rejected=false;try{System::builtinProviders().at("LinearEquation").compile(system,{"reaction",{"C"}},{"reaction","LinearEquation"});}catch(const std::runtime_error& e){rejected=std::string(e.what()).find("Unsupported")!=std::string::npos;}
        require(rejected,"Nonlinear explicit equation masqueraded as implemented implicit solve");
        system.registry.add({"spoof",FormulaExpr::symbol("C",SF::Test::symbol("p").signature()),FormulaExpr::symbol("C",SF::Test::symbol("p").signature())});
        rejected=false;try{System::compileExecutionProgram(system,System::ExecutionProgram{System::EquationCall{"spoof",{"C"}}},{{"spoof","DirectEvaluation"}},System::builtinProviders());}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Same-name logical symbol overrode STATE signature");
        auto other=state;other.id=other.storageKey="partner";system.state.add(other);
        auto otherDiffusion=FormulaExpr::op(Math::Op::Diffusion,{FormulaExpr::constantValue(.03,Math::signature(0,{0,2,-1,0,0,0,0})),other.expression()});
        system.registry.add({"wrongDiffusion",nonlinear.lhs,rhs+otherDiffusion});
        rejected=false;try{System::builtinProviders().at("ScalarDiffusionCentral2").compile(system,{"wrongDiffusion",{"C"}},{"wrongDiffusion","ScalarDiffusionCentral2"});}catch(const std::runtime_error& e){rejected=std::string(e.what()).find("only diffuses its declared target")!=std::string::npos;}
        require(rejected,"Selected diffusion kernel silently consumed another variable's operator");
        const Math::OperatorRule twice{"test.module.twice",[](const std::vector<Math::Signature>& a){if(a.size()!=1 || a[0].rank)throw std::runtime_error("scalar required");return a[0];}};
        const Math::OperatorRule conflict{"test.module.twice",[](const std::vector<Math::Signature>& a){return a.at(0);}};
        auto extension=FormulaExpr::op(twice,{c},"module.operator");System::Equation ext{"module",c,extension};
        System::EquationRegistry definitions;definitions.add(ext);
        rejected=false;try{definitions.add({"conflict",c,FormulaExpr::op(conflict,{c})});}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Conflicting mathematical extension rules shared an identity");
        System::FormulaOperatorProvider implementation;implementation.id="test.module.twice.values";implementation.mathematicalOperator=Math::Op::Extension;implementation.extensionName=twice.id;implementation.extensionRule=twice.result;
        implementation.matches=[](const auto& n){return n.arguments().size()==1 && n.arguments()[0].kind()==FormulaExpr::Kind::Symbol;};
        implementation.compileValue=[](const auto&){return System::FormulaValueKernel{[](const System::FormulaValues& v,int cell,int component){return 2*v.boundReads.at(0)(cell,component);}};};
        System::FormulaOperatorCatalog module;module.add(implementation);std::vector<std::string> used;
        const auto kernel=System::compileFormulaValue(ext,extension,module,used,{"C"});System::FormulaValues input;input.boundReads.push_back([](int,int){return 3.;});
        require(kernel(input,0,0)==6 && used==std::vector<std::string>{implementation.id},"Module extension did not bind its existing catalog implementation");
        implementation.id="test.module.conflict";implementation.extensionRule=conflict.result;
        rejected=false;try{module.add(implementation);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Numerics accepted a second mathematical rule for the same extension");
        std::cout<<"Formula structure, module extension and independent Central2 linear/value verification passed\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
