#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include "solver/linearAlgebra/SF_globalDofSystem.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        if(parallel.size()!=2) throw std::runtime_error("requires two ranks");
        SF::Execution::Runtime runtime(&parallel.coordinator());
        using namespace SF::LinearAlgebra;
        // One marker spans both velocity owners: [2 0 .5; 0 3 .5; .5 .5 0].
        const auto u0=GlobalDofId::make(GlobalDofSpace::Velocity,0);
        const auto u1=GlobalDofId::make(GlobalDofSpace::Velocity,1);
        const auto lambda=GlobalDofId::make(GlobalDofSpace::Constraint,0);
        const int rank=parallel.rank();
        StaticDistributedNumbering numbering(rank==0?0:1,rank==0?0:2,3,{{u0,0},{u1,1},{lambda,2}});
        SF::FDM::LinearSolverConfig config;config.relativeTolerance=1e-11;config.absoluteTolerance=1e-13;
        GlobalDofSystem equations;
        GlobalDofRow primal(rank==0?u0:u1);primal.add(rank==0?u0:u1,rank==0?2:3);primal.add(lambda,.5);
        primal.setRightHandSide(rank==0?3.5:7.5);equations.rows.push_back(primal);
        if(rank==1) {GlobalDofRow constraint(lambda);constraint.add(u0,.5);constraint.add(u1,.5);constraint.setRightHandSide(1.5);equations.rows.push_back(constraint);}
        DistributedLinearSystem solver(config,numbering);auto solved=solver.solve(equations);
        if(!solved.converged || solved.solution.size()!=equations.rows.size()) throw std::runtime_error("KKT solve did not converge");
        if(std::abs(solved.solution[0]-(rank==0?1:2))>1e-9 || (rank==1 && std::abs(solved.solution[1]-3)>1e-9)) throw std::runtime_error("KKT primal/multiplier differs from exact constrained solution");
        std::vector<std::int64_t> ids{0};std::vector<double> value{runtime.ownsCanonicalEntity(0)?3.:0.};runtime.copyCanonicalEntities(ids,value);
        if(value[0]!=3) throw std::runtime_error("lambda COPY failed");
        if(parallel.isRoot()) std::cout<<"minimal KKT passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
