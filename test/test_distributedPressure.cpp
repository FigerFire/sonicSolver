#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include "solver/linearAlgebra/SF_globalDofSystem.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        if (argc>1 && std::string(argv[1])=="--fail-on-one-rank") {
            if (parallel.rank()==1)
                throw std::runtime_error("injected failure before a peer collective");
            runtime.globalMaximum(1.0);
            return 2;
        }
        if (parallel.size()!=2) throw std::runtime_error("test requires two ranks");
        const double dt=runtime.globalMinimum(parallel.rank()==1?0.125:0.5);
        if (dt!=0.125) throw std::runtime_error("nonzero-rank dt restriction was lost");
        // Local ratios produce max(1/100,2/2)=1; the global ratio is 2/100.
        const double delta=runtime.globalMaximum(parallel.rank()==0?1.0:2.0);
        const double scale=runtime.globalMaximum(parallel.rank()==0?100.0:2.0);
        if (delta/scale!=0.02) throw std::runtime_error("global ratio reduced local ratios");
        const std::int64_t aboveDoubleExactness=(std::int64_t{1}<<53)+3;
        if (runtime.globalMinimum(aboveDoubleExactness+parallel.rank())
                !=aboveDoubleExactness
            || runtime.globalMaximum(aboveDoubleExactness+parallel.rank())
                !=aboveDoubleExactness+1
            || runtime.globalSum(aboveDoubleExactness+parallel.rank())
                !=2*aboveDoubleExactness+1) {
            throw std::runtime_error(
                "GlobalDof integer reductions lost identity above 2^53");
        }
        const auto range=runtime.allocateDistributedIndices(2);
        using namespace SF::LinearAlgebra;
        const auto id=[](std::int64_t row) {
            return GlobalDofId::make(GlobalDofSpace::Pressure,100+row);
        };
        std::vector<std::pair<GlobalDofId,std::int64_t>> entries;
        for (std::int64_t row=0;row<range.total;++row) entries.push_back({id(row),row});
        StaticDistributedNumbering numbering(range.first,range.last,range.total,entries);
        SF::FDM::LinearSolverConfig config;
        config.absoluteTolerance=1e-12;
        config.relativeTolerance=1e-10;
        DistributedLinearSystem solver(config,numbering);
        // Rank 0 has exactly zero RHS; rank 1 does not. No local early return.
        GlobalDofSystem system;
        for (std::int64_t row=range.first;row<=range.last;++row) {
            GlobalDofRow equation(id(row));
            equation.add(id(row),2.0);
            equation.setRightHandSide(row<2?0.0:2.0*(row+1));
            system.rows.push_back(equation);
        }
        auto result=solver.solve(system);
        if (!result.converged || result.solution.size()!=2)
            throw std::runtime_error("distributed nonzero solve failed");
        for (int local=0;local<2;++local) {
            const auto row=range.first+local;
            if (std::abs(result.solution[local]-(row<2?0.0:row+1.0))>1e-10)
                throw std::runtime_error("local zero RHS skipped distributed solve");
        }
        for (auto& row:system.rows) row.setRightHandSide(0.0);
        result=solver.solve(system);
        if (!result.converged || result.iterations!=0)
            throw std::runtime_error("global exact-zero RHS not recognized");
        for (auto& row:system.rows) row.setRightHandSide(parallel.rank()==1?1e-15:0.0);
        result=solver.solve(system);
        if (!result.converged || result.iterations!=0)
            throw std::runtime_error("global near-zero residual not checked against configured tolerance");
        if (parallel.isRoot()) std::cout << "distributed pressure primitives passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
