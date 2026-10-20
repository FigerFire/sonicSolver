#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include "solver/linearAlgebra/SF_distributedConstraintLayout.h"
#include <set>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        std::vector<SF::GlobalConstraintDofId> ids;
        for (int i=0;i<13;++i) ids.push_back(SF::GlobalConstraintDofId::fromMarker(SF::GlobalMarkerId::fromSurfacePrimitive(i)));
        auto layout=SF::LinearAlgebra::DistributedConstraintLayout::build(ids,runtime);
        std::set<std::int64_t> rows;
        for(auto id:ids) {
            const auto row=layout.row(id);rows.insert(row);
            if(runtime.globalMinimum(row)!=runtime.globalMaximum(row)) throw std::runtime_error("replica row identity diverged");
            if(runtime.globalSum(std::int64_t(runtime.ownsCanonicalEntity(id.value())))!=1) throw std::runtime_error("constraint has multiple owners");
        }
        if(layout.globalSize()!=13 || rows.size()!=13 || *rows.begin()!=0 || *rows.rbegin()!=12) throw std::runtime_error("constraint rows not unique and contiguous");
        if(parallel.isRoot()) std::cout<<"constraint layout passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
