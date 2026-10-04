#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        if (parallel.size()!=4) throw std::runtime_error("requires four ranks");
        // Every entity's owner is a consumer, but other consumers are sparse
        // and ordered differently. Exact COPY must preserve signed values.
        std::vector<std::int64_t> ids;
        for (std::int64_t id=11;id>=0;--id)
            if (runtime.ownsCanonicalEntity(id) || (id+parallel.rank())%3==0) ids.push_back(id);
        std::vector<double> values;
        for (auto id:ids) values.push_back(runtime.ownsCanonicalEntity(id)?-0.125*(id+1):0.0);
        runtime.copyCanonicalEntities(ids,values);
        for (std::size_t i=0;i<ids.size();++i)
            if (values[i]!=-0.125*(ids[i]+1)) throw std::runtime_error("canonical COPY differs from exact owner value");
        if (parallel.isRoot()) std::cout<<"sparse COPY passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
