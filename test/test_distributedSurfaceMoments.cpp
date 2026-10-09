#include "SF_surfaceMomentFixture.h"
#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        verifySurfaceMoments(parallel.rank(),parallel.size(),[&](std::vector<double>& data){runtime.globalSum(data);});
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
