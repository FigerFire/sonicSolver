#include "SF_surfaceProjectionFixture.h"
#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        verifySurfaceProjection(parallel.rank(),parallel.size(),[&](std::vector<double>& values){runtime.globalSum(values);});
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
