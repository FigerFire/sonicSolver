#include "infrastructure/mpi/SF_haloExchange.h"
#include "infrastructure/mpi/backend/SF_mpiBackend.h"
#include "infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.h"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc,char** argv) {
    try {
        SF::Parallel::Backend::MPIBackend backend(argc,argv,true);
        if (backend.size()!=2) throw std::runtime_error("requires two ranks");
        SF::MeshCommunication::HaloExchangePlan plan;
        plan.blockOwnerRanks={0,1};
        SF::Parallel::HaloExchange halo;
        halo.configure(&plan,backend.rank(),&backend);
        SF::Field field;
        field.setup(1,1,1,0,5);
        std::vector<double> values((size_t)field.TotalSize(),1.0);
        SF::Parallel::ScalarBlockValues view{backend.rank(),&field,&values,"fixture"};

        // A single rank's malformed scalar payload must reject both peers
        // before exchangeNeighbours and retain the origin rank in diagnostics.
        if (backend.rank()==1) values.pop_back();
        bool scalarRejected=false;
        try { halo.exchangeScalarValues({view}); }
        catch (const std::runtime_error& error) {
            const std::string message=error.what();
            scalarRejected=message.find("scalar preflight on rank ")!=std::string::npos
                && message.find(backend.rank()==1
                    ? "component count" : "peer rank")!=std::string::npos;
        }
        if (!backend.allRanksAgree(scalarRejected))
            throw std::runtime_error("rank-local scalar failure was not collective");

        values.assign((size_t)field.TotalSize(),1.0);
        bool mappingRejected=false;
        try {
            halo.exchangeScalarValues(backend.rank()==1
                ? std::vector<SF::Parallel::ScalarBlockValues>{view,view}
                : std::vector<SF::Parallel::ScalarBlockValues>{view});
        } catch (const std::runtime_error& error) {
            mappingRejected=std::string(error.what()).find("scalar preflight")
                !=std::string::npos;
        }
        if (!backend.allRanksAgree(mappingRejected))
            throw std::runtime_error("duplicate halo block mapping was not rejected");

        SF::MeshCommunication::HaloInterfaceFluxSyncGroup group;
        group.faces.emplace_back(); // invalid block id on both ranks
        plan.interfaceFluxSyncGroups.push_back(group);
        std::vector<SF::MeshBlockField> blocks(2);
        for (auto& block:blocks) block.field.setup(1,1,1,0,5);
        bool faceRejected=false;
        try {
            halo.assembleCanonicalInterfaceFluxes(blocks,{nullptr,nullptr},
                                                  {nullptr,nullptr});
        } catch (const std::runtime_error& error) {
            faceRejected=std::string(error.what()).find("canonical face preflight on rank ")
                !=std::string::npos;
        }
        if (!backend.allRanksAgree(faceRejected))
            throw std::runtime_error("canonical face mismatch was not collective");
        if (backend.rank()==0) std::cout<<"halo failure contract passed\n";
    } catch (const std::exception& error) {
        std::cerr<<"halo failure contract: "<<error.what()<<'\n';
        return 1;
    }
}
