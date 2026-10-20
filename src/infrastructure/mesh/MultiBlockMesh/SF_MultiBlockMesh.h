/// @file SF_MultiBlockMesh.h
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include "core/field/SF_field.h"
#include "SF_cellFaceMesh.h"
#include "SF_communicationPlan.h"
#include "SF_meshConfig.h"
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace SF {

namespace Physics {
namespace Multiphase {
struct MultiPhaseConfig;
} // namespace Multiphase
} // namespace Physics

/// @brief 分解前附着在源zone物理点上的IBM几何分类。
struct RawIBMPointData {
    int cellType = 0;
    int ghostLayer = 0;
    double signedDistance = 0.0;
    bool hasGeometry = false;
    std::array<double, 3> wallPoint{0.0, 0.0, 0.0};
    std::array<double, 3> imagePoint{0.0, 0.0, 0.0};
    std::array<double, 3> wallNormal{0.0, 0.0, 0.0};
    std::array<double, 3> wallVelocity{0.0, 0.0, 0.0};
};

struct RawMeshBlock {
    std::string name;
    std::string sourceFile;
    int sourceZoneId = -1;
    int ownerRank = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int ng = 0;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<int> globalPointIds;
    std::vector<RawIBMPointData> ibmPointData;
    std::map<std::string, std::vector<int>> pointSets;
    CellFaceMesh topology;
};

struct MeshBlockField {
    std::string name;
    std::string sourceFile;
    int sourceZoneId = -1;
    int ownerRank = 0;
    std::vector<int> globalPointIds;
    std::vector<RawIBMPointData> ibmPointData;
    CellFaceMesh topology;
    Field field;
};

/// @brief 一个MPI rank拥有的复合网格partition。
struct MeshPartition {
    int rank = -1;
    std::vector<int> patchIds;
};

class MultiBlockMesh {
public:
    using CompositePreprocessor =
        std::function<bool(std::vector<RawMeshBlock>&)>;

    bool loadFiles(const std::string& caseDir,
                   const std::vector<std::string>& meshFiles,
                   const MeshRuntimeConfig& config,
                   const CompositePreprocessor& preprocessor = {});

    bool empty() const { return blocks_.empty(); }
    bool isSingleBlock() const { return blocks_.size() == 1; }
    size_t size() const { return blocks_.size(); }

    MeshBlockField& block(size_t i) { return blocks_.at(i); }
    const MeshBlockField& block(size_t i) const { return blocks_.at(i); }

    std::vector<MeshBlockField>& blocks() { return blocks_; }
    const std::vector<MeshBlockField>& blocks() const { return blocks_; }
    const std::vector<MeshPartition>& partitions() const { return partitions_; }
    const MeshPartition& partition(size_t i) const { return partitions_.at(i); }
    size_t partitionCount() const { return partitions_.size(); }
    const MeshCommunication::HaloExchangePlan& haloExchangePlan() const {
        return haloPlan_;
    }

private:
    std::vector<MeshBlockField> blocks_;
    std::vector<MeshPartition> partitions_;
    MeshCommunication::HaloExchangePlan haloPlan_;

    static std::string resolvePath(const std::string& caseDir,
                                   const std::string& meshFile);
    static std::string blockNameFromPath(const std::string& path,
                                         size_t fallbackId);
    static bool readSFMFile(const std::string& path,
                            std::vector<RawMeshBlock>& out);
    static bool appendSetFile(const RawMeshBlock& setFile, RawMeshBlock& target);
    static bool validateRawBlock(const RawMeshBlock& raw);
};

} // namespace SF
