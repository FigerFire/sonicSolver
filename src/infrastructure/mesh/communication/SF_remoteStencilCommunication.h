#pragma once

/// @file SF_remoteStencilCommunication.h
/// @brief IBM 等几何模块请求远端 stencil 数据的 mesh 通信抽象。

#include <string>
#include <vector>

namespace SF::MeshCommunication {

struct RemoteStencilRequest {
    int requesterBlock = -1;
    int donorBlock = -1;
    std::vector<int> donorIndices;
    int components = 1;
    std::string fieldName;
};

struct RemoteStencilData {
    std::vector<double> values;
    int components = 1;
};

/// @brief 与 MPI 无关的远端样本请求；实现由 mesh/infrastructure 装配。
class IRemoteStencilCommunication {
public:
    virtual ~IRemoteStencilCommunication() = default;
    virtual RemoteStencilData gather(
        const RemoteStencilRequest& request) const = 0;
};

} // namespace SF::MeshCommunication
