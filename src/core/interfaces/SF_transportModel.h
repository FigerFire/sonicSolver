#pragma once

/// @file SF_transportModel.h
/// @brief 黏性、湍流及界面物性向求解器提供的统一输运闭合接口。

#include "core/field/SF_field.h"

#include <string>
#include <vector>

namespace SF::FDM {

/// @brief 求解器可见的有效输运闭合接口。
class ITransportModel {
public:
    virtual ~ITransportModel() = default;
    virtual void applyBoundary(const Field& field) = 0;
    virtual void correct(const Field& field, double dt) = 0;
    virtual double dynamicViscosity(const Field& field,
                                    int i, int j, int k,
                                    double laminarMu) const = 0;

    /// @brief correct() 前需要新鲜 halo 的模型自管字段。
    virtual std::vector<std::string> distributedReadFields() const {
        return {};
    }
    /// @brief correct()/边界处理会写入的模型自管字段。
    virtual std::vector<std::string> distributedWriteFields() const {
        return {};
    }
    /// @brief 上述模型 stencil 的实际 halo 深度；由模型而非注册点声明。
    virtual int distributedHaloDepth() const { return 0; }
};

} // namespace SF::FDM
