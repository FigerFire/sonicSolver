#pragma once
/// @file SF_immersedTurbulenceBoundary.h
/// @brief 稳定注入的壁面数值端口，直接写模型原数组，不刷新流动或调度时间。
#include "core/field/SF_field.h"
#include <vector>
namespace SF::FDM {
class IImmersedTurbulenceBoundary {
public:
    virtual ~IImmersedTurbulenceBoundary()=default;
    /// @brief 在已经 boundary-ready 的速度上施加 SST k/omega/mu_t ghost 条件。
    virtual void applySST(const Field& field,double laminarMu,std::vector<double>& k,
                          std::vector<double>& omega,std::vector<double>& muT) const=0;
};
}
