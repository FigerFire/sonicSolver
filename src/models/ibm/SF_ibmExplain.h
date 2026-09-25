#pragma once

/// @file SF_ibmExplain.h
/// @brief IBM 算法解释快照：由 application inspection 持有，不属于 resolved
/// numerical system 的 authority。
///
/// `ResolvedSimulationSystem` 只聚合 WHAT/ORDER/HOW/RUNTIME 的数学身份；
/// IBM 的 `explain` 元数据是模型自己的展示信息，由 `ImmersedAlgorithmDescriptor`
/// 直接派生，因此单独作为 value object 放在 application 快照中。

#include "SF_immersedSystem.h"

#include <string>

namespace SF::IBM {

/// @brief 只用于 `explain` 输出的 IBM 算法摘要，不携带数值状态。
struct ImmersedExplain {
    std::string algorithm;
    std::string reference;
    std::string support;
    std::string representation;
    std::string enforcement;
    std::string solid;
    std::string functional;

    bool empty() const {
        return algorithm.empty() && reference.empty() && support.empty()
            && representation.empty() && enforcement.empty() && solid.empty()
            && functional.empty();
    }
};

/// @brief 从解析好的数学描述派生解释摘要。
inline ImmersedExplain explainSnapshot(
        const FDM::ImmersedAlgorithmDescriptor& immersed) {
    ImmersedExplain result;
    result.algorithm = immersed.id;
    result.reference = immersed.referenceName;
    result.support = FDM::toString(immersed.support);
    result.representation = FDM::toString(immersed.representation);
    result.enforcement = FDM::toString(immersed.enforcement);
    result.solid = FDM::toString(immersed.solid);
    result.functional = immersed.variational.stationaryFunctional;
    return result;
}

/// @brief 渲染 `explain` 中的 IBM 段；空快照返回空字符串。
inline std::string renderExplain(const ImmersedExplain& explain) {
    if (explain.empty()) return std::string();
    std::string output = "\nIBM\n"
        "  algorithm      : " + explain.algorithm + "\n"
        "  reference      : " + explain.reference + "\n"
        "  support        : " + explain.support + "\n"
        "  representation : " + explain.representation + "\n"
        "  enforcement    : " + explain.enforcement + "\n"
        "  solid          : " + explain.solid + "\n";
    if (!explain.functional.empty()) {
        output += "  functional     : " + explain.functional + "\n";
    }
    return output;
}

} // namespace SF::IBM
