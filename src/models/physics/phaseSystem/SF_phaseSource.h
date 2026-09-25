#pragma once

/// @file SF_phaseSource.h
/// @brief 与方程装配解耦的 Eulerian 相源项接口和注册表。

#include "SF_phaseSystem.h"

#include <memory>
#include <string>
#include <vector>

namespace SF::Physics::PhaseSystems {

/// @brief 一个可插拔相方程贡献模块。
class PhaseEquationSource {
public:
    virtual ~PhaseEquationSource() = default;
    /// @brief 稳定模块名，用于启动日志和错误定位。
    virtual std::string name() const = 0;
    /// @brief 向统一相方程贡献容器累加源项。
    virtual void add(
        const PhaseSystem& system,
        double dt,
        PhaseEquationSources& sources) const = 0;
};

/// @brief 外部相方程源项注册表；不识别具体物理模块。
class PhaseSourceRegistry {
public:
    /// @brief 注册一个已配置的 source provider。
    void add(std::unique_ptr<PhaseEquationSource> source);
    /// @brief 依注册顺序累加所有外部源项。
    void assemble(PhaseSystem& system, double dt) const;
    /// @brief 返回注册模块名，供日志显示。
    std::vector<std::string> names() const;

private:
    std::vector<std::unique_ptr<PhaseEquationSource>> sources_;
};

} // namespace SF::Physics::PhaseSystems
