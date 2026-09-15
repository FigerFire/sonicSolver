#pragma once

/// @file SF_compressible.h
/// @brief 可压缩守恒方程定义及其离散绑定入口。

#include "solver/equation/SF_expression.h"

#include <memory>

namespace SF {
class Field;
class FluxField;
class Residual;
namespace FDM {
struct SolverConfig;
class ITransportModel;
}
namespace Physics::EquationSet { class Model; }
namespace Equation::Compressible {

/// @brief 一次空间项装配所需的非拥有上下文。
struct AssemblyContext {
    const FDM::SolverConfig& config;
    double timeStep = 0.0;
    const FDM::ITransportModel* transport = nullptr;
    const Physics::EquationSet::Model* thermodynamics = nullptr;
};

/// @brief 可压缩质量、动量和能量方程组。
class System {
public:
    System();

    /// @brief 返回不含格式名称的方程声明。
    const Equation::System& definition() const { return definition_; }

    /// @brief 清空旧残差，开始当前 stage 的方程装配。
    void begin(FluxField& fluxField, Residual& residual) const;
    /// @brief 离散所有对流通量项。
    void convection(Field& field, FluxField& fluxField, Residual& residual,
                    const AssemblyContext& context) const;
    /// @brief 离散扩散与体源项。
    void diffusionAndSources(
        Field& field, Residual& residual, const AssemblyContext& context) const;
    /// @brief 单场便捷入口：依次执行 begin、convection、diffusion/source。
    void assemble(Field& field, FluxField& fluxField, Residual& residual,
                  const AssemblyContext& context) const;

private:
    Equation::System definition_;
};

} // namespace Equation::Compressible
} // namespace SF
