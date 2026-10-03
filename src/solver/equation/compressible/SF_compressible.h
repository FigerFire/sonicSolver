#pragma once

/// @file SF_compressible.h
/// @brief 可压缩守恒方程定义及其离散绑定入口。

#include "core/interfaces/SF_termKernel.h"
#include "SF_configTypes.h"


namespace SF {
class Field;
class FluxField;
class Residual;
namespace FDM {
struct SolverConfig;
class ITransportModel;
}
namespace Physics::FluidStateModel { class Model; }
namespace Equation::Compressible {

/// @brief 一次空间项装配所需的非拥有上下文。
struct AssemblyContext {
    const FDM::TermRecipe* convection = nullptr;
    const FDM::TermRecipe* diffusion = nullptr;
    const std::vector<SF::System::ConservativeSourceKernel>* sources = nullptr;
    double timeStep = 0.0;
    FDM::IBMBoundaryScheme ibmBoundary = FDM::IBMBoundaryScheme::LowOrder;
    int ilwOrder = 0;
    double dynamicViscosity = 0.0;
    double prandtl = 0.72;
    double idealGasGamma = 1.4;
    double idealGasConstant = 287.05;
    const FDM::ITransportModel* transport = nullptr;
    const Physics::FluidStateModel::Model* thermodynamics = nullptr;
};

/// @brief 可压缩质量、动量和能量方程组。
class System {
public:
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

};

} // namespace SF::Equation::Compressible
} // namespace SF
