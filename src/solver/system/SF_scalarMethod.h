#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::System {
namespace ScalarOps {
inline constexpr const char* Dt="scalar.dt.fixed";
inline constexpr const char* Begin="scalar.snapshot";
inline constexpr const char* Stage="scalar.stage";
inline constexpr const char* Commit="scalar.publish";
inline constexpr const char* Provider="equation.scalar-central2";
inline constexpr const char* TransportProvider="equation.scalar-transport";
inline constexpr const char* CapacityProvider="equation.capacity-transport";
inline constexpr const char* ClosureProvider="closure.sst-nonlocal";
inline constexpr const char* SharedFluxProvider="equation.scalar-shared-flux";
}
/// Immutable AST lowering output. No physical storage or temporal coefficients.
struct CompiledScalarEquation {
    std::string target;
    /// @brief 每个空间 occurrence 独立绑定；空集合表示局部反应。
    struct Diffusion {std::string occurrence,input; double coefficient=0; int components=1; std::string coefficientState;};
    std::vector<Diffusion> diffusions;
    int components=1;
    bool conservativeTransport=false;
    bool sharedMassFlux=false;
    bool primitiveStorage=false;

    std::string diffusivityState;
    std::vector<std::string> sources;
    /// @brief 借用所选版本的输入，不要求每个读取对象都有演化方程。
    std::vector<StateUse> valueReads;
    /// @brief 显式选择节点系数的算术面值；不猜测不连续物性的平均规则。
    bool arithmeticFaces=false;
    FormulaValueKernel rhs;
    std::string normalized;
    bool nonlinear=false;
    std::vector<std::string> lowering;
};
/// @brief 原非局部闭合的明确数值契约，不参与时间积分。
struct CompiledStageClosure {std::vector<std::string> targets;};
/// @brief 针对HOW目标的局部仿射数值实现；WHAT仍保存原等式。
struct CompiledAffineEquation {
    std::string target;
    std::vector<StateUse> valueReads;
    WorkspaceBoundary boundaryPreparation=WorkspaceBoundary::None;
    FormulaValueKernel coefficient,remainder;
};
/// @brief 串行标量块的不可变AST lowering；行标签来自唯一主扩散目标。
struct CompiledLinearBlock {
    std::vector<std::string> equations,targets,rowTargets;
    std::vector<StateUse> valueReads;
    std::vector<FormulaLinearKernel> rows;
    FDM::LinearSolverConfig linear;
};
void addScalarMethod(ProviderRegistry& registry);
}
