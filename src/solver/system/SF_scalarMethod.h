#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::System {
/// Immutable AST lowering output. No physical storage or temporal coefficients.
struct CompiledScalarEquation {
    std::string target;
    /// @brief 每个空间 occurrence 独立绑定；空集合表示局部反应。
    struct Diffusion {std::string occurrence,input; double coefficient=0; int components=1; std::string coefficientState;};
    std::vector<Diffusion> diffusions;
    int components=1;
    bool conservativeTransport=false;
    /// @brief 显式 HOW Step 的一次 ForwardEuler 求值/发布；不拥有 clock。
    bool sequential=false;
    bool fixedTime=false;double relaxation=1.0;
    /// @brief 显式metric实现及原集合边界/bounds选项；非热模型身份。
    std::optional<FDM::ScalarBoundaryOptions> metric;
    std::string volumeFlux,fluxVelocity;
    bool sharedMassFlux=false;
    bool primitiveStorage=false;

    std::string diffusivityState;
    std::vector<std::string> sources;
    /// @brief 借用所选版本的输入，不要求每个读取对象都有演化方程。
    std::vector<StateUse> valueReads;
    /// @brief 真正读取相邻节点的操作数；局部代数输入没有此需求。
    std::vector<StateUse> stencilReads;
    /// @brief 显式选择节点系数的算术面值；不猜测不连续物性的平均规则。
    bool arithmeticFaces=false;
    FormulaValueKernel rhs;
    std::string normalized;
    bool nonlinear=false;
    std::vector<std::string> lowering;
};
std::optional<CompiledEquationCall> compileScalarEquation(const ExecutableEquationSystem&,const EquationCall&,const NumericalBinding&);
}
