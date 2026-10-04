#pragma once
/// @file SF_eulerianAssembly.h
/// @brief 由 AST/STATE/WHICH 编译的双欧拉实现契约；不是另一份数学方程。
#include "core/system/SF_solveProgram.h"
#include <memory>

namespace SF::System {
/// @brief 已校验的专用数值关系；不按 ddt/div 等平面项分类。
enum class EulerianAssemblyRelation {
    PhaseMassBalance, ReferencePhaseClosure, PhaseMomentumBalance,
    SharedPressureCorrection, PhaseMomentumCorrection, SharedPressureUpdate,
    PairedFaceFluxCorrection, PhaseEnthalpyBalance
};
/// @brief 已支持且由 AST 显式声明的源项扩展。
enum class EulerianSourceExtension { Gravity, MRF, WallHeat };
/// @brief 编译产物只冻结实现/相槽/输出绑定，数值算术仍在原专用内核。
struct CompiledEulerianAssemblyContract {
    EquationRef equation;
    Target target;
    std::string method;
    std::string provider;
    EulerianAssemblyRelation relation;
    std::optional<std::size_t> phaseSlot;
    std::string phaseName;
    std::vector<EulerianSourceExtension> extensions;
};
using EulerianAssemblyContractPtr=std::shared_ptr<const CompiledEulerianAssemblyContract>;
/// @brief 只读取冻结的 provider 编译产物；不存在时返回 nullptr。
inline const CompiledEulerianAssemblyContract* eulerianAssemblyContract(const CompiledEquationCall& call) {
    const auto* pointer=std::any_cast<EulerianAssemblyContractPtr>(&call.providerContract);
    return pointer ? pointer->get() : nullptr;
}
const char* toString(EulerianAssemblyRelation relation);
const char* toString(EulerianSourceExtension extension);
/// @brief 启动时把冻结相槽与实际 PhaseSystem 对齐；不解释 AST，也不进入数值循环。
void validateEulerianAssemblyBindings(const CompiledExecutionProgram& program,
    const std::vector<std::string>& phaseNames,std::size_t referencePhase);
}
