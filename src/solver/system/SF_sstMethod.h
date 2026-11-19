#pragma once
/// @brief SST 数学实现的不可变输入输出契约；不拥有物理场或执行过程。
#include "SF_methodObjects.h"
#include "solver/discretization/SF_sstClosure.h"
namespace SF::System {
/// @brief 原非局部闭合的明确数值契约，不参与时间积分。
struct CompiledSstClosure {
    std::vector<std::string> targets,outputs;
    std::string density,momentum,viscosity;
    /// @brief 数学操作对应的实际结果索引；不由字段身份推断。
    std::vector<int> outputRoles;
    bool eddyViscosity=false;
    bool packedVelocity=false;
    Discretization::SstParameters parameters;
};
CompiledEquationCall compileSstClosure(const ExecutableEquationSystem&,const EquationCall&,const NumericalBinding&);
}
