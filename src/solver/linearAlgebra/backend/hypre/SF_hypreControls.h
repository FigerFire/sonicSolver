#pragma once
/// @file SF_hypreControls.h
/// @brief 用户声明的 HYPRE AMG 参数映射；预条件周期固定一次、零求解容差。
#include "SF_configTypes.h"
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#include <stdexcept>

namespace SF::LinearAlgebra {
inline void hypreCheck(int code,const char* operation) {
    if(!code)return;
    char description[256]={}; HYPRE_DescribeError(code,description);
    HYPRE_ClearAllErrors();
    throw std::runtime_error(std::string(operation)+": "+description);
}
inline void configureAMGCycle(HYPRE_Solver amg,const FDM::LinearSolverConfig& c) {
    hypreCheck(HYPRE_BoomerAMGSetPrintLevel(amg,0),"AMG print");
    hypreCheck(HYPRE_BoomerAMGSetMaxIter(amg,1),"AMG cycle");
    hypreCheck(HYPRE_BoomerAMGSetTol(amg,0),"AMG tolerance");
    hypreCheck(HYPRE_BoomerAMGSetCoarsenType(amg,c.amgCoarsenType),"AMG coarsening");
    hypreCheck(HYPRE_BoomerAMGSetRelaxType(amg,c.amgRelaxType),"AMG relaxation");
    hypreCheck(HYPRE_BoomerAMGSetNumSweeps(amg,c.amgSweeps),"AMG sweeps");
    hypreCheck(HYPRE_BoomerAMGSetMaxLevels(amg,c.amgMaxLevels),"AMG levels");
    hypreCheck(HYPRE_BoomerAMGSetStrongThreshold(amg,c.amgStrongThreshold),"AMG strength");
}
} // namespace SF::LinearAlgebra
