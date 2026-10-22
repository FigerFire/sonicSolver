#pragma once

/// @file SF_pressureCoupling.h
/// Pressure coupling contributes mathematics, generic HOW and independent WHICH.
/// Single-fluid and Eulerian counts/order live only in their native HOW.

#include "SF_couplingStatus.h"
#include "core/system/SF_solveProgram.h"
#include "core/config/types/SF_pressureConfigTypes.h"

#include <string>
#include <vector>

namespace SF::System {
class SystemCompositionBuilder;

/// @brief 已有预测器的表示；不从 equation/occurrence 名称推断存储或算法。
enum class PressurePredictorKind { Velocity, ConservativeMomentum };
/// @brief 方法输入引用同一 WHAT/HOW 对象；只在 composition 时使用，不拥有状态。
struct PressureCouplingInputs {
    EquationCall predictor;
    PressurePredictorKind representation=PressurePredictorKind::Velocity;
    Target predictorTarget;
    std::string predictorOccurrence;
    std::string constraint;
    std::vector<EquationCall> fusedBalances;
    std::string replacedCommit;
    EquationCall correction,velocityUpdate,pressureUpdate,fluxUpdate,publication;
    EquationCall relaxation,fluxConsistency,convergence;
};
/// @brief 内置方案的声明默认；显式 C++ 组装可以提供其他 equation/occurrence 引用。
PressureCouplingInputs pressureCouplingInputs(EquationCall predictor,PressurePredictorKind);


namespace PressureMathematics {
inline Math::Signature signature(std::string_view name,bool conservative=false) {
    if(name=="rho")return Math::signature(0,{1,-3,0,0,0,0,0});
    if(name=="U")return Math::signature(1,{0,1,-1,0,0,0,0});
    if(name=="rhoU")return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name=="rhoE" || name=="p" || name=="pPrime")return Math::signature(0,{1,-1,-2,0,0,0,0});
    if(name=="soundSpeed")return Math::signature(0,{0,1,-1,0,0,0,0});
    if(name=="dt")return Math::signature(0,{0,0,1,0,0,0,0});
    if(name=="momentumRelaxation" || name=="pressureRelaxation" || name=="fluxValidity" || name=="converged")return Math::signature(0,{});
    if(name=="rAU" || name=="faceResponse")return Math::signature(0,{-1,3,1,0,0,0,0});
    if(name=="phi" || name=="correctedFlux")return Math::signature(0,{0,3,-1,0,0,0,0});
    if(name=="iterate" || name=="laggedIterate")return conservative
        ? Math::product({signature("rho"),signature("rhoU"),signature("rhoE")})
        : Math::product({signature("U"),signature("p"),signature("phi")});
    throw std::runtime_error("Unknown pressure mathematical symbol: "+std::string(name));
}
inline FormulaExpr symbol(const std::string& id,bool conservative=false) {return FormulaExpr::symbol(id,signature(id,conservative));}
inline FormulaExpr op(std::string_view name,std::vector<FormulaExpr> args,std::string occurrence={}) {
    using S=Math::Signature;
    static const Math::OperatorRule rules[]={
        {"faceAverage",[](const std::vector<S>& a){if(a.size()!=1 || a[0].rank)throw std::runtime_error("Face average requires scalar");return a[0];}},
        {"faceVolumeDivergence",[](const std::vector<S>& a){if(a.size()!=1)throw std::runtime_error("Face divergence needs one volume rate");Math::requireCompatible(a[0],signature("phi"),"face rate");return Math::signature(0,{0,0,-1,0,0,0,0});}},
        {"pressureFlux",[](const std::vector<S>& a){if(a.size()!=2)throw std::runtime_error("Pressure flux needs response/p");Math::requireCompatible(a[0],signature("faceResponse"),"pressure response");Math::requireCompatible(a[1],signature("pPrime"),"pressure flux");return signature("phi");}},
        {"invalidateDerivedFlux",[](const std::vector<S>& a){if(a.size()!=1)throw std::runtime_error("Flux invalidation needs momentum");Math::requireCompatible(a[0],signature("rhoU"),"flux invalidation");return signature("fluxValidity");}},
        {"totalEnergyFromPressure",[](const std::vector<S>& a){if(a.size()!=3)throw std::runtime_error("Energy closure needs rho/rhoU/p");for(unsigned i=0;i<3;++i)Math::requireCompatible(a[i],signature(i==0?"rho":i==1?"rhoU":"p"),"energy closure");return signature("rhoE");}},
        {"relax",[](const std::vector<S>& a){if(a.size()!=2 || !a[0].tupleSize)throw std::runtime_error("Relaxation requires typed iterate products");Math::requireCompatible(a[0],a[1],"relax");return a[0];}},
        {"residualConvergence",[](const std::vector<S>& a){if(a.size()!=2 || !a[0].tupleSize)throw std::runtime_error("Convergence requires typed iterate products");Math::requireCompatible(a[0],a[1],"convergence");return signature("converged");}},
        {"consistentFlux",[](const std::vector<S>& a){if(a.size()!=2)throw std::runtime_error("Flux consistency needs velocity/momentum and pressure");Math::requireCompatible(a[1],signature("p"),"flux pressure");if(Math::compatible(a[0],signature("U")))return signature("phi");Math::requireCompatible(a[0],signature("rhoU"),"flux momentum");return signature("fluxValidity");}}
    };
    for(const auto& r:rules)if(name==r.id)return FormulaExpr::op(r,std::move(args),std::move(occurrence));
    throw std::runtime_error("Unknown pressure mathematical relation: "+std::string(name));
}
}

/// @brief 既有 conservative Corrector 的数学关系；同一描述用于注册与 provider 验证。
std::vector<Equation> conservativePressureRelations();

/// @brief 匹配 requirements，不做任何 contribution。
CouplingReport matchPressureCoupling(
    const CouplingPresetRequest& request,
    const RawEquationSystem& raw,
    const PressureCouplingInputs* inputs=nullptr);

/// @brief 把 active preset 的 formulation/plan 贡献写入 composition。
///
/// 只有 status==Active 时才写入 transformation 请求；shared Eulerian 使用自身原生 HOW；
/// inactive/invalid/unsupported 一律不写，并由调用方按状态报告或 fail。
/// @return 同一个 report（带上 derived operations/equations/providers）。
CouplingReport contributePressureCoupling(
    SystemCompositionBuilder& system,
    const CouplingPresetRequest& request,
    const RawEquationSystem& raw,
    const PressureCouplingInputs* inputs=nullptr);

/// @brief 从 case 输入解析 coupling preset 请求。
CouplingPresetRequest couplingRequestFrom(
    const FDM::PressureCouplingConfig& coupling,
    bool explicitlyRegistered);

/// Single-fluid coupling contributes HOW directly. The
/// EquationMethods, not this tree, supply equation-level numerical leaves.
void applyPressureExecution(ExecutionProgram& program,
                            const CouplingPresetRequest& request,
                            const PressureCouplingInputs& inputs);
/// Independent WHICH contribution; never mutates HOW.
std::vector<NumericalBinding> pressureNumerics(const CouplingPresetRequest& request,
    const PressureCouplingInputs& inputs);

} // namespace SF::System
