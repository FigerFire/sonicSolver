#pragma once

/// @file SF_levelSetSystemContribution.h
/// @brief Level-set-owned mathematical system contribution.

#include "core/system/SF_systemContribution.h"
#include "core/system/SF_formula.h"

namespace SF::Physics::InterfaceModels::LevelSetContribution {
namespace Mathematics {
inline Math::Signature signature(std::string_view name) {
    if(name=="phi" || name=="phi0")return Math::signature(0,{0,1,0,0,0,0,0});
    if(name=="U")return Math::signature(1,{0,1,-1,0,0,0,0});
    if(name=="interfaceCurvature")return Math::signature(0,{0,-1,0,0,0,0,0});
    if(name=="interfaceNormal")return Math::signature(1,{});
    if(name=="interfaceForce")return Math::signature(1,{1,-2,-2,0,0,0,0});
    if(name=="interfacePower")return Math::signature(0,{1,-1,-3,0,0,0,0});
    if(name=="pLiquidMinusGas")return Math::signature(0,{1,-1,-2,0,0,0,0});
    throw std::runtime_error("Unknown level-set mathematical symbol: "+std::string(name));
}
inline System::FormulaExpr symbol(const std::string& name) {return System::FormulaExpr::symbol(name,signature(name));}
inline System::FormulaExpr op(std::string_view name,std::vector<System::FormulaExpr> args) {
    using S=Math::Signature;
    static const Math::OperatorRule rules[]={
        {"advect",[](const std::vector<S>& a){if(a.size()!=2)throw std::runtime_error("advect needs U/phi");return Math::resultSignature(Math::Op::Dot,{a[0],Math::resultSignature(Math::Op::Grad,{a[1]})});}},
        {"ddtau",[](const std::vector<S>& a){if(a.size()!=1)throw std::runtime_error("ddtau needs phi");Math::requireCompatible(a[0],signature("phi"),"pseudo time has length units");return Math::signature(0,{});}},
        {"gradientNorm",[](const std::vector<S>& a){auto r=Math::resultSignature(Math::Op::Grad,a);r.rank=0;return r;}},
        {"smoothedSign",[](const std::vector<S>& a){if(a.size()!=1)throw std::runtime_error("sign needs phi");Math::requireCompatible(a[0],signature("phi"),"smoothedSign");return Math::signature(0,{});}},
        {"normalizedGradient",[](const std::vector<S>& a){auto r=Math::resultSignature(Math::Op::Grad,a);r.dimensions=Math::Dimensions::dimensionless();return r;}},
        {"meanCurvature",[](const std::vector<S>& a){if(a.size()!=1)throw std::runtime_error("curvature needs phi");Math::requireCompatible(a[0],signature("phi"),"meanCurvature");return signature("interfaceCurvature");}},
        {"regularizedDelta",[](const std::vector<S>& a){if(a.size()!=2)throw std::runtime_error("delta needs phi/width");Math::requireCompatible(a[0],signature("phi"),"delta phi");Math::requireCompatible(a[1],a[0],"delta width");return Math::signature(0,{0,-1,0,0,0,0,0});}}
    };
    for(const auto& rule:rules)if(name==rule.id)return System::FormulaExpr::op(rule,std::move(args));
    throw std::runtime_error("Unknown level-set mathematical operation");
}
}
}

namespace SF::Physics::InterfaceModels::LevelSetContribution {



struct Spec {
    bool ghostFluid = false;
    int advectionOrder = 5;
    int reinitializationOrder = 5;
    int reinitializationSteps = 0;
    double pseudoTimeStep = 0.0;
    double wenoEpsilon = 1e-6;
    double wenoPower = 2.0;
    double signSmoothingFactor = 1.0;
    bool continuousSurfaceForce = false;
    double surfaceTension = 0.0;
    double interfaceThickness = 0.0;
};

void contribute(System::SystemContribution& system, const Spec& spec);

} // namespace SF::Physics::InterfaceModels::LevelSetContribution
