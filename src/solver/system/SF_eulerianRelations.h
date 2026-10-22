#pragma once
/// @file SF_eulerianRelations.h
/// @brief 现有双欧拉内核的数学关系；不含时间步或耦合顺序。
#include "core/system/SF_equationIR.h"
namespace SF::System {
namespace EulerianMathematics {
inline Math::Signature signature(std::string name) {
    name=name.substr(0,name.find('.'));
    if(name=="rho" || name=="phaseMass" || name=="faceDensity")return Math::signature(0,{1,-3,0,0,0,0,0});
    if(name=="U" || name=="predictedVelocity")return Math::signature(1,{0,1,-1,0,0,0,0});
    if(name=="momentum")return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name=="enthalpy" || name=="p" || name=="pPrime")return Math::signature(0,{1,-1,-2,0,0,0,0});
    if(name=="T")return Math::signature(0,{0,0,0,1,0,0,0});
    if(name=="h")return Math::signature(0,{0,2,-2,0,0,0,0});
    if(name=="alpha" || name=="referencePhaseRemainder" || name=="pressureRelaxation")return Math::signature(0,{});
    if(name=="effectivePhaseViscosity" || name=="effectivePhaseConductivity")return Math::signature(0,{1,-1,-1,0,0,0,0});
    if(name=="phaseMassFlux")return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name=="phaseMomentumFlux")return Math::signature(2,{1,-1,-2,0,0,0,0});
    if(name=="phaseEnthalpyFlux")return Math::signature(1,{1,0,-3,0,0,0,0});
    if(name=="phaseMassSources")return Math::signature(0,{1,-3,-1,0,0,0,0});
    if(name=="phaseMomentumSources")return Math::signature(1,{1,-2,-2,0,0,0,0});
    if(name=="phaseEnthalpySources" || name=="materialPressureRate")return Math::signature(0,{1,-1,-3,0,0,0,0});
    if(name=="momentumDiagonal")return Math::signature(0,{1,-3,-1,0,0,0,0});
    if(name=="volumeFaceFlux")return Math::signature(0,{0,3,-1,0,0,0,0});
    if(name=="massFaceFlux")return Math::signature(0,{1,0,-1,0,0,0,0});
    if(name=="facePressureResponse")return Math::signature(0,{-1,5,1,0,0,0,0});
    if(name=="mixtureCompressibilityOverDt")return Math::signature(0,{-1,1,1,0,0,0,0});
    if(name=="sumFaceAlphaSquaredOverDiagonal")return Math::signature(0,{-1,3,1,0,0,0,0});
    if(name=="sumMassSourceOverDensity" || name=="pressureHistory")return Math::signature(0,{0,0,-1,0,0,0,0});
    if(name=="mixtureVolumeFlux")return Math::signature(1,{0,1,-1,0,0,0,0});
    throw std::runtime_error("Unknown Eulerian mathematical symbol: "+name);
}
inline FormulaExpr symbol(const std::string& name){return FormulaExpr::symbol(name,signature(name));}
inline FormulaExpr op(std::string_view name,std::vector<FormulaExpr> args) {
    static const Math::OperatorRule pair{"pairedFaceFlux",[](const std::vector<Math::Signature>& a){
        if(a.size()!=2)throw std::runtime_error("Paired face flux needs two rates");
        Math::requireCompatible(a[0],signature("volumeFaceFlux"),"volume rate");Math::requireCompatible(a[1],signature("massFaceFlux"),"mass rate");
        return Math::product({a[0],a[1]});
    }};
    static const Math::OperatorRule gradient{"faceGradient",[](const std::vector<Math::Signature>& a){
        if(a.size()!=1)throw std::runtime_error("Face gradient needs pressure");Math::requireCompatible(a[0],signature("pPrime"),"face gradient");
        auto result=Math::resultSignature(Math::Op::Grad,a);result.rank=0;return result;
    }};
    if(name==pair.id)return FormulaExpr::op(pair,std::move(args));
    if(name==gradient.id)return FormulaExpr::op(gradient,std::move(args));
    return FormulaExpr::op(Math::resolveOperator(name),std::move(args));
}
}
std::vector<Equation> eulerianPhaseRelations(const std::string& phase,bool reference);
std::vector<Equation> eulerianPressureRelations();
}
