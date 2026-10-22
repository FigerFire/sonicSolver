#pragma once
/// @brief 既有测试的显式数学声明；不是生产reader的兼容或默认路径。
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_pressureCoupling.h"
#include "models/turbulence/SF_transportMathematics.h"
namespace SF::Test {
inline System::FormulaExpr scalar(const std::string& id) {return System::FormulaExpr::symbol(id,Math::signature(0,{}));}
inline System::FormulaExpr symbol(const std::string& id) {
    const System::BuiltinStateCatalog builtin;
    if(builtin.contains(id))return builtin.expression(id);
    if(id=="rhoC")return System::FormulaExpr::symbol(id,Math::signature(0,{1,-3,0,0,0,0,0}));
    if(id=="D" || id=="materialD")return System::FormulaExpr::symbol(id,Math::signature(0,{0,2,-1,0,0,0,0}));
    if(id=="mu_t")return Turbulence::mathematicalSymbol(id);
    if(id=="Sc")return System::FormulaExpr::symbol(id,Math::signature(0,{1,-3,-1,0,0,0,0}));
    if(id=="gravity" || id=="MRF" || id=="customForce" || id=="unauthorizedSource" || id=="unknownSource")return System::FormulaExpr::symbol(id,Math::signature(1,{1,-2,-2,0,0,0,0}));
    if(id=="wallHeat")return System::FormulaExpr::symbol(id,Math::signature(0,{1,-1,-3,0,0,0,0}));
    if(id=="flux" || id=="scalarFlux")return System::FormulaExpr::symbol(id,Math::signature(1,{0,1,-1,0,0,0,0}));
    if(id=="userFluxWithoutProvider")return System::FormulaExpr::symbol(id,Math::signature(1,{1,-2,-1,0,0,0,0}));
    if(id=="unimplementedFlux")return System::FormulaExpr::symbol(id,Math::signature(2,{1,-1,-2,0,0,0,0}));
    if(id=="S")return System::FormulaExpr::symbol(id,Math::signature(0,{0,0,-1,0,0,0,0}));
    if(id=="pPrime")return System::PressureMathematics::symbol(id);
    if(id.rfind("SST.",0)==0 || id.rfind("k.",0)==0 || id.rfind("omega.",0)==0 || id.rfind("epsilon.",0)==0 || id.rfind("kOmegaSST.",0)==0 || id.rfind("kEpsilon.",0)==0)
        return Turbulence::mathematicalSymbol(id);
    if(id.find("eddyDynamicViscosity.")!=std::string::npos)return System::FormulaExpr::symbol(id,Math::signature(0,{1,-1,-1,0,0,0,0}));
    if(id=="linearScalarForce" || id=="linearScalarWork") {
        const auto result=id=="linearScalarForce"?Math::signature(1,{1,-2,-2,0,0,0,0}):Math::signature(0,{1,-1,-3,0,0,0,0});
        return System::FormulaExpr::symbol(id,Math::function(result,{Math::signature(0,{1,-3,0,0,0,0,0}),Math::signature(0,{0,1,-2,0,0,0,0}),Math::signature(0,{})}));
    }
    // Explicit dimensionless declarations for these metadata-only fixtures.
    for(const auto* declared:{"A","B","C","Y_O2","tracer","wrongRelaxation","V","value0","value1","value2","otherValue","U*","a","b","e","alpha.liquid","m","q","x","y","z","packed","unrelated","c","s","kappa","axis","userScalar","unused","missing"})
        if(id==declared)return scalar(id);
    throw std::runtime_error("Test fixture lacks an explicit mathematical declaration: "+id);
}
inline System::StateSymbol state(System::StateSymbol value) {
    if(!value.signature.dimensions.known && !value.signature.tupleSize) {
        value.signature=symbol(value.id).signature();
        if(value.signature.rank==0 && value.shape!=System::ValueShape::Scalar) {
            value.signature.rank=value.shape==System::ValueShape::Vector?1:2;
            if(value.shape==System::ValueShape::Vector)value.signature.extent=value.components;
        }
    }
    return value;
}
inline System::FormulaExpr replaceAt(System::FormulaExpr expression,std::vector<std::size_t> path,System::FormulaExpr replacement) {
    if(path.empty())return replacement;
    const auto i=path.front();path.erase(path.begin());auto children=expression.arguments();
    children.at(i)=replaceAt(children.at(i),std::move(path),std::move(replacement));return expression.withArguments(std::move(children));
}

}
