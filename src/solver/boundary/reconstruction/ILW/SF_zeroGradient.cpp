/// @file SF_zeroGradient.cpp
/// @brief ILW ZeroGradient：从当前内点发布实边界与高阶 Taylor ghost。
#include "SF_zeroGradient.h"
#include "SF_boundaryClosure.h"
namespace SF::Boundary::ILW::ZeroGradient {
void applyScalar(Field& field,int i,int j,int k,int vIdx,int accuracyOrder,int axis) {
    BoundaryClosure::writeScalarTaylor(field,i,j,k,axis,accuracyOrder,
        BoundaryClosure::ScalarConstraint::NormalGradient,0.,
        BoundaryClosure::FieldVariableGetter{vIdx},
        [&](int a,int b,int c,double value) { field(a,b,c,vIdx)=value; },"ZeroGradient");
}
void applyVector3(Field& field,int i,int j,int k,int vIdx,int accuracyOrder,int axis) {
    for (int comp=0;comp<3;++comp) {
        auto get=[&](const Field& f,int a,int b,int c) {
            return vIdx==RU ? BoundaryClosure::VelocityComponentGetter{comp}(f,a,b,c)
                           : f(a,b,c,vIdx+comp);
        };
        auto set=[&](int a,int b,int c,double value) {
            field(a,b,c,vIdx+comp)=value*(vIdx==RU ? BoundaryClosure::velocityDensity(field,a,b,c) : 1.);
        };
        BoundaryClosure::writeScalarTaylor(field,i,j,k,axis,accuracyOrder,
            BoundaryClosure::ScalarConstraint::NormalGradient,0.,get,set,"ZeroGradientVector3");
    }
}
} // namespace SF::Boundary::ILW::ZeroGradient
