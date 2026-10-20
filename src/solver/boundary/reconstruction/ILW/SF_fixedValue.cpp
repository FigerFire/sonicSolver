/// @file SF_fixedValue.cpp
/// @brief ILW FixedValue：从当前内点发布实边界与高阶 Taylor ghost。
#include "SF_fixedValue.h"
#include "SF_boundaryClosure.h"
namespace SF::Boundary::ILW::FixedValue {
void applyScalar(Field& field,int i,int j,int k,double bcValue, int vIdx,int accuracyOrder,int axis) {
    BoundaryClosure::writeScalarTaylor(field,i,j,k,axis,accuracyOrder,
        BoundaryClosure::ScalarConstraint::Value,bcValue,
        BoundaryClosure::FieldVariableGetter{vIdx},
        [&](int a,int b,int c,double value) { field(a,b,c,vIdx)=value; },"FixedValue");
}
void applyVector3(Field& field,int i,int j,int k,const Vector3& bcValue, int vIdx,int accuracyOrder,int axis) {
    const double fixed[3]={bcValue.x,bcValue.y,bcValue.z};
    for (int comp=0;comp<3;++comp) {
        auto get=[&](const Field& f,int a,int b,int c) {
            return vIdx==RU ? BoundaryClosure::VelocityComponentGetter{comp}(f,a,b,c)
                           : f(a,b,c,vIdx+comp);
        };
        auto set=[&](int a,int b,int c,double value) {
            field(a,b,c,vIdx+comp)=value*(vIdx==RU ? BoundaryClosure::velocityDensity(field,a,b,c) : 1.);
        };
        BoundaryClosure::writeScalarTaylor(field,i,j,k,axis,accuracyOrder,
            BoundaryClosure::ScalarConstraint::Value,fixed[comp],get,set,"FixedValueVector3");
    }
}
} // namespace SF::Boundary::ILW::FixedValue
