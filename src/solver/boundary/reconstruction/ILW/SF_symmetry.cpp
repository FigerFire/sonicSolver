/// @file SF_symmetry.cpp
/// @brief 对称面：标量/切向速度为偶延拓，法向速度为奇延拓。
#include "SF_symmetry.h"
#include "SF_boundaryClosure.h"
namespace SF::Boundary::ILW::Symmetry {
void applyScalar(Field& field,int i,int j,int k,int vIdx,int accuracyOrder,int axis) {
    BoundaryClosure::writeScalarTaylor(field,i,j,k,axis,accuracyOrder,
        BoundaryClosure::ScalarConstraint::Even,0.,BoundaryClosure::FieldVariableGetter{vIdx},
        [&](int a,int b,int c,double value) { field(a,b,c,vIdx)=value; },"symmetry");
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
            comp==axis ? BoundaryClosure::ScalarConstraint::Odd : BoundaryClosure::ScalarConstraint::Even,
            0.,get,set,"symmetryVector3");
    }
}
} // namespace SF::Boundary::ILW::Symmetry
