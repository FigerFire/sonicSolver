/// @file SF_linear.cpp
/// @brief 普通物理边界条件的线性重建分派器。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_linear.h"

#include "SF_empty.h"
#include "SF_fixedValue.h"
#include "SF_symmetry.h"
#include "SF_zeroGradient.h"

namespace SF {
namespace Boundary {
namespace Algebraic {

void setFixedValueScalar(Field& field, int i, int j, int k,
                         int axis, double bcValue, int vIdx) {
    FixedValue::applyScalar(field, i, j, k, axis, bcValue, vIdx);
}

void setFixedValueVector3(Field& field, int i, int j, int k,
                          int axis, const SF::Vector3& bcValue, int vIdx) {
    FixedValue::applyVector3(field, i, j, k, axis, bcValue, vIdx);
}

void setZeroGradientScalar(Field& field, int i, int j, int k,
                           int axis, int vIdx) {
    ZeroGradient::applyScalar(field, i, j, k, axis, vIdx);
}

void setZeroGradientVector3(Field& field, int i, int j, int k,
                            int axis, int vIdx) {
    ZeroGradient::applyVector3(field, i, j, k, axis, vIdx);
}

void setSymmetryScalar(Field& field, int i, int j, int k,
                       int axis, int vIdx) {
    Symmetry::applyScalar(field, i, j, k, axis, vIdx);
}

void setSymmetryVector3(Field& field, int i, int j, int k,
                        int axis, int vIdx) {
    Symmetry::applyVector3(field, i, j, k, axis, vIdx);
}

void setEmptyScalar(Field& field, int i, int j, int k, int vIdx, int axis) {
    Empty::applyScalar(field, i, j, k, vIdx, axis);
}

void setEmptyVector3(Field& field, int i, int j, int k, int vIdx, int axis) {
    Empty::applyVector3(field, i, j, k, vIdx, axis);
}

} // namespace Algebraic
} // namespace Boundary
} // namespace SF
