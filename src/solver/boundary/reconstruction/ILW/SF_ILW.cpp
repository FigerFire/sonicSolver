/// @file SF_ILW.cpp
/// @brief ILW 边界重建方法族的类型调度入口。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#include "SF_ILW.h"

#include "SF_empty.h"
#include "SF_fixedValue.h"
#include "SF_symmetry.h"
#include "SF_zeroGradient.h"

namespace SF {
namespace Boundary {
namespace ILW {

void setFixedValueScalar(Field& field, int i, int j, int k,
                         double bcValue, int vIdx, int accuracyOrder, int axis) {
    FixedValue::applyScalar(
        field, i, j, k, bcValue, vIdx, accuracyOrder, axis);
}

void setFixedValueVector3(Field& field, int i, int j, int k,
                          const SF::Vector3& bcValue, int vIdx,
                          int accuracyOrder, int axis) {
    FixedValue::applyVector3(
        field, i, j, k, bcValue, vIdx, accuracyOrder, axis);
}

void setZeroGradientScalar(Field& field, int i, int j, int k, int vIdx,
                           int accuracyOrder, int axis) {
    ZeroGradient::applyScalar(field, i, j, k, vIdx, accuracyOrder, axis);
}

void setZeroGradientVector3(Field& field, int i, int j, int k, int vIdx,
                            int accuracyOrder, int axis) {
    ZeroGradient::applyVector3(field, i, j, k, vIdx, accuracyOrder, axis);
}

void setSymmetryScalar(Field& field, int i, int j, int k, int vIdx,
                       int accuracyOrder, int axis) {
    Symmetry::applyScalar(field, i, j, k, vIdx, accuracyOrder, axis);
}

void setSymmetryVector3(Field& field, int i, int j, int k, int vIdx,
                        int accuracyOrder, int axis) {
    Symmetry::applyVector3(field, i, j, k, vIdx, accuracyOrder, axis);
}

void setEmptyScalar(Field& field, int i, int j, int k, int vIdx, int axis) {
    Empty::applyScalar(field, i, j, k, vIdx, axis);
}

void setEmptyVector3(Field& field, int i, int j, int k, int vIdx, int axis) {
    Empty::applyVector3(field, i, j, k, vIdx, axis);
}

} // namespace ILW
} // namespace Boundary
} // namespace SF
