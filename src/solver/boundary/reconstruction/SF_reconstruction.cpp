/// @file SF_reconstruction.cpp
/// @brief 根据用户选择创建线性或 ILW 边界重建器。

#include "reconstruction/SF_reconstruction.h"

#include "reconstruction/ILW/SF_ILW.h"
#include "reconstruction/linear/SF_linear.h"
#include "SF_numericsPolicy.h"

#include <stdexcept>

namespace SF::Boundary::Reconstruction {
namespace {

[[noreturn]] void unsupported(Kind kind) {
    const char* name = kind == Kind::Polynomial ? "Polynomial" : "MLS";
    throw std::runtime_error(
        std::string("Boundary reconstruction ") + name
        + " is declared but has no migrated Field adapter. Select Linear "
          "or explicitly configure ILW; no implicit fallback is applied.");
}

} // namespace

Selection fromILWSetting(bool enabled, int order) {
    if (!enabled) return {Kind::Linear, 0};
    if (!FDM::isSupportedILWOrder(order) || order == 0) {
        throw std::runtime_error(
            "ILW boundary reconstruction requires explicit order 3/5/7/9.");
    }
    return {Kind::ILW, order};
}

void applyScalar(Field& field, int i, int j, int k,
                 int axis, int variable, double prescribedValue,
                 Law::Kind law, Selection selection) {
    if (selection.kind == Kind::Polynomial || selection.kind == Kind::MLS) {
        unsupported(selection.kind);
    }
    const bool ilw = selection.kind == Kind::ILW;
    switch (law) {
        case Law::Kind::FixedValue:
            if (ilw) ILW::setFixedValueScalar(
                field,i,j,k,prescribedValue,variable,selection.order);
            else Algebraic::setFixedValueScalar(
                field,i,j,k,axis,prescribedValue,variable);
            return;
        case Law::Kind::ZeroGradient:
            if (ilw) ILW::setZeroGradientScalar(
                field,i,j,k,variable,selection.order);
            else Algebraic::setZeroGradientScalar(
                field,i,j,k,axis,variable);
            return;
        case Law::Kind::Symmetry:
            if (ilw) ILW::setSymmetryScalar(
                field,i,j,k,variable,selection.order);
            else Algebraic::setSymmetryScalar(
                field,i,j,k,axis,variable);
            return;
        case Law::Kind::Empty:
            if (ilw) ILW::setEmptyScalar(field,i,j,k,variable,axis);
            else Algebraic::setEmptyScalar(field,i,j,k,variable,axis);
            return;
    }
    throw std::runtime_error("Unknown scalar boundary law.");
}

void applyVector(Field& field, int i, int j, int k,
                 int axis, int variable, const Vector3& prescribedValue,
                 Law::Kind law, Selection selection) {
    if (selection.kind == Kind::Polynomial || selection.kind == Kind::MLS) {
        unsupported(selection.kind);
    }
    const bool ilw = selection.kind == Kind::ILW;
    switch (law) {
        case Law::Kind::FixedValue:
            if (ilw) ILW::setFixedValueVector3(
                field,i,j,k,prescribedValue,variable,selection.order);
            else Algebraic::setFixedValueVector3(
                field,i,j,k,axis,prescribedValue,variable);
            return;
        case Law::Kind::ZeroGradient:
            if (ilw) ILW::setZeroGradientVector3(
                field,i,j,k,variable,selection.order);
            else Algebraic::setZeroGradientVector3(
                field,i,j,k,axis,variable);
            return;
        case Law::Kind::Symmetry:
            if (ilw) ILW::setSymmetryVector3(
                field,i,j,k,variable,selection.order);
            else Algebraic::setSymmetryVector3(
                field,i,j,k,axis,variable);
            return;
        case Law::Kind::Empty:
            if (ilw) ILW::setEmptyVector3(field,i,j,k,variable,axis);
            else Algebraic::setEmptyVector3(field,i,j,k,variable,axis);
            return;
    }
    throw std::runtime_error("Unknown vector boundary law.");
}

} // namespace SF::Boundary::Reconstruction
