/// @file SF_scalarField.cpp
/// @brief canonical 场数据、索引与状态存储实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.02-----------*/

#include "SF_scalarField.h"

namespace SF {

void ScalarField::setupLike(const Field& field,
                            const std::string& name,
                            double fillValue) {
    name_ = name;
    mx_ = field.MX();
    my_ = field.MY();
    mz_ = field.MZ();
    ng_ = field.NG();
    totalSize_ = field.TotalSize();
    values_.assign((size_t)totalSize_, fillValue);
}

bool ScalarField::isCompatibleWith(const Field& field) const {
    return mx_ == field.MX()
        && my_ == field.MY()
        && mz_ == field.MZ()
        && ng_ == field.NG()
        && totalSize_ == field.TotalSize()
        && (int)values_.size() == totalSize_;
}

} // namespace SF
