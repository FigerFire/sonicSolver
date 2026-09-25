/// @file SF_flux.cpp
/// @brief canonical 面通量容器与访问实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#include "SF_flux.h"

#include <algorithm>

namespace SF {

void FluxField::setup(std::size_t faceCount, int nVar) {
    faceCount_ = faceCount;
    nVar_ = nVar;
    values_.assign(faceCount_ * (std::size_t)nVar_, 0.0);
}

void FluxField::clear() {
    std::fill(values_.begin(), values_.end(), 0.0);
}

double& FluxField::operator()(std::size_t faceId, int var) {
    return values_[index(faceId, var)];
}

double FluxField::operator()(std::size_t faceId, int var) const {
    return values_[index(faceId, var)];
}

std::size_t FluxField::index(std::size_t faceId, int var) const {
    return faceId * (std::size_t)nVar_ + (std::size_t)var;
}

} // namespace SF
