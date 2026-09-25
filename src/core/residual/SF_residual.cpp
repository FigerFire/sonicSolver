/// @file SF_residual.cpp
/// @brief GlobalDof/局部残差容器与累加实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#include "SF_residual.h"

#include <algorithm>

namespace SF {

void Residual::setup(int mx, int my, int mz, int nVar) {
    mx_ = mx;
    my_ = my;
    mz_ = mz;
    nVar_ = nVar;
    const std::size_t total =
        (std::size_t)mx_ * (std::size_t)my_ * (std::size_t)mz_
        * (std::size_t)nVar_;
    resX_.assign(total, 0.0);
    resY_.assign(total, 0.0);
    resZ_.assign(total, 0.0);
    source_.assign(total, 0.0);
    local_.assign(total, 0.0);
    global_.assign(total, 0.0);
    globalMask_.assign(
        (std::size_t)mx_ * (std::size_t)my_ * (std::size_t)mz_, 0);
}

void Residual::clear() {
    std::fill(resX_.begin(), resX_.end(), 0.0);
    std::fill(resY_.begin(), resY_.end(), 0.0);
    std::fill(resZ_.begin(), resZ_.end(), 0.0);
    std::fill(local_.begin(), local_.end(), 0.0);
    clearGlobal();
    clearSource();
}

void Residual::clearSource() {
    std::fill(source_.begin(), source_.end(), 0.0);
}

double& Residual::x(int i, int j, int k, int v) {
    return resX_[index(i, j, k, v)];
}

double Residual::x(int i, int j, int k, int v) const {
    return resX_[index(i, j, k, v)];
}

double& Residual::y(int i, int j, int k, int v) {
    return resY_[index(i, j, k, v)];
}

double Residual::y(int i, int j, int k, int v) const {
    return resY_[index(i, j, k, v)];
}

double& Residual::z(int i, int j, int k, int v) {
    return resZ_[index(i, j, k, v)];
}

double Residual::z(int i, int j, int k, int v) const {
    return resZ_[index(i, j, k, v)];
}

double& Residual::source(int i, int j, int k, int v) {
    return source_[index(i, j, k, v)];
}

double Residual::source(int i, int j, int k, int v) const {
    return source_[index(i, j, k, v)];
}

void Residual::stageLocal(
        int i, int j, int k, int v, double value) {
    local_[index(i, j, k, v)] = value;
}

double Residual::local(int i, int j, int k, int v) const {
    return local_[index(i, j, k, v)];
}

void Residual::setGlobal(
        int i, int j, int k, int v, double value) {
    global_[index(i, j, k, v)] = value;
    globalMask_[pointIndex(i, j, k)] = 1;
}

bool Residual::hasGlobal(int i, int j, int k) const {
    return globalMask_[pointIndex(i, j, k)] != 0;
}

double Residual::global(int i, int j, int k, int v) const {
    return global_[index(i, j, k, v)];
}

void Residual::clearGlobal() {
    std::fill(global_.begin(), global_.end(), 0.0);
    std::fill(globalMask_.begin(), globalMask_.end(), 0);
}

std::size_t Residual::index(int i, int j, int k, int v) const {
    return (((std::size_t)k * (std::size_t)my_ + (std::size_t)j)
            * (std::size_t)mx_ + (std::size_t)i)
            * (std::size_t)nVar_ + (std::size_t)v;
}

std::size_t Residual::pointIndex(int i, int j, int k) const {
    return ((std::size_t)k * (std::size_t)my_ + (std::size_t)j)
        * (std::size_t)mx_ + (std::size_t)i;
}

} // namespace SF
