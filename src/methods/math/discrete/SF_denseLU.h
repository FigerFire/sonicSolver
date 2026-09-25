#pragma once
/// @file SF_denseLU.h
/// @brief 小型稠密矩阵的带部分主元 LU；奇异时报告错误，不修改主元。
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::Math {
class DenseLU {
public:
    void factor(std::vector<double> values,int n) {
        if(n<=0 || values.size()!=static_cast<size_t>(n*n))
            throw std::runtime_error("DenseLU invalid dimensions.");
        n_=n; lu_=std::move(values); pivots_.resize(n);
        for(double v:lu_)if(!std::isfinite(v))
            throw std::runtime_error("DenseLU non-finite matrix.");
        for(int k=0;k<n_;++k) {
            int pivot=k;
            for(int i=k+1;i<n_;++i)
                if(std::abs(lu_[i*n_+k])>std::abs(lu_[pivot*n_+k]))pivot=i;
            if(lu_[pivot*n_+k]==0.0)
                throw std::runtime_error("DenseLU singular pivot at "+std::to_string(k));
            pivots_[k]=pivot;
            for(int j=0;j<n_;++j)std::swap(lu_[k*n_+j],lu_[pivot*n_+j]);
            for(int i=k+1;i<n_;++i) {
                lu_[i*n_+k]/=lu_[k*n_+k];
                for(int j=k+1;j<n_;++j)
                    lu_[i*n_+j]-=lu_[i*n_+k]*lu_[k*n_+j];
            }
        }
    }
    std::vector<double> solve(std::vector<double> x) const {
        if(n_<=0 || x.size()!=static_cast<size_t>(n_))
            throw std::runtime_error("DenseLU invalid right-hand side.");
        for(int k=0;k<n_;++k)std::swap(x[k],x[pivots_[k]]);
        for(int i=0;i<n_;++i)
            for(int j=0;j<i;++j)x[i]-=lu_[i*n_+j]*x[j];
        for(int i=n_-1;i>=0;--i) {
            for(int j=i+1;j<n_;++j)x[i]-=lu_[i*n_+j]*x[j];
            x[i]/=lu_[i*n_+i];
            if(!std::isfinite(x[i]))throw std::runtime_error("DenseLU non-finite solution.");
        }
        return x;
    }
private:
    int n_=0;
    std::vector<double> lu_;
    std::vector<int> pivots_;
};
} // namespace SF::Math
