/// @file SF_riemann.h
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <cmath>

namespace SF {
namespace Riemann {

/// Steger-Warming 正负通量分裂 (物理空间)
void splitStegerWarming(const double* Q, const double* n,
                        double gamma, double* f_pos, double* f_neg);

/// Lax-Friedrichs 正负通量分裂 (物理空间)
/// f^± = 0.5 * (F_phys(Q)·n̂ ± λ_max * Q)
/// λ_max = |u·n̂| + a
void splitLaxFriedrichs(const double* Q, const double* n,
                        double gamma, double* f_pos, double* f_neg);

/// 构造 Euler 方程的左右特征矩阵
/// Q_avg[5]: 界面平均守恒变量
/// n[3]: 界面单位法向量
/// L[25]: 左特征矩阵 (行优先)
/// R[25]: 右特征矩阵 (行优先)
void buildCharacteristicMatrix(const double* Q_avg, const double* n,
                               double gamma, double* L, double* R);

/// 5×5 矩阵 × 5×1 矢量
inline void gemv5(const double* M, const double* V, double* R) {
    for (int i = 0; i < 5; ++i)
        R[i] = M[i*5+0]*V[0] + M[i*5+1]*V[1] + M[i*5+2]*V[2]
             + M[i*5+3]*V[3] + M[i*5+4]*V[4];
}

} // namespace Riemann
} // namespace SF
