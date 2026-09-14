/// @file SF_boundary.cpp
/// @brief 物理边界共用的镜像索引与度量复制工具。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_boundary.h"

namespace SF {
namespace Boundary {

void copyMetrics(Field& field, int i, int j, int k,
                 int ii, int jj, int kk) {
    field.XiX(i, j, k) = field.XiX(ii, jj, kk);
    field.XiY(i, j, k) = field.XiY(ii, jj, kk);
    field.XiZ(i, j, k) = field.XiZ(ii, jj, kk);
    field.EtX(i, j, k) = field.EtX(ii, jj, kk);
    field.EtY(i, j, k) = field.EtY(ii, jj, kk);
    field.EtZ(i, j, k) = field.EtZ(ii, jj, kk);
    field.ZeX(i, j, k) = field.ZeX(ii, jj, kk);
    field.ZeY(i, j, k) = field.ZeY(ii, jj, kk);
    field.ZeZ(i, j, k) = field.ZeZ(ii, jj, kk);
    field.Jac(i, j, k) = field.Jac(ii, jj, kk);
}

void getMirrorIJK(int i, int j, int k,
                  int ng, int nx, int ny, int nz,
                  int& ii, int& jj, int& kk) {
    ii = i;
    jj = j;
    kk = k;
    if (i < ng) ii = 2 * ng - 1 - i;
    else if (i >= nx + ng) ii = 2 * (nx + ng) - 1 - i;
    if (j < ng) jj = 2 * ng - 1 - j;
    else if (j >= ny + ng) jj = 2 * (ny + ng) - 1 - j;
    if (k < ng) kk = 2 * ng - 1 - k;
    else if (k >= nz + ng) kk = 2 * (nz + ng) - 1 - k;
}

} // namespace Boundary
} // namespace SF
