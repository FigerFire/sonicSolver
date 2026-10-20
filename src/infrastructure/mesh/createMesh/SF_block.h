/// @file SF_block.h
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include "SF_edges.h"

namespace SF {

/// TFI (Transfinite Interpolation) 曲面+体插值
/// @param xi, eta, zeta  逻辑坐标 [0,1]
/// @param v              8 个顶点
/// @param edges          12 条边 (OOP 风格指针)，可为 nullptr (退化直线)
Vector3 calculateTFI(double xi, double eta, double zeta,
                     const Vector3 v[8], Edge* edges[12]);

} // namespace SF
