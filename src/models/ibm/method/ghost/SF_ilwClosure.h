/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#pragma once

/// @file SF_ilwClosure.h
/// @brief IBM ghost-cell 的 ILW 几何前处理与高阶闭合接口。

#include "SF_field.h"
#include "SF_ibmConfig.h"
#include "SF_ilwStorage.h"
#include "models/ibm/topology/SF_ibmTopology.h"
#include "methods/numerics/structured/SF_structured.h"

namespace SF {
namespace IBM {
namespace GhostILW {

/// @brief IBM/ILW几何前处理统计。
struct PreprocessStats {
    int candidates = 0;
    int ghostCandidates = 0;
    int solidCandidates = 0;
    int built = 0;
    int reusedPlans = 0;
    int directBuiltPlans = 0;
    int insufficientFluidSamples = 0;
    int curvatureFailures = 0;
    int fitFailures = 0;
    int maxOrderLessThanTwo = 0;
    int maxOrderLessThanFour = 0;
    int maxOrderLessThanSix = 0;
    int maxOrderLessThanEight = 0;
    double planGeometrySeconds = 0.0;
    double planSampleSeconds = 0.0;
    double planSpacingSeconds = 0.0;
    double planProbeSeconds = 0.0;
    double planWallFitSeconds = 0.0;
    double planCurvatureSeconds = 0.0;
    double planHigherFitSeconds = 0.0;
};

/// @brief 在IBM setup阶段预计算ILW局部几何矩阵和导数提取权重。
///
/// 该函数只使用几何、单元拓扑和已缓存的样本索引，不读取流场状态；
/// 时间推进时ILW只需按预计算权重代入当前流场变量。
///
/// @param field 结构网格场，只读取IBM拓扑与单元分类。
/// @param storage IBM模块拥有的ILW样本与计划存储。
/// @param geometry IBM拓扑几何（ghost层、壁面/镜像点、法向、壁面速度）。
/// @param requestedMaxOrder 希望预计算的最高Taylor阶；ILW3/5/7/9对应2/4/6/8。
/// @param minimizeWallFitAmplification 是否在增长模板中选择权重放大最小的同阶壁面拟合。
/// @return 前处理统计信息。
PreprocessStats preprocessIBMGeometry(const Field& field,
                                      Storage& storage,
                                      const IBMGeometry& geometry,
                                      const IBMRuntimeConfig& config,
                                      int requestedMaxOrder = 4,
                                      bool minimizeWallFitAmplification = false);

/// @brief 构造一阶Euler滑移壁面ghost状态（静止壁面）。
/// @param qFluid 流体侧守恒状态 `[rho,rho*u,rho*v,rho*w,E]`。
/// @param wallNormal 从固体指向流体的壁面单位法向。
/// @brief 为一个IBM ghost/solid单元构造高阶ILW守恒状态。
/// @param field 结构网格场。
/// @param storage IBM模块拥有的ILW样本与计划存储。
/// @param geometry IBM拓扑几何（ghost层、壁面/镜像点、法向、壁面速度）。
/// @param ibmI IBM ghost/solid单元i索引。
/// @param ibmJ IBM ghost/solid单元j索引。
/// @param ibmK IBM ghost/solid单元k索引。
/// @param d 当前计算方向；张量ILW路径只用于保持调用点方向语义。
/// @param requestedOrder ILW Taylor最高阶数；ILW3/5/7/9对应2/4/6/8。
/// @param config application 层冻结的 IBM/ILW 与热力学配置。
/// @param qGhost 输出ghost守恒状态。
/// @return 成功构造物理有效ghost状态时返回true。
bool buildGhostStateForCell(const Field& field,
                            const Storage& storage,
                            const IBMGeometry& geometry,
                            int ibmI, int ibmJ, int ibmK,
                            Math::Dir d,
                            int requestedOrder,
                            const IBMRuntimeConfig& config,
                            double qGhost[5]);

} // namespace GhostILW
} // namespace IBM
} // namespace SF
