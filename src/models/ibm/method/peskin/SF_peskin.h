#pragma once

/// @file SF_peskin.h
/// @brief Peskin 原始显式浸没边界法方法族入口。
///
/// 当前实现 Bhalla Algorithm 1 的 prescribed-solid 特化：上一时间层表面
/// 乘子经规则化传播进入本步动量，新的乘子由速度误差显式更新。

#include "method/SF_method.h"
