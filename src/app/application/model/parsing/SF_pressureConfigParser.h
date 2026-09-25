#pragma once

#include "SF_parserCommon.h"
#include "core/config/SF_configTypes.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {
// 线性求解与压力工作流的校验已移入 core/config/types 的 typed 值对象：
// 它们是值对象自校验，不属于解析层，solver 也不能反向依赖 application。
} // namespace SF::FDM
