/// @file SF_init.h
/// @brief 初值设置与 legacy 配置迁移实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include "SF_field.h"
#include "SF_configTypes.h"

namespace SF {
    namespace Init {
        /// @brief 根据显式值对象初始化 legacy 五变量守恒场。
        void setupAll(SF::Field& field,
                      const FDM::InitialConditionConfig& config,
                      double idealGasGamma,
                      double idealGasConstant);

    }
}
