/// @file SF_physics.h
/// @brief 物理模型注册、选择与统一调度实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "gravity/SF_gravity.h"
#include "mrf/SF_mrf.h"
#include "multiphase/SF_multiphase.h"
#include "EOS/SF_EOS.h"
#include "EOS/SF_perfectGasEOS.h"
#include "EOS/SF_stiffenedGasEOS.h"
#include "EOS/SF_taitEOS.h"
#include "EOS/SF_pengRobinsonEOS.h"
#include "EOS/SF_tabulatedEOS.h"
#include "SF_fluidStateModel.h"
#include "fluidStateModel/SF_singleFluidStateModel.h"
#include "fluidStateModel/SF_homogeneousMultiphaseStateModel.h"
#include "fluidStateModel/SF_factory.h"
