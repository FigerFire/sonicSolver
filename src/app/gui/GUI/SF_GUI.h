#pragma once

/// @file SF_GUI.h
/// @brief GUI 控件层调度入口。
///
/// 其他模块只需要 include 这个文件，就能拿到 GUI 层主要控件。
/// 目前包含主窗口 MainWindow 和通用 true/false 开关 SwitchButton。
///
/// 实现文件拆分：
/// - SF_GUI.cpp：主窗口总装配、菜单、信号连接、全局样式。
/// - SF_mainWindowTree.cpp：左侧算例设置树、右键菜单、配置项编辑。
/// - SF_mainWindowLayout.cpp：右侧工作区、中间视图、底部输出/终端/日志。

#include "SF_mainWindow.h"
#include "SF_switchButton.h"
