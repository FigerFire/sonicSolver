#pragma once

/// @file SF_interactorStyle.h
/// @brief ParaView 风格 VTK 鼠标交互。

#ifdef SONIC_GUI_WITH_VTK

#include <vtkInteractorStyle.h>
#include <vtkSmartPointer.h>

namespace SF::GUI {

/// @brief 创建左键旋转、右键平移的 VTK 相机交互样式。
vtkSmartPointer<vtkInteractorStyle> createParaViewInteractorStyle();

} // namespace SF::GUI

#endif
