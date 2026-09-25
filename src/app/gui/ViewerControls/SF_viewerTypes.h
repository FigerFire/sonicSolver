#pragma once

/// @file SF_viewerTypes.h
/// @brief 3D 预览控件之间共享的轻量枚举。

namespace SF::GUI {

/// @brief 顶部工具栏当前服务的主工作区。
enum class ViewerWorkspaceMode {
    Geometry,
    Mesh
};

/// @brief 切片法向轴；坐标值沿该轴解释。
enum class ViewerSliceAxis {
    X,
    Y,
    Z
};

/// @brief VTK 预览的显示模式。
enum class ViewerDisplayMode {
    /// @brief 按当前物理标量渲染外表面。
    PhysicalSurface,
    /// @brief 只显示网格线框。
    Mesh
};

/// @brief 六个轴向正视角。
enum class ViewerCameraDirection {
    PositiveX,
    NegativeX,
    PositiveY,
    NegativeY,
    PositiveZ,
    NegativeZ
};

/// @brief 几何草图工具栏中的占位工具。
enum class ViewerGeometryTool {
    Vertex,
    Line,
    CenterArc,
    LengthDimension,
    AngleDimension,
    TrimSegment
};

} // namespace SF::GUI
