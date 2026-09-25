/// @file SF_editor.h
/// @brief GUI 几何编辑与参数面板实现。

#pragma once

#include "SF_project.h"

#include <QWidget>

class QLabel;

namespace SF::GUI {

class GeometryCanvas;

/// @brief 用线段和三点圆弧构造单块二维截面的几何编辑器。
class GeometryEditor : public QWidget {
    Q_OBJECT

public:
    explicit GeometryEditor(QWidget* parent = nullptr);

    /// @brief 从 blockMeshDict 或兼容几何配置文件加载并显示几何轮廓。
    /// @param path `blockMeshDict` 或兼容几何配置文件。
    /// @param error 失败时返回诊断信息。
    /// @return 是否成功解析并显示。
    bool loadFile(const QString& path, QString* error = nullptr);

    /// @brief 切到选择/占位工具；用于顶部几何工具栏。
    void setSelectMode();
    /// @brief 切到线段绘制工具。
    void setLineMode();
    /// @brief 切到圆弧绘制工具。
    void setArcMode();

signals:
    /// @brief 请求将当前草图转换并保存为网格参数。
    /// @param sketch 当前二维草图。
    void saveRequested(const GeometrySketch& sketch);

private:
    GeometryCanvas* canvas_ = nullptr;
    QLabel* statusLabel_ = nullptr;
};

} // namespace SF::GUI
