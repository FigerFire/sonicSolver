/// @file SF_viewModeSwitch.h
/// @brief GUI 主窗口布局、控件或视图切换实现。

#pragma once

#include <QAbstractButton>

namespace SF::GUI {

/// @brief GEO/MESH 双态键槽开关；选中状态表示网格与结果视图。
///
/// unchecked：显示 GEO，MainWindow 切到 GeometryEditor。
/// checked：显示 MESH，MainWindow 切到 VTKView。
class ViewModeSwitch final : public QAbstractButton {
    Q_OBJECT

public:
    /// @brief 创建 GEO/MESH 开关。
    explicit ViewModeSwitch(QWidget* parent = nullptr);
    /// @brief 告诉 Qt 这个控件默认需要多大空间。
    QSize sizeHint() const override;

protected:
    /// @brief 自己画键槽、滑块和 GEO/MESH 文字。
    void paintEvent(QPaintEvent* event) override;
};

} // namespace SF::GUI
