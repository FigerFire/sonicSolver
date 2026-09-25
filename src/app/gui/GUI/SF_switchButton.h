/// @file SF_switchButton.h
/// @brief GUI 主窗口布局、控件或视图切换实现。

#pragma once

#include <QAbstractButton>

namespace SF::GUI {

/// @brief GUI 层使用的 Apple 风格布尔开关按钮。
///
/// 它继承 QAbstractButton，所以也有 checked/toggled 信号。
/// 在 MainWindow::rebuildTree() 里，sf 文件中的 true/false 会显示成这个开关。
class SwitchButton : public QAbstractButton {
    Q_OBJECT

public:
    /// @brief 创建开关，并设置成可点击、可聚焦、固定大小。
    explicit SwitchButton(QWidget* parent = nullptr);
    /// @brief 告诉 Qt 这个控件默认需要多大空间。
    QSize sizeHint() const override;

protected:
    /// @brief 自己画开关轨道和圆形滑块。
    void paintEvent(QPaintEvent* event) override;
};

} // namespace SF::GUI
