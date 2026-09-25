/// @file SF_switchButton.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_switchButton.h"

#include <QPainter>

namespace SF::GUI {

SwitchButton::SwitchButton(QWidget* parent)
    : QAbstractButton(parent)
{
    // QAbstractButton 默认可以按，但不一定是“开关”。
    // setCheckable(true) 后，它才有 checked/unchecked 两种状态。
    setCheckable(true);
    // 鼠标悬停显示手形，告诉用户这是可点击控件。
    setCursor(Qt::PointingHandCursor);
    // 允许键盘 Tab 聚焦，方便以后做无鼠标操作。
    setFocusPolicy(Qt::StrongFocus);
    // 用 sizeHint() 返回的尺寸固定控件大小。
    setFixedSize(sizeHint());
}

QSize SwitchButton::sizeHint() const
{
    // 左侧配置树 true/false 开关的默认大小。
    return {42, 24};
}

void SwitchButton::paintEvent(QPaintEvent*)
{
    // paintEvent 是 Qt 自绘控件的核心：控件需要重画时 Qt 会调用它。
    QPainter painter(this);
    // 抗锯齿让圆角和圆形滑块更平滑。
    painter.setRenderHint(QPainter::Antialiasing);

    // track 是外层键槽。rect() 是整个控件区域，adjusted() 留出 1px 边距。
    const QRectF track = rect().adjusted(1.0, 1.0, -1.0, -1.0);
    painter.setPen(Qt::NoPen);
    // isChecked() 来自 QAbstractButton。true 用绿色，false 用灰色。
    painter.setBrush(isChecked() ? QColor(113, 133, 111)
                                 : QColor(181, 171, 154));
    painter.drawRoundedRect(track, track.height() / 2.0, track.height() / 2.0);

    // knob 是白色圆形滑块。checked 时靠右，否则靠左。
    const qreal diameter = track.height() - 4.0;
    const qreal x = isChecked()
        ? track.right() - diameter - 2.0
        : track.left() + 2.0;
    painter.setBrush(QColor(255, 252, 246));
    painter.drawEllipse(QRectF(x, track.top() + 2.0, diameter, diameter));
}

} // namespace SF::GUI
