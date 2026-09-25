/// @file SF_viewModeSwitch.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_viewModeSwitch.h"

#include <QPainter>

namespace SF::GUI {

ViewModeSwitch::ViewModeSwitch(QWidget* parent)
    : QAbstractButton(parent)
{
    // 和 SwitchButton 一样，这也是一个 checked/unchecked 双态按钮。
    // unchecked 表示 GEO，checked 表示 MESH。
    setCheckable(true);
    setCursor(Qt::PointingHandCursor);
    setFocusPolicy(Qt::StrongFocus);
    setFixedSize(sizeHint());
}

QSize ViewModeSwitch::sizeHint() const
{
    // 顶部 GEO/MESH 开关比普通 true/false 开关宽，因为要显示文字。
    return {90, 38};
}

void ViewModeSwitch::paintEvent(QPaintEvent*)
{
    // 这个控件完全由 QPainter 自绘，不依赖 Qt 默认按钮外观。
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // 外层键槽。
    const QRectF track = rect().adjusted(1.0, 1.0, -1.0, -1.0);
    painter.setPen(QPen(QColor(176, 165, 146, 170), 1.0));
    painter.setBrush(QColor(255, 252, 246, 205));
    painter.drawRoundedRect(track, track.height() / 2.0,
                            track.height() / 2.0);

    // 圆形滑块位置：checked 时在右侧，unchecked 时在左侧。
    const qreal diameter = track.height() - 8.0;
    const qreal knobX = isChecked()
        ? track.right() - diameter - 4.0
        : track.left() + 4.0;
    const QRectF knob(knobX, track.top() + 4.0, diameter, diameter);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(92, 112, 92));
    painter.drawEllipse(knob);

    // 开关文字。为了表达“另一侧可切换到什么”，
    // checked(MESH) 时左侧显示 MESH，unchecked(GEO) 时右侧显示 GEO。
    QFont labelFont = font();
    labelFont.setBold(true);
    labelFont.setPointSizeF(9.0);
    painter.setFont(labelFont);
    painter.setPen(QColor(84, 77, 66));
    if (isChecked()) {
        painter.drawText(QRectF(track.left() + 10.0, track.top(),
                                track.width() - diameter - 12.0,
                                track.height()),
                         Qt::AlignCenter, QStringLiteral("MESH"));
    } else {
        painter.drawText(QRectF(knob.right() + 2.0, track.top(),
                                track.right() - knob.right() - 8.0,
                                track.height()),
                         Qt::AlignCenter, QStringLiteral("GEO"));
    }
}

} // namespace SF::GUI
