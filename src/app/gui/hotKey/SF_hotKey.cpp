/// @file SF_hotKey.cpp
/// @brief GUI 快捷键注册与分派实现。

#include "SF_hotKey.h"

#include <QAbstractSpinBox>
#include <QApplication>
#include <QComboBox>
#include <QEvent>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QTextEdit>
#include <QWidget>

#include <utility>

namespace SF::GUI {
namespace {

/// @brief 判断当前焦点是否位于会接收文字输入的控件内部。
bool isEditableWidget(QWidget* widget)
{
    for (QWidget* current = widget; current; current = current->parentWidget()) {
        if (qobject_cast<QLineEdit*>(current)
            || qobject_cast<QTextEdit*>(current)
            || qobject_cast<QPlainTextEdit*>(current)
            || qobject_cast<QAbstractSpinBox*>(current)
            || qobject_cast<QComboBox*>(current)) {
            return true;
        }
    }
    return false;
}

} // namespace

ModuleHotKeyRouter::ModuleHotKeyRouter(QWidget* owner,
                                       QVector<ModuleHotKey> hotKeys,
                                       ModuleHotKeyHandler handler)
    : QObject(owner),
      owner_(owner),
      hotKeys_(std::move(hotKeys)),
      handler_(std::move(handler))
{
    qApp->installEventFilter(this);
}

ModuleHotKeyRouter::~ModuleHotKeyRouter()
{
    if (qApp) qApp->removeEventFilter(this);
}

bool ModuleHotKeyRouter::eventFilter(QObject* watched, QEvent* event)
{
    Q_UNUSED(watched);
    if (event->type() != QEvent::KeyPress || !owner_ || !handler_) {
        return false;
    }

    auto* keyEvent = static_cast<QKeyEvent*>(event);
    if (keyEvent->isAutoRepeat()) return false;
    if (keyEvent->modifiers()
        & (Qt::ControlModifier | Qt::AltModifier | Qt::MetaModifier)) {
        return false;
    }

    bool found = false;
    const ProjectModuleKind moduleKind =
        moduleForKey(keyEvent->key(), &found);
    if (!found || shouldIgnoreFocus()) return false;

    handler_(moduleKind);
    keyEvent->accept();
    return true;
}

bool ModuleHotKeyRouter::shouldIgnoreFocus() const
{
    if (!owner_->isActiveWindow()) return true;
    if (QApplication::activeModalWidget()
        || QApplication::activePopupWidget()) {
        return true;
    }

    QWidget* focus = QApplication::focusWidget();
    if (!focus) return false;
    if (focus != owner_ && !owner_->isAncestorOf(focus)) return true;
    return isEditableWidget(focus);
}

ProjectModuleKind ModuleHotKeyRouter::moduleForKey(int key, bool* found) const
{
    for (const ModuleHotKey& hotKey : hotKeys_) {
        if (hotKey.key == key) {
            if (found) *found = true;
            return hotKey.moduleKind;
        }
    }
    if (found) *found = false;
    return ProjectModuleKind::Geometry;
}

ModuleHotKeyRouter* installDefaultModuleHotKeys(
    QWidget* owner,
    ModuleHotKeyHandler handler)
{
    return new ModuleHotKeyRouter(
        owner,
        {
            {Qt::Key_G, ProjectModuleKind::Geometry,
             QStringLiteral("G：切换到几何模块")},
            {Qt::Key_M, ProjectModuleKind::Multiphase,
             QStringLiteral("M：切换到多相流模块")},
            {Qt::Key_S, ProjectModuleKind::Solver,
             QStringLiteral("S：切换到求解器模块")}
        },
        std::move(handler));
}

} // namespace SF::GUI
