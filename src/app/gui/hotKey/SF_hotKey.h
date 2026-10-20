#pragma once

/// @file SF_hotKey.h
/// @brief sonicGui 热键模块调度入口。
///
/// 热键层只负责把键盘输入翻译成 GUI 模块选择请求，不直接操作树节点、
/// 不读写 case 文件，也不启动求解器。真正的模块展开与视图切换继续由
/// MainWindow::focusModule(...) 完成。

#include "SF_project.h"

#include <QObject>
#include <QVector>

#include <functional>

class QWidget;

namespace SF::GUI {

/// @brief 一个模块热键绑定。
struct ModuleHotKey {
    /// @brief Qt::Key_* 键值，例如 Qt::Key_G。
    int key = 0;
    /// @brief 触发后要切换到的 GUI 业务模块。
    ProjectModuleKind moduleKind = ProjectModuleKind::Geometry;
    /// @brief 给代码阅读者看的中文说明，不参与显示。
    QString description;
};

/// @brief 模块热键触发后的回调。
using ModuleHotKeyHandler = std::function<void(ProjectModuleKind kind)>;

/// @brief 安装在 QApplication 上的模块热键路由器。
///
/// 它会在主窗口激活且焦点不在可编辑控件中时处理单键热键。
/// 这样按 G/M/S 可以快速切换模块，同时不会影响终端输入框或配置编辑框打字。
class ModuleHotKeyRouter final : public QObject {
public:
    /// @brief 创建热键路由器。
    /// @param owner 热键所属主窗口；只有该窗口激活时才响应。
    /// @param hotKeys 要注册的模块热键列表。
    /// @param handler 命中热键后的模块切换回调。
    explicit ModuleHotKeyRouter(QWidget* owner,
                                QVector<ModuleHotKey> hotKeys,
                                ModuleHotKeyHandler handler);
    /// @brief 析构时从 QApplication 移除事件过滤器。
    ~ModuleHotKeyRouter() override;

protected:
    /// @brief 过滤键盘事件并触发模块切换。
    bool eventFilter(QObject* watched, QEvent* event) override;

private:
    bool shouldIgnoreFocus() const;
    ProjectModuleKind moduleForKey(int key, bool* found) const;

    QWidget* owner_ = nullptr;
    QVector<ModuleHotKey> hotKeys_;
    ModuleHotKeyHandler handler_;
};

/// @brief 注册 sonicGui 默认模块热键。
///
/// 当前只注册：
/// - G：几何
/// - M：网格
/// - S：求解器
///
/// IBM 和 MRF 暂不注册热键。
ModuleHotKeyRouter* installDefaultModuleHotKeys(
    QWidget* owner,
    ModuleHotKeyHandler handler);

} // namespace SF::GUI
