/// @file SF_GUI.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_mainWindow.h"

#include "SF_Geometry.h"
#include "SF_hotKey.h"
#include "SF_toolbar.h"
#include "SF_viewModeSwitch.h"
#include "SF_VTKView.h"

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QIcon>
#include <QKeySequence>
#include <QMenu>
#include <QMenuBar>
#include <QPainter>
#include <QProxyStyle>
#include <QStyle>
#include <QStyleOption>
#include <QStackedWidget>
#include <QTabWidget>
#include <QTreeWidget>

// GUI 控件层总装配文件：资源初始化、模块图标映射、菜单、信号连接和全局样式。
// 左侧配置树实现见 SF_mainWindowTree.cpp。
// 右侧工作区和底部输出实现见 SF_mainWindowLayout.cpp。

/// @brief 初始化 Qt Resource 中的 GUI 图标。
///
/// Qt 的 .qrc 资源在静态库里时，有时需要显式 Q_INIT_RESOURCE，
/// 否则运行时 QIcon(":/sonic/icons/...") 可能找不到图片。
static void initializeGuiResources()
{
    Q_INIT_RESOURCE(SF_resources);
}

namespace SF::GUI {
namespace {

/// @brief 把 QTreeWidget 默认的小三角展开符号改成 STAR-CCM 类似的 +/- 方框。
///
/// QProxyStyle 是 Qt 的“样式代理”。这里只拦截 PE_IndicatorBranch，
/// 也就是树节点左边的展开/收起符号；其他控件仍然交给原始 style 绘制。
class PlusMinusTreeStyle final : public QProxyStyle {
public:
    using QProxyStyle::QProxyStyle;

    void drawPrimitive(PrimitiveElement element,
                       const QStyleOption* option,
                       QPainter* painter,
                       const QWidget* widget = nullptr) const override
    {
        if (element != PE_IndicatorBranch
            || !(option->state & State_Children)) {
            QProxyStyle::drawPrimitive(element, option, painter, widget);
            return;
        }

        const QRect box(option->rect.center().x() - 6,
                        option->rect.center().y() - 6, 12, 12);
        painter->save();
        painter->setRenderHint(QPainter::Antialiasing);
        painter->setPen(QPen(QColor(133, 124, 109), 1.2));
        painter->setBrush(QColor(248, 243, 234));
        painter->drawRoundedRect(box, 2.5, 2.5);
        painter->drawLine(box.left() + 3, box.center().y(),
                          box.right() - 3, box.center().y());
        if (!(option->state & State_Open)) {
            painter->drawLine(box.center().x(), box.top() + 3,
                              box.center().x(), box.bottom() - 3);
        }
        painter->restore();
    }
};


} // namespace

QIcon MainWindow::moduleIcon(ProjectModuleKind kind)
{
    // GUI 图标都通过 src/app/gui/Resources/SF_resources.qrc 统一打包。
    // 这里把业务模块枚举映射到稳定的 :/sonic/icons/... 运行时路径。
    switch (kind) {
        case ProjectModuleKind::Geometry:
            return QIcon(QStringLiteral(":/sonic/icons/geometry.png"));
        case ProjectModuleKind::Mesh:
            return QIcon(QStringLiteral(":/sonic/icons/mesh.png"));
        case ProjectModuleKind::Solver:
            return QIcon(QStringLiteral(":/sonic/icons/solver.png"));
        case ProjectModuleKind::Turbulence:
            return QIcon(QStringLiteral(":/sonic/icons/turbulence.png"));
        case ProjectModuleKind::Multiphase:
        case ProjectModuleKind::PhaseChange:
            return QIcon(QStringLiteral(":/sonic/icons/turbulence.png"));
        case ProjectModuleKind::InitialCondition:
            return QIcon(QStringLiteral(":/sonic/icons/initial-condition.png"));
        case ProjectModuleKind::BoundaryCondition:
            return QIcon(QStringLiteral(":/sonic/icons/boundary-condition.png"));
        case ProjectModuleKind::IBM:
            return QIcon(QStringLiteral(":/sonic/icons/ibm.png"));
        case ProjectModuleKind::MRF:
            return QIcon(QStringLiteral(":/sonic/icons/mrf.png"));
        case ProjectModuleKind::Output:
            return QIcon(QStringLiteral(":/sonic/icons/output.png"));
        case ProjectModuleKind::Parallel:
            return QIcon(QStringLiteral(":/sonic/icons/parallel.png"));
    }
    return QIcon(QStringLiteral(":/sonic/icons/logo.png"));
}

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    // 构造顺序很重要：
    // 1. 先加载 qrc 图标资源。
    // 2. 创建控件和菜单。
    // 3. 套 stylesheet 和树展开样式。
    // 4. 最后 connect 信号，避免控件还没创建就连接。
    initializeGuiResources();
    setWindowIcon(QIcon(QStringLiteral(":/sonic/icons/logo.png")));
    buildUi();
    buildMenus();
    applyStyle();
    configTree_->setStyle(new PlusMinusTreeStyle(configTree_->style()));
    connectUi();
    installDefaultModuleHotKeys(this, [this](ProjectModuleKind kind) {
        focusModule(kind);
    });

    setWindowTitle(QStringLiteral("sonicGui"));
    resize(1440, 920);
    setMinimumSize(1060, 680);
#ifdef Q_OS_MACOS
    setUnifiedTitleAndToolBarOnMac(true);
#endif
}

void MainWindow::buildMenus()
{
    // 顶部菜单栏：大部分菜单 action 都不直接做事，只 emit 一个 signal。
    // 真正的业务动作在 Controller 层完成。
    QMenu* fileMenu = menuBar()->addMenu(QStringLiteral("文件"));

    // 文件 -> 设置工作文件夹...
    // 和左侧小文件夹按钮发同一个 workingDirectoryRequested()。
    QAction* workspaceAction =
        fileMenu->addAction(
            QIcon(QStringLiteral(":/sonic/icons/set-directory.png")),
            QStringLiteral("设置工作文件夹..."));
    workspaceAction->setShortcut(QKeySequence::Open);
    connect(workspaceAction, &QAction::triggered,
            this, &MainWindow::workingDirectoryRequested);

    // 文件 -> 保存全部
    // Controller 收到 saveAllRequested() 后写回全部 case 字典。
    QAction* saveAction = fileMenu->addAction(QStringLiteral("保存全部"));
    saveAction->setShortcut(QKeySequence::Save);
    connect(saveAction, &QAction::triggered,
            this, &MainWindow::saveAllRequested);

    fileMenu->addSeparator();
    // 文件 -> 导入 SFM 网格...
    QAction* sfmAction = fileMenu->addAction(
        moduleIcon(ProjectModuleKind::Mesh),
        QStringLiteral("导入 SFM 网格..."));
    connect(sfmAction, &QAction::triggered,
            this, &MainWindow::importSfmRequested);
    // 文件 -> 导入 IBM STL...
    QAction* stlAction = fileMenu->addAction(
        moduleIcon(ProjectModuleKind::IBM),
        QStringLiteral("导入 IBM STL..."));
    connect(stlAction, &QAction::triggered,
            this, &MainWindow::importStlRequested);
    // 文件 -> 打开 PVD/VTK 结果...
    QAction* vtkAction =
        fileMenu->addAction(
            moduleIcon(ProjectModuleKind::Mesh),
            QStringLiteral("打开 PVD/VTK 结果..."));
    connect(vtkAction, &QAction::triggered,
            this, &MainWindow::importVtkRequested);

    fileMenu->addSeparator();
    QAction* quitAction = fileMenu->addAction(QStringLiteral("退出"));
    quitAction->setShortcut(QKeySequence::Quit);
    connect(quitAction, &QAction::triggered, qApp, &QApplication::quit);

    QMenu* moduleMenu = menuBar()->addMenu(QStringLiteral("模块"));
    // 模块菜单和左侧树右键菜单共用这三个动作。
    addOptionalModuleActions(moduleMenu);

    QMenu* runMenu = menuBar()->addMenu(QStringLiteral("运行"));
    // 运行 -> 运行当前算例。previewOneStep=false 表示正常求解。
    QAction* runAction = runMenu->addAction(
        moduleIcon(ProjectModuleKind::Solver),
        QStringLiteral("运行当前算例"));
    runAction->setShortcut(QKeySequence(QStringLiteral("Ctrl+R")));
    connect(runAction, &QAction::triggered, this,
            [this] { emit runRequested(false); });
    // 运行 -> 计算一步并预览。previewOneStep=true 表示只做一步预览。
    QAction* previewAction = runMenu->addAction(
        moduleIcon(ProjectModuleKind::Solver),
        QStringLiteral("计算一步并预览"));
    previewAction->setShortcut(QKeySequence(QStringLiteral("Ctrl+Shift+R")));
    connect(previewAction, &QAction::triggered, this,
            [this] { emit runRequested(true); });
    // 运行 -> 停止。只有求解器正在跑时才由 setSolverRunning(true) 启用。
    stopAction_ = runMenu->addAction(QStringLiteral("停止"));
    stopAction_->setEnabled(false);
    connect(stopAction_, &QAction::triggered, this, &MainWindow::stopRequested);

    QMenu* outputMenu = menuBar()->addMenu(QStringLiteral("输出"));
    QAction* newOutputTaskAction = outputMenu->addAction(
        moduleIcon(ProjectModuleKind::Output), QStringLiteral("新建任务"));
    connect(newOutputTaskAction, &QAction::triggered,
            this, &MainWindow::requestCreateOutputTask);
    submitOutputTaskAction_ = outputMenu->addAction(
        moduleIcon(ProjectModuleKind::Solver), QStringLiteral("提交当前任务"));
    connect(submitOutputTaskAction_, &QAction::triggered, this, [this] {
        const QString taskName = currentOutputTaskName();
        if (!taskName.isEmpty()) emit outputTaskSubmitRequested(taskName);
    });
    QAction* residualAction = outputMenu->addAction(
        moduleIcon(ProjectModuleKind::Output), QStringLiteral("残差显示"));
    connect(residualAction, &QAction::triggered, this, [this] {
        const QString taskName = currentOutputTaskName();
        if (!taskName.isEmpty()) emit outputResidualRequested(taskName);
    });

    QMenu* viewMenu = menuBar()->addMenu(QStringLiteral("视图"));
    // 视图 -> 几何 / 网格与结果，和顶部 GEO/MESH 开关控制同一个 viewerStack_。
    QAction* geometryAction = viewMenu->addAction(QStringLiteral("几何"));
    connect(geometryAction, &QAction::triggered,
            this, &MainWindow::showGeometryView);
    QAction* resultAction = viewMenu->addAction(QStringLiteral("网格与结果"));
    connect(resultAction, &QAction::triggered,
            this, &MainWindow::showVisualizationView);
    QAction* bottomAction = viewMenu->addAction(QStringLiteral("底部面板"));
    bottomAction->setCheckable(true);
    bottomAction->setChecked(true);
    connect(bottomAction, &QAction::toggled, bottomTabs_, &QWidget::setVisible);

    // “已启用”菜单不是固定内容，updateEnabledMenu() 会动态填充所有 true 开关。
    enabledMenu_ = menuBar()->addMenu(QStringLiteral("已启用"));
    enabledMenu_->setEnabled(false);
}

void MainWindow::connectUi()
{
    // 左侧树右键 -> showModuleContextMenu()，根据点中的 section 决定能新增什么。
    connect(configTree_, &QTreeWidget::customContextMenuRequested,
            this, &MainWindow::showModuleContextMenu);
    // 左侧树双击 -> 直接编辑第二列属性值。
    connect(configTree_, &QTreeWidget::itemDoubleClicked, this,
            [this](QTreeWidgetItem* item) { editTreeItem(item); });
    connect(configTree_, &QTreeWidget::itemChanged,
            this, &MainWindow::commitTreeItemEdit);
    // 左侧树单击 -> 同步顶部模块选择器；如果进入几何/网格模块，也同步中间视图。
    connect(configTree_, &QTreeWidget::itemClicked, this,
            [this](QTreeWidgetItem* item) {
        if (!item) return;
        if (handleOutputTreeClick(item)) return;
        QTreeWidgetItem* top = item;
        while (top->parent()) top = top->parent();
        if (!top->data(0, ModuleRole).isValid()) return;
        const int moduleIndex = top->data(0, ModuleRole).toInt();
        if (moduleIndex < 0 || moduleIndex >= modules_.size()) return;
        syncModuleSelector(moduleIndex);
        const ProjectModuleKind kind = modules_.at(moduleIndex).kind;
        if (kind == ProjectModuleKind::Geometry) {
            showGeometryView();
        } else if (kind == ProjectModuleKind::Mesh) {
            showVisualizationView();
        }
    });
    // GEO/MESH 开关 -> 切换 viewerStack_ 当前页面。
    connect(viewModeSwitch_, &QAbstractButton::toggled,
            this, [this](bool showMesh) {
        if (!showMesh) {
            showGeometryView();
        } else {
            showVisualizationView();
        }
    });
    connect(viewerToolbar_, &ViewerToolbar::geometryToolChanged,
            this, [this](ViewerGeometryTool tool) {
        switch (tool) {
            case ViewerGeometryTool::Vertex:
                geometryEditor_->setSelectMode();
                appendLog(QStringLiteral("几何工具：顶点（界面占位）"));
                break;
            case ViewerGeometryTool::Line:
                geometryEditor_->setLineMode();
                appendLog(QStringLiteral("几何工具：直线"));
                break;
            case ViewerGeometryTool::CenterArc:
                geometryEditor_->setArcMode();
                appendLog(QStringLiteral("几何工具：以中心画圆弧（暂复用圆弧绘制）"));
                break;
            case ViewerGeometryTool::LengthDimension:
                geometryEditor_->setSelectMode();
                appendLog(QStringLiteral("几何工具：标识长度（界面占位）"));
                break;
            case ViewerGeometryTool::AngleDimension:
                geometryEditor_->setSelectMode();
                appendLog(QStringLiteral("几何工具：标识角度（界面占位）"));
                break;
            case ViewerGeometryTool::TrimSegment:
                geometryEditor_->setSelectMode();
                appendLog(QStringLiteral("几何工具：剪裁线段（界面占位）"));
                break;
        }
    });
    connect(viewerToolbar_, &ViewerToolbar::solidExtrudeRequested,
            this, [this](double length) {
        appendLog(QStringLiteral("实体操作：拉伸长度 %1（界面占位）")
                      .arg(length));
    });
    connect(viewerToolbar_, &ViewerToolbar::solidRevolveRequested,
            this, [this](double angle) {
        appendLog(QStringLiteral("实体操作：旋转扫描角度 %1（界面占位）")
                      .arg(angle));
    });
    // 高亮只跟随“当前选择”，鼠标划过不会改变模型显示；点击空白处清除选择即可复位。
    connect(configTree_, &QTreeWidget::itemSelectionChanged, this, [this] {
        const QList<QTreeWidgetItem*> selected = configTree_->selectedItems();
        handleTreeHover(selected.isEmpty() ? nullptr : selected.constFirst());
    });
    connect(this, &MainWindow::geometryHoverBlock,
            vtkView_, &VTKView::highlightBlock);
    connect(this, &MainWindow::geometryHoverPatch,
            vtkView_, &VTKView::highlightPatch);
    connect(this, &MainWindow::geometryHoverFace,
            vtkView_, &VTKView::highlightFace);
    connect(this, &MainWindow::geometryHoverClear,
            vtkView_, &VTKView::clearHighlight);
    // VTKView 发来的状态信息 -> 底部日志页。
    connect(vtkView_, &VTKView::statusMessage,
            this, &MainWindow::appendLog);
    // 几何编辑器要求保存 -> MainWindow 继续转发给 Controller。
    connect(geometryEditor_, &GeometryEditor::saveRequested,
            this, &MainWindow::geometrySaveRequested);
}

void MainWindow::addOptionalModuleActions(QMenu* menu)
{
    // 这三个动作只发 addModuleRequested(kind)，真正创建文件和修改 system/sonicDict
    // 在 Controller/Project 层处理。
    QAction* ibm = menu->addAction(
        moduleIcon(ProjectModuleKind::IBM), QStringLiteral("新增 IBM"));
    connect(ibm, &QAction::triggered, this, [this] {
        emit addModuleRequested(ProjectModuleKind::IBM);
    });
    QAction* mrf = menu->addAction(
        moduleIcon(ProjectModuleKind::MRF), QStringLiteral("新增 MRF"));
    connect(mrf, &QAction::triggered, this, [this] {
        emit addModuleRequested(ProjectModuleKind::MRF);
    });
    QAction* parallel = menu->addAction(
        moduleIcon(ProjectModuleKind::Parallel),
        QStringLiteral("新增 Parallel"));
    connect(parallel, &QAction::triggered, this, [this] {
        emit addModuleRequested(ProjectModuleKind::Parallel);
    });
    QAction* multiphase = menu->addAction(
        moduleIcon(ProjectModuleKind::Multiphase),
        QStringLiteral("创建多相流 (M)"));
    connect(multiphase, &QAction::triggered, this, [this] {
        emit addModuleRequested(ProjectModuleKind::Multiphase);
    });
}

void MainWindow::applyStyle()
{
    // 这里集中写 Qt stylesheet。
    // objectName 对应 buildUi() 里 setObjectName() 设置的名字，例如 #sidebar。
    // 如果你想改某个按钮/区域的颜色，先看它有没有 objectName。
    setWindowOpacity(1.0);
    qApp->setStyleSheet(QStringLiteral(R"(
        QMainWindow, QDialog {
            background: #e9e2d7;
            color: #3d382f;
        }
        QMenuBar {
            background: rgba(250, 247, 240, 242);
            color: #3d382f;
            border-bottom: 1px solid rgba(145, 132, 111, 90);
        }
        QMenuBar::item:selected, QMenu::item:selected {
            background: #71856f;
            color: white;
            border-radius: 7px;
        }
        QMenu {
            background: rgba(252, 249, 242, 248);
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 120);
            border-radius: 12px;
            padding: 7px;
        }
        QSplitter::handle { background: #cfc4b4; }
        QSplitter::handle:horizontal { width: 1px; }
        QSplitter::handle:vertical { height: 1px; }
        #rightWorkspaceSplitter::handle:vertical {
            height: 5px;
            background: #d2c8b8;
            border-top: 1px solid rgba(145, 132, 111, 90);
            border-bottom: 1px solid rgba(145, 132, 111, 90);
        }
        #sidebar {
            background: rgba(247, 242, 233, 225);
            border-right: 1px solid rgba(145, 132, 111, 90);
        }
        #sectionTitle {
            color: #3d382f;
            font-size: 13px;
            font-weight: 700;
            letter-spacing: 1px;
        }
        #panelTitle { color: #3d382f; font-weight: 600; }
        #mutedLabel { color: #817767; }
        #workspaceButton {
            min-width: 32px;
            max-width: 32px;
            min-height: 32px;
            max-height: 32px;
            padding: 4px;
            border-radius: 10px;
        }
        QTreeWidget {
            background: transparent;
            color: #4a443a;
            border: none;
            outline: none;
        }
        QTreeWidget::item { min-height: 28px; border-radius: 6px; }
        QTreeWidget::item:hover { background: rgba(255, 255, 255, 100); }
        QTreeWidget::item:selected {
            background: #71856f;
            color: white;
        }
        QLineEdit, QPlainTextEdit {
            background: rgba(255, 252, 246, 205);
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 110);
            border-radius: 10px;
            padding: 7px;
            selection-background-color: #71856f;
        }
        #moduleSelector {
            background: rgba(255, 252, 246, 205);
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 110);
            border-radius: 11px;
            padding: 6px 12px;
            font-weight: 600;
        }
        #moduleSelector:hover {
            background: rgba(255, 255, 255, 235);
            border-color: rgba(112, 132, 108, 170);
        }
        #moduleSelector::menu-indicator {
            image: none;
            width: 0px;
        }
        QPushButton, QToolButton {
            background: rgba(255, 252, 246, 190);
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 105);
            border-radius: 9px;
            padding: 6px 11px;
        }
        QPushButton:hover, QToolButton:hover {
            background: rgba(255, 255, 255, 235);
        }
        QPushButton:pressed, QToolButton:checked {
            background: #71856f;
            color: white;
        }
        #primaryButton {
            background: #71856f;
            color: white;
            border-color: #71856f;
            font-weight: 600;
        }
        #viewerToolbar, #viewModeBar {
            background: #f8f3ea;
            border-bottom: 1px solid rgba(145, 132, 111, 90);
        }
        #viewerToolbar QToolButton {
            background: #fffaf2;
            min-width: 28px;
            max-width: 28px;
            min-height: 28px;
            max-height: 28px;
            padding: 2px;
            border-radius: 7px;
        }
        #viewerToolbar QToolButton:hover {
            background: #ffffff;
        }
        #viewerToolbar QToolButton:checked {
            background: #71856f;
            color: white;
        }
        #viewerToolbar QToolButton:disabled {
            background: rgba(255, 250, 242, 110);
            color: rgba(74, 68, 58, 95);
            border-color: rgba(145, 132, 111, 55);
        }
        #viewerToolbar #moduleSelector {
            min-width: 120px;
            max-width: 220px;
            min-height: 32px;
            max-height: 32px;
            padding: 4px 10px;
            border-radius: 8px;
            text-align: left;
        }
        #viewerScalarSelector {
            background: #fffaf2;
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 100);
            border-radius: 7px;
            padding: 2px 6px;
        }
        #viewerTimeLabel {
            color: #4a443a;
            font-weight: 600;
            padding-left: 4px;
        }
        #playbackRateLabel {
            color: #655d51;
            font-size: 11px;
            font-weight: 600;
        }
        #playbackRateSlider::groove:horizontal {
            height: 5px;
            background: #d5cbbb;
            border-radius: 2px;
        }
        #playbackRateSlider::sub-page:horizontal {
            background: #71856f;
            border-radius: 2px;
        }
        #playbackRateSlider::handle:horizontal {
            width: 13px;
            margin: -5px 0;
            background: #fffaf2;
            border: 1px solid #71856f;
            border-radius: 6px;
        }
        #slicePlanePanel {
            background: #fffaf2;
            border: 1px solid rgba(145, 132, 111, 120);
            border-radius: 8px;
        }
        #slicePlanePanel QLabel {
            color: #4a443a;
            font-weight: 600;
        }
        #slicePlanePanel QDoubleSpinBox {
            background: #ffffff;
            color: #3d382f;
            border: 1px solid rgba(145, 132, 111, 110);
            border-radius: 6px;
            padding: 1px 3px;
        }
        #slicePlanePanel QPushButton {
            padding: 1px 4px;
            border-radius: 6px;
        }
        #viewerPlaceholder {
            color: #817767;
            background: qradialgradient(
                cx:0.5, cy:0.42, radius:0.85,
                stop:0 rgba(255, 252, 245, 230),
                stop:1 rgba(220, 210, 194, 230));
        }
        QTabWidget::pane {
            border: none;
            background: rgba(239, 232, 220, 220);
        }
        #bottomOutputTabs {
            background: #eee7dc;
        }
        #bottomOutputTabs::pane {
            border-top: 1px solid rgba(145, 132, 111, 90);
            background: #f5efe5;
        }
        #bottomOutputTabs::tab-bar {
            alignment: left;
        }
        #bottomOutputTabs QTabBar {
            background: #eee7dc;
        }
        #bottomTerminalPage {
            background: #f5efe5;
        }
        #bottomConsole {
            background: #fffdf8;
            color: #3d382f;
            border: none;
            border-radius: 0px;
            padding: 7px;
            selection-background-color: #71856f;
        }
        #bottomCommandInput {
            background: #fffdf8;
            color: #3d382f;
            border: 1px solid #aa9f8f;
            border-radius: 6px;
            padding: 7px;
            selection-background-color: #71856f;
        }
        #bottomRunButton {
            background: #fffdf8;
            color: #3d382f;
            border: 1px solid #aa9f8f;
            border-radius: 6px;
            padding: 7px 12px;
        }
        QTabBar::tab {
            background: rgba(232, 224, 211, 220);
            color: #766d5f;
            padding: 8px 15px;
            border: none;
        }
        QTabBar::tab:selected {
            color: #3d382f;
            background: rgba(255, 252, 246, 220);
            border-top: 2px solid #71856f;
        }
        QScrollBar:vertical { background: transparent; width: 10px; }
        QScrollBar::handle:vertical {
            background: rgba(120, 110, 95, 125);
            border-radius: 4px;
            min-height: 24px;
        }
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
            height: 0;
        }
    )"));
}

} // namespace SF::GUI
