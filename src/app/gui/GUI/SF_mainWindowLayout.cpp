/// @file SF_mainWindowLayout.cpp
/// @brief GUI 主窗口布局、控件或视图切换实现。

#include "SF_mainWindow.h"

#include "SF_Geometry.h" // GeometryEditor 几何编辑器
#include "SF_toolbar.h" // ViewerToolbar 顶部工具栏
#include "SF_viewModeSwitch.h" // GEO/MESH 切换开关
#include "SF_VTKView.h"  // VTKView 可视化窗口      

#include <QFileDialog> // 文件/文件夹选择框
#include <QHeaderView> // 表头控件，比如树的
#include <QHBoxLayout> // 水平布局
#include <QIcon> // 图标
#include <QLabel> // 文本标签
#include <QLineEdit> // 单行输入框
#include <QMenu> // 菜单
#include <QMessageBox> // 弹窗
#include <QPlainTextEdit> // 多行纯文本框
#include <QPushButton> // 普通按钮
#include <QSignalBlocker> // 临时屏蔽信号
#include <QSplitter> // 可拖动分割器
#include <QStackedWidget> // 叠放页面控件
#include <QTabWidget> // 标签页控件
#include <QToolButton> // 工具按钮
#include <QTreeWidget> // 树形控件
#include <QVBoxLayout> // 垂直布局

#include <utility>

namespace SF::GUI {

// 右侧工作区文件：顶部模块栏、GEO/MESH 切换、中间视图、底部输出/终端/日志。

// 构建主窗口的 UI 结构，所有控件都在这里创建。
void MainWindow::buildUi()
{
    // 整体结构：
    // QSplitter(水平)
    //   左：sidebar，放“算例设置”和配置树
    //   右：rightWorkspace
    //       上：modeBar，放模块下拉和 GEO/MESH 开关
    //       中：viewerStack，几何编辑器或 VTK 视图二选一
    //       下：bottomTabs，求解器输出/终端/日志

    // 最外层分割器：左右分割，左侧 sidebar，右侧 rightWorkspace。
    // 左侧 sidebar 固定宽度，右侧 rightWorkspace 可拉伸。
    // auto*表示编译器自动推导指针类型，->表示指针调用，对象调用是.
    auto* workspaceSplitter = new QSplitter(Qt::Horizontal, this);
    workspaceSplitter->setChildrenCollapsible(false);// 设置分割器的子控件不可折叠到最小

    // 左侧 sidebar：放“算例设置”标题、工作目录按钮、当前工作目录文字、配置树。
    auto* sidebar = new QWidget(workspaceSplitter); // 创建一个 QWidget 作为左侧 sidebar 的容器
    sidebar->setObjectName(QStringLiteral("sidebar"));// 设置对象名称，便于样式表和调试
    // 左侧 sidebar 的布局：垂直布局，间距 9，边距 12。
    auto* sidebarLayout = new QVBoxLayout(sidebar); // 创建一个垂直布局管理器，管理 sidebar 内的控件
    sidebarLayout->setContentsMargins(12, 12, 12, 12); // 设置布局的边距为12
    sidebarLayout->setSpacing(9); // 设置布局的间距为9

    // 左侧标题行：放“算例设置”文字和右侧小文件夹按钮。
    auto* titleRow = new QWidget(sidebar); // 创建一个 QWidget 作为标题行的容器
    auto* titleLayout = new QHBoxLayout(titleRow); // 创建一个水平布局管理器，管理 titleRow 内的控件
    titleLayout->setContentsMargins(0, 0, 0, 0);// 设置布局的边距为0
    titleLayout->setSpacing(7);// 设置布局的间距为7

    auto* explorerTitle = new QLabel(QStringLiteral("算例设置"), titleRow); // 创建一个 QLabel 显示“算例设置”文字，父控件为 titleRow
    explorerTitle->setObjectName(QStringLiteral("sectionTitle")); // 设置对象名称，便于样式表和调试
    titleLayout->addWidget(explorerTitle); // 将 explorerTitle 添加到 titleLayout 中，使其显示在标题行的左侧
    titleLayout->addStretch(); // 添加一个弹性空间，使右侧的小文件夹按钮靠右显示

    // “算例设置”右侧的小文件夹按钮。
    // 点它时不直接打开文件夹，而是发 workingDirectoryRequested()，
    // Controller 收到后再调用 chooseWorkingDirectory()。
    auto* workspaceButton = new QToolButton(titleRow);
    workspaceButton->setObjectName(QStringLiteral("workspaceButton"));
    workspaceButton->setIcon(
        QIcon(QStringLiteral(":/sonic/icons/set-directory.png")));
    workspaceButton->setIconSize(QSize(20, 20));
    workspaceButton->setToolTip(QStringLiteral("设置工作文件夹"));
    workspaceButton->setAutoRaise(true);
    titleLayout->addWidget(workspaceButton);
    sidebarLayout->addWidget(titleRow);
    connect(workspaceButton, &QToolButton::clicked,
            this, &MainWindow::workingDirectoryRequested);

    // 当前工作目录显示文字；setWorkspace() 会更新它。
    workspaceLabel_ = new QLabel(QStringLiteral("请从“文件”设置工作文件夹"),
                                 sidebar);
    workspaceLabel_->setObjectName(QStringLiteral("mutedLabel"));
    workspaceLabel_->setWordWrap(true);
    workspaceLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    sidebarLayout->addWidget(workspaceLabel_);

    // 左侧配置树。第一列显示 key/section/module，第二列显示 value 或开关。
    // 树节点和真实配置项的关系由 TreeRole 存在 item->data() 里。
    configTree_ = new QTreeWidget(sidebar);
    configTree_->setColumnCount(2);
    configTree_->setHeaderHidden(true);
    configTree_->setRootIsDecorated(true);
    configTree_->setIndentation(15);
    configTree_->setIconSize(QSize(18, 18));
    configTree_->setContextMenuPolicy(Qt::CustomContextMenu);
    configTree_->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    configTree_->setTextElideMode(Qt::ElideRight);
    // Patch 名称允许用户自定义；第一列给出足够宽度，数值控件则集中在第二列。
    const QFontMetrics metrics(configTree_->font());
    configTree_->setColumnWidth(
        0, metrics.horizontalAdvance(QStringLiteral("customPatchName")) + 42);
    configTree_->header()->setSectionResizeMode(0, QHeaderView::Fixed);
    configTree_->header()->setStretchLastSection(true);
    sidebarLayout->addWidget(configTree_, 1);

    auto* rightWorkspace = new QWidget(workspaceSplitter);
    auto* rightLayout = new QVBoxLayout(rightWorkspace);
    rightLayout->setContentsMargins(0, 0, 0, 0);
    rightLayout->setSpacing(0);

    viewerToolbar_ = new ViewerToolbar(rightWorkspace);

    // 顶部模块下拉按钮，类似“当前正在看哪个设置模块”。
    // 菜单内容在 setProject() 里根据 modules_ 动态生成。
    moduleSelector_ = new QToolButton(viewerToolbar_);
    moduleSelector_->setObjectName(QStringLiteral("moduleSelector"));
    moduleSelector_->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    moduleSelector_->setPopupMode(QToolButton::InstantPopup);
    moduleSelector_->setIconSize(QSize(22, 22));
    moduleSelector_->setMinimumWidth(100);
    moduleSelectorMenu_ = new QMenu(moduleSelector_);
    moduleSelector_->setMenu(moduleSelectorMenu_);

    // GEO/MESH 开关：放在模块下拉旁边，作为当前模块视图的快速切换。
    // 它不直接保存任何配置，只切换中间 viewerStack_ 当前页面。
    viewModeSwitch_ = new ViewModeSwitch(viewerToolbar_);
    viewModeSwitch_->setToolTip(
        QStringLiteral("在几何编辑与网格/结果显示之间切换"));
    viewerToolbar_->setNavigationControls(moduleSelector_, viewModeSwitch_);
    rightLayout->addWidget(viewerToolbar_);

    auto* rightSplitter = new QSplitter(Qt::Vertical, rightWorkspace);
    rightSplitter->setObjectName(QStringLiteral("rightWorkspaceSplitter"));
    rightSplitter->setHandleWidth(5);
    rightSplitter->setChildrenCollapsible(false);

    // viewerStack_ 是“叠放页面”：同一块区域里只能显示一个页面。
    // showGeometryView() 选择 geometryEditor_，showVisualizationView() 选择 vtkView_。
    viewerStack_ = new QStackedWidget(rightSplitter);
    geometryEditor_ = new GeometryEditor(viewerStack_);
    vtkView_ = new VTKView(viewerToolbar_, viewerStack_);
    viewerStack_->addWidget(geometryEditor_);
    viewerStack_->addWidget(vtkView_);

    bottomTabs_ = new QTabWidget(rightSplitter);
    bottomTabs_->setObjectName(QStringLiteral("bottomOutputTabs"));
    bottomTabs_->setDocumentMode(true);
    bottomTabs_->setAttribute(Qt::WA_StyledBackground, true);
    // “求解器输出”页：显示 sonicSolver 进程 stdout/stderr 的核心日志。
    solverConsole_ = new QPlainTextEdit(bottomTabs_);
    solverConsole_->setObjectName(QStringLiteral("bottomConsole"));
    solverConsole_->setReadOnly(true);
    solverConsole_->setMaximumBlockCount(10000);
    bottomTabs_->addTab(solverConsole_, QStringLiteral("输出"));

    auto* terminalPage = new QWidget(bottomTabs_);
    terminalPage->setObjectName(QStringLiteral("bottomTerminalPage"));
    terminalPage->setAttribute(Qt::WA_StyledBackground, true);
    auto* terminalLayout = new QVBoxLayout(terminalPage);
    terminalLayout->setContentsMargins(0, 0, 0, 0);
    terminalLayout->setSpacing(0);
    terminalConsole_ = new QPlainTextEdit(terminalPage);
    terminalConsole_->setObjectName(QStringLiteral("bottomConsole"));
    terminalConsole_->setReadOnly(true);
    terminalConsole_->setMaximumBlockCount(5000);
    terminalLayout->addWidget(terminalConsole_, 1);
    auto* commandRow = new QHBoxLayout;
    commandRow->setContentsMargins(8, 6, 8, 8);
    commandRow->setSpacing(8);
    // “终端”页的输入框。按回车或点“运行”都会触发 runTerminal。
    terminalInput_ = new QLineEdit(terminalPage);
    terminalInput_->setObjectName(QStringLiteral("bottomCommandInput"));
    terminalInput_->setPlaceholderText(QStringLiteral("在工作目录运行命令"));
    auto* runCommandButton = new QPushButton(QStringLiteral("运行"), terminalPage);
    runCommandButton->setObjectName(QStringLiteral("bottomRunButton"));
    commandRow->addWidget(terminalInput_, 1);
    commandRow->addWidget(runCommandButton);
    terminalLayout->addLayout(commandRow);
    bottomTabs_->addTab(terminalPage, QStringLiteral("终端"));

    // “日志”页：显示 GUI 自己的状态提示，例如导入成功、刷新结果等。
    logConsole_ = new QPlainTextEdit(bottomTabs_);
    logConsole_->setObjectName(QStringLiteral("bottomConsole"));
    logConsole_->setReadOnly(true);
    logConsole_->setMaximumBlockCount(5000);
    bottomTabs_->addTab(logConsole_, QStringLiteral("日志"));
    // 底部标签横向显示，保持输出/终端/日志切换和普通桌面软件一致。
    bottomTabs_->setTabPosition(QTabWidget::North);

    // 加入分割器，分别上下切割 viewerStack_ 和 bottomTabs_。
    rightSplitter->addWidget(viewerStack_);
    rightSplitter->addWidget(bottomTabs_);
    // 设置上下的尺寸
    rightSplitter->setSizes({650, 230});
    rightSplitter->setStretchFactor(0, 1);
    rightLayout->addWidget(rightSplitter, 1);

    workspaceSplitter->addWidget(sidebar);
    workspaceSplitter->addWidget(rightWorkspace);
    workspaceSplitter->setSizes({330, 1110});
    workspaceSplitter->setStretchFactor(1, 1);
    setCentralWidget(workspaceSplitter);

    // 底部终端输入框的触发逻辑。
    // 注意：这里不直接运行 shell，只发 terminalCommandRequested(command)，
    // Controller 层负责启动 QProcess。
    auto runTerminal = [this] {
        const QString command = terminalInput_->text().trimmed();
        if (command.isEmpty()) return;
        terminalConsole_->appendPlainText(QStringLiteral("$ %1").arg(command));
        terminalInput_->clear();
        emit terminalCommandRequested(command);
    };
    connect(runCommandButton, &QPushButton::clicked, this, runTerminal);
    connect(terminalInput_, &QLineEdit::returnPressed, this, runTerminal);
}

void MainWindow::showGeometryView()
{
    // Geometry 与 Mesh 共用同一 VTK 上下文；Geometry 只切换为参数化 blockMesh 预览。
    // QSignalBlocker 防止 setChecked(false) 反过来再次触发 toggled 信号。
    viewerStack_->setCurrentWidget(vtkView_);
    vtkView_->showGeometryScene();
    viewerToolbar_->setWorkspaceMode(ViewerWorkspaceMode::Geometry);
    const QSignalBlocker blocker(viewModeSwitch_);
    viewModeSwitch_->setChecked(false);
    viewModeSwitch_->update();
}

void MainWindow::showVisualizationView()
{
    // 外部调用或菜单选择“网格与结果”时，切到 VTKView 并刷新 result。
    viewerStack_->setCurrentWidget(vtkView_);
    vtkView_->showResultScene();
    viewerToolbar_->setWorkspaceMode(ViewerWorkspaceMode::Mesh);
    const QSignalBlocker blocker(viewModeSwitch_);
    viewModeSwitch_->setChecked(true);
    viewModeSwitch_->update();
    vtkView_->activateMeshDefaults();
    vtkView_->refreshResult();
}

void MainWindow::setSolverRunning(bool running)
{
    // Controller 在启动/结束求解器时调用这里，控制“停止”动作能否点击。
    solverRunning_ = running;
    stopAction_->setEnabled(running);
    if (submitOutputTaskAction_) submitOutputTaskAction_->setEnabled(!running);
    rebuildTree();
}

void MainWindow::clearSolverOutput()
{
    // 每次新运行前清空旧输出，并自动切到底部第一个 tab。
    solverConsole_->clear();
    bottomTabs_->setCurrentWidget(solverConsole_);
}

void MainWindow::appendSolverOutput(const QString& text)
{
    // sonicSolver stdout/stderr 的文本最终会被 Controller 转发到这里。
    if (!text.isEmpty()) solverConsole_->appendPlainText(text.trimmed());
}

void MainWindow::appendTerminalOutput(const QString& text)
{
    // 底部终端命令的输出显示在这里。
    if (!text.isEmpty()) terminalConsole_->appendPlainText(text.trimmed());
}

void MainWindow::appendLog(const QString& text)
{
    // GUI 自身日志，例如“已导入 SFM”“PVD 解析失败”等。
    logConsole_->appendPlainText(text);
}

bool MainWindow::loadVisualization(const QString& path, QString* error)
{
    // 文件菜单“打开 PVD/VTK 结果”最终会到这里。
    focusModule(ProjectModuleKind::Mesh);
    showVisualizationView();
    return vtkView_->loadFile(path, error);
}

void MainWindow::setVisualizationLoadedHandler(
    std::function<void(const QString&)> handler)
{
    vtkView_->setVisualizationLoadedHandler(std::move(handler));
}

bool MainWindow::loadGeometry(const QString& path, QString* error)
{
    // blockMeshDict 是 Geometry 的唯一数据源；VTK 表示只是即时派生预览。
    return vtkView_->loadBlockMeshGeometry(path, error);
}

void MainWindow::watchResultDirectory(const QString& directory)
{
    // 打开项目后监听 result 目录，计算产生新 PVD 时可自动刷新。
    vtkView_->watchResultDirectory(directory);
}

void MainWindow::refreshVisualization()
{
    // 例如初始输出完成后，Controller 调用这里刷新 VTKView。
    vtkView_->refreshResult();
}

bool MainWindow::vtkAvailable() const
{
    // GUI 可以在没有 VTK 的构建中运行，只是中间预览不可用。
    return vtkView_->vtkAvailable();
}

QString MainWindow::chooseWorkingDirectory(const QString& initialDirectory)
{
    // 只负责弹 Qt 文件夹选择框，不保存选择结果。
    // 返回路径给 Controller，由 Controller 决定怎么打开项目。
    return QFileDialog::getExistingDirectory(
        this, QStringLiteral("设置 sonicSolver 工作文件夹"),
        initialDirectory,
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
}

QString MainWindow::chooseSfmFile(const QString& initialDirectory)
{
    // 弹出导入 SFM 的文件选择框。
    return QFileDialog::getOpenFileName(
        this, QStringLiteral("导入 SFM 网格"), initialDirectory,
        QStringLiteral("Sonic Fluid Mesh (*.sfm *.SFM)"));
}

QString MainWindow::chooseVtkFile(const QString& initialDirectory)
{
    // 弹出 VTK/PVD 结果文件选择框。
    return QFileDialog::getOpenFileName(
        this, QStringLiteral("打开 VTK 网格或结果"), initialDirectory,
        QStringLiteral("VTK 结果 (*.pvd *.vts *.vtu *.vtm)"));
}

QString MainWindow::chooseStlFile(const QString& initialDirectory)
{
    // 弹出 IBM STL 文件选择框。
    return QFileDialog::getOpenFileName(
        this, QStringLiteral("导入 IBM STL"), initialDirectory,
        QStringLiteral("STL geometry (*.stl *.STL)"));
}

void MainWindow::showWarning(const QString& title, const QString& message)
{
    // Controller 出错时用这个统一弹警告。
    QMessageBox::warning(this, title, message);
}

void MainWindow::showInformation(const QString& title,
                                 const QString& message)
{
    // Controller 需要普通提示时用这个。
    QMessageBox::information(this, title, message);
}

} // namespace SF::GUI
