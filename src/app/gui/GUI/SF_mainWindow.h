/// @file SF_mainWindow.h
/// @brief GUI 主窗口布局、控件或视图切换实现。

#pragma once

#include "SF_project.h"

#include <QMainWindow>

#include <functional>
#include <QStringList>
#include <QVector>

class QAction;
class QIcon;
class QLabel;
class QLineEdit;
class QMenu;
class QPlainTextEdit;
class QPushButton;
class QStackedWidget;
class QTabWidget;
class QToolButton;
class QTreeWidget;
class QTreeWidgetItem;

namespace SF::GUI {

class GeometryEditor;
class SwitchButton;
class ViewModeSwitch;
class ViewerToolbar;
class VTKView;

/// @brief sonicGui 的主窗口，也就是你看到的整个 GUI 外壳。
///
/// 这个类只负责“显示”和“发信号”，不直接改 case 文件，也不直接跑求解器。
/// 典型流程是：
/// 1. 用户点击按钮或菜单。
/// 2. MainWindow 发出一个 signal，例如 workingDirectoryRequested()。
/// 3. Controller 层接收 signal，真正执行打开文件夹、保存配置、启动求解器等操作。
/// 4. Controller 再调用 setProject()/appendLog() 等函数把结果显示回窗口。
class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    /// @brief 构造主窗口，内部会调用 buildUi()/buildMenus()/connectUi()。
    explicit MainWindow(QWidget* parent = nullptr);

    /// @brief 更新左侧工作目录文字，例如显示 sodCase。
    void setWorkspace(const QString& displayName, const QString& absolutePath);
    /// @brief 把 Controller 解析好的配置文档和模块列表显示到左侧树。
    void setProject(const QVector<ConfigDocument>& documents,
                    const QVector<ProjectModule>& modules);
    /// @brief 更新输出任务列表，任务来自 result/*.pvd 和当前 #job。
    void setOutputTasks(const QVector<OutputTaskInfo>& tasks);
    /// @brief 设置求解器是否正在运行；运行时“停止”菜单会变可点。
    void setSolverRunning(bool running);
    /// @brief 清空底部“求解器输出”页。
    void clearSolverOutput();
    /// @brief 在底部“求解器输出”页追加文本。
    void appendSolverOutput(const QString& text);
    /// @brief 在底部“终端”页追加文本。
    void appendTerminalOutput(const QString& text);
    /// @brief 在底部“日志”页追加文本。
    void appendLog(const QString& text);
    /// @brief 让中间 VTK 视图打开一个 .pvd/.vts/.vtu/.vtm 文件。
    bool loadVisualization(const QString& path, QString* error = nullptr);
    /// @brief 让中间几何编辑器打开 geometry/createMesh 配置。
    bool loadGeometry(const QString& path, QString* error = nullptr);
    /// @brief 监听 result 目录；有新的 PVD 时由 VTKView 自动刷新。
    void watchResultDirectory(const QString& directory);
    /// @brief 手动刷新 VTK 结果视图。
    void refreshVisualization();
    /// @brief 当前构建是否带 VTK 支持。
    bool vtkAvailable() const;
    /// @brief 切到几何编辑页面，并同步 GEO/MESH 开关。
    void showGeometryView();
    /// @brief 切到网格/结果显示页面，并同步 GEO/MESH 开关。
    void showVisualizationView();
    /// @brief 打开 case 后默认聚焦网格模块并同步顶部模块下拉。
    void focusMeshModule();
    /// @brief 注册右侧实际显示文件变化回调。
    void setVisualizationLoadedHandler(
        std::function<void(const QString&)> handler);
    /// @brief 聚焦左侧输出模块，不改变当前 3D/几何视图。
    void showOutputModule();

    /// @brief 弹出“选择工作文件夹”对话框。
    QString chooseWorkingDirectory(const QString& initialDirectory);
    /// @brief 弹出“导入 SFM 网格”对话框。
    QString chooseSfmFile(const QString& initialDirectory);
    /// @brief 弹出“打开 VTK/PVD 结果”对话框。
    QString chooseVtkFile(const QString& initialDirectory);
    /// @brief 弹出“导入 IBM STL”对话框。
    QString chooseStlFile(const QString& initialDirectory);
    /// @brief 弹出警告窗口。
    void showWarning(const QString& title, const QString& message);
    /// @brief 弹出普通提示窗口。
    void showInformation(const QString& title, const QString& message);

signals:
    /// @brief 用户想设置工作目录。来源：文件菜单、左侧标题旁的小文件夹按钮。
    void workingDirectoryRequested();
    /// @brief 用户点击“保存全部”。
    void saveAllRequested();
    /// @brief 用户点击“导入 SFM 网格”。
    void importSfmRequested();
    /// @brief 用户请求按当前 Mesh backend 生成网格。
    void generateMeshRequested();
    /// @brief 用户请求校验当前 OpenFOAM polyMesh。
    void checkMeshRequested();
    /// @brief 用户点击“打开 PVD/VTK 结果”。
    void importVtkRequested();
    /// @brief 用户点击“导入 IBM STL”。
    void importStlRequested();
    /// @brief 用户右键或菜单新增 IBM/MRF/Parallel 模块。
    void addModuleRequested(ProjectModuleKind kind);
    /// @brief 从初始/边界条件根节点新增 0/ 字段文件。
    void fieldAddRequested(ProjectModuleKind kind, const QString& fieldName);
    /// @brief 几何编辑器要求保存草图。
    void geometrySaveRequested(const GeometrySketch& sketch);
    void geometryHoverBlock(int blockId);
    void geometryHoverPatch(int patchId);
    void geometryHoverFace(int patchId, int faceId);
    void geometryHoverClear();
    /// @brief 用户点击运行；previewOneStep=true 表示只计算一步预览。
    void runRequested(bool previewOneStep);
    /// @brief 用户在输出模块中新建任务。
    void outputTaskCreateRequested(const QString& taskName);
    /// @brief 用户请求提交某个输出任务计算。
    void outputTaskSubmitRequested(const QString& taskName);
    /// @brief 用户选中某个任务，要求加载对应 PVD。
    void outputTaskSelected(const QString& taskName);
    /// @brief 用户请求查看某个任务的残差。
    void outputResidualRequested(const QString& taskName);
    /// @brief 用户点击停止。
    void stopRequested();
    /// @brief 用户改了一个已有配置项的值。
    /// @param documentIndex 对应 documents_ 里的第几个配置文档。
    /// @param entryIndex 对应该文档里的第几个配置项。
    /// @param value 新值，Controller 会负责写回文件。
    void configValueEdited(int documentIndex,
                           int entryIndex,
                           const QString& value);
    /// @brief 用户通过右键菜单新增配置行。
    /// @param newSection true 表示新建一个重复 section，例如新的 [blocks]。
    void configAssignmentsRequested(int documentIndex,
                                    const QString& section,
                                    const QStringList& assignments,
                                    bool newSection);
    /// @brief 用户请求新增或编辑 blockMesh patch 及其四边形面。
    void blockMeshPatchAddRequested(int documentIndex, const QString& name,
                                    const QString& type);
    void blockMeshPatchFaceAddRequested(int documentIndex,
                                        const QString& patchName,
                                        const QString& face);
    void blockMeshPatchRenameRequested(int documentIndex,
                                       const QString& oldName,
                                       const QString& newName);
    void blockMeshPatchRemoveRequested(int documentIndex, const QString& name);
    void configEntryRemoveRequested(int documentIndex, int entryIndex);
    /// @brief 用户在底部终端输入了一条命令。
    void terminalCommandRequested(const QString& command);

private:
    /// @brief QTreeWidgetItem 上自定义保存的数据角色。
    ///
    /// Qt 的树节点只有显示文字还不够，所以用 setData()/data() 把“这个节点
    /// 对应哪个配置文档、哪个配置项、哪个模块”也藏在节点里。
    enum TreeRole {
        /// @brief 节点属于 documents_ 中的第几个 ConfigDocument。
        DocumentRole = Qt::UserRole + 1,
        /// @brief 节点属于 ConfigDocument::entries() 中的第几个 ConfigEntry。
        EntryRole,
        /// @brief 节点属于 modules_ 中的第几个 ProjectModule。
        ModuleRole,
        /// @brief 节点代表哪个 section，例如 vertices、solver、velocity。
        SectionRole,
        /// @brief 输出任务名。
        OutputTaskRole,
        /// @brief 输出任务动作：new/submit/residual。
        OutputActionRole,
        /// @brief 属性树中用于比较的原始显示值，避免重复提交未变更内容。
        OriginalValueRole,
        /// @brief blockMesh patch 的名称，用于右键新增/删除 face 与 patch。
        BlockIDRole,
        PatchIDRole,
        FaceIDRole,
        PatchNameRole
    };

    /// @brief 创建所有可见控件：左树、顶部栏、中间视图、底部输出页。
    void buildUi();
    /// @brief 创建 macOS 顶部菜单：文件、模块、运行、视图、已启用。
    void buildMenus();
    /// @brief 设置全局 stylesheet，也就是颜色、圆角、边框、字体风格。
    void applyStyle();
    /// @brief 把控件信号和 MainWindow 槽/信号连起来。
    void connectUi();
    /// @brief 根据 documents_ 和 modules_ 重建左侧配置树。
    void rebuildTree();
    /// @brief 双击树节点时编辑该配置项。
    void editTreeItem(QTreeWidgetItem* item);
    /// @brief 提交树内第二列的直接编辑。
    void commitTreeItemEdit(QTreeWidgetItem* item, int column);
    /// @brief 点击输出模块里的任务或动作节点。
    bool handleOutputTreeClick(QTreeWidgetItem* item);
    void handleTreeHover(QTreeWidgetItem* item);
    /// @brief 弹出输入框并创建新的输出任务。
    void requestCreateOutputTask();
    /// @brief 当前激活的输出任务名。
    QString currentOutputTaskName() const;
    /// @brief 顶部模块下拉选择后，展开并滚动到对应模块。
    void focusModule(int moduleIndex);
    /// @brief 按业务模块类型展开左侧树；热键和菜单都复用这个入口。
    void focusModule(ProjectModuleKind kind);
    /// @brief 左侧树选中模块后，同步顶部模块按钮文字和图标。
    void syncModuleSelector(int moduleIndex);
    /// @brief 更新菜单栏“已启用”，列出所有 true 的开关项。
    void updateEnabledMenu();
    /// @brief 从“已启用”菜单跳回对应配置树节点。
    void selectEntry(int documentIndex, int entryIndex);
    /// @brief 左侧树右键菜单：新增顶点、block、边界条件、IBM/MRF 等。
    void showModuleContextMenu(const QPoint& position);
    /// @brief 给某个菜单加入“新增 IBM/MRF/Parallel”三项。
    void addOptionalModuleActions(QMenu* menu);
    /// @brief 某些 type 字段只能从固定枚举里选，这里返回候选列表。
    QStringList enumChoices(const ConfigEntry& entry) const;
    /// @brief 把配置文档里的内部字段名翻译成界面显示名称。
    QString entryLabel(const ConfigEntry& entry) const;
    /// @brief 根据模块类型返回左侧树和顶部菜单使用的图标。
    static QIcon moduleIcon(ProjectModuleKind kind);

    /// @brief 当前打开项目的配置文档。
    QVector<ConfigDocument> documents_;
    /// @brief 左侧树顶层模块：几何、网格、求解器、湍流模型等。
    QVector<ProjectModule> modules_;
    /// @brief 输出模块的任务列表，由 Controller 从 result/*.pvd 收集。
    QVector<OutputTaskInfo> outputTasks_;
    /// @brief 求解器是否正在运行，用于禁用任务提交按钮。
    bool solverRunning_ = false;
    /// @brief 重建树时屏蔽 itemChanged，避免把显示刷新误当成用户输入。
    bool rebuildingTree_ = false;

    /// @brief 左侧配置树；用户主要在这里查看和编辑 case 配置。
    QTreeWidget* configTree_ = nullptr;
    /// @brief 左侧显示当前工作目录名称的文字。
    QLabel* workspaceLabel_ = nullptr;
    /// @brief 顶部“设置”旁的模块下拉按钮。
    QToolButton* moduleSelector_ = nullptr;
    /// @brief moduleSelector_ 弹出的菜单。
    QMenu* moduleSelectorMenu_ = nullptr;
    /// @brief 顶部 GEO/MESH 双态开关。
    ViewModeSwitch* viewModeSwitch_ = nullptr;
    /// @brief 顶部三行可视化工具栏。
    ViewerToolbar* viewerToolbar_ = nullptr;
    /// @brief 当前展开/选中的模块序号。
    int currentModuleIndex_ = 0;
    /// @brief 中间几何编辑器页面。
    GeometryEditor* geometryEditor_ = nullptr;
    /// @brief 中间 VTK 网格/结果显示页面。
    VTKView* vtkView_ = nullptr;
    /// @brief 中间页面栈；在 geometryEditor_ 和 vtkView_ 之间切换。
    QStackedWidget* viewerStack_ = nullptr;
    /// @brief 底部 tab：求解器输出、终端、日志。
    QTabWidget* bottomTabs_ = nullptr;
    /// @brief 底部“求解器输出”文本框。
    QPlainTextEdit* solverConsole_ = nullptr;
    /// @brief 底部“日志”文本框。
    QPlainTextEdit* logConsole_ = nullptr;
    /// @brief 底部“终端”输出框。
    QPlainTextEdit* terminalConsole_ = nullptr;
    /// @brief 底部“终端”命令输入框。
    QLineEdit* terminalInput_ = nullptr;
    /// @brief 菜单栏“已启用”菜单。
    QMenu* enabledMenu_ = nullptr;
    /// @brief 运行菜单里的“停止”动作；求解器运行时才可用。
    QAction* stopAction_ = nullptr;
    QAction* submitOutputTaskAction_ = nullptr;
};

} // namespace SF::GUI
