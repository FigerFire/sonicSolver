#pragma once

/// @file SF_controller.h
/// @brief Controller 层调度入口与协调器定义。

#include "SF_project.h"

#include <QObject>
#include <QProcess>

namespace SF::GUI {

class MainWindow;

/// @brief 协调项目文件、界面、求解器进程与可视化。
class GuiController : public QObject {
    Q_OBJECT

public:
    GuiController(MainWindow* view,
                  QString sourceRoot,
                  QObject* parent = nullptr);

    /// @brief 初始化空白项目界面。
    void initialize();

private:
    enum class SolverProcessPurpose {
        Solve,
        PreviewStep,
        InitialOutput,
        MeshGeneration,
        MeshCheck
    };

    void connectView();
    void connectProcesses();
    void requestWorkingDirectory();
    void openWorkspace(const QString& directory);
    void refreshProjectView();
    void generateInitialOutput();
    bool canGenerateInitialOutput(QString* reason = nullptr) const;
    void editConfigValue(int documentIndex,
                         int entryIndex,
                         const QString& value);
    void saveAll();
    void addOptionalModule(ProjectModuleKind kind);
    void generateMesh();
    void checkMesh();
    void requestImportSfm();
    void requestImportVtk();
    void requestImportStl();
    void saveGeometry(const GeometrySketch& sketch);

    void runSolver(bool previewOneStep);
    void createOutputTask(const QString& taskName);
    void submitOutputTask(const QString& taskName);
    void selectOutputTask(const QString& taskName);
    void showOutputResidual(const QString& taskName);
    void handleVisualizationLoaded(const QString& path);
    void refreshOutputTasks();
    QVector<OutputTaskInfo> collectOutputTasks() const;
    void stopSolver();
    void startSolverProcess(const QString& program,
                            const QStringList& arguments,
                            const QString& workingDirectory,
                            SolverProcessPurpose purpose);
    void appendSolverProcessOutput();
    void handleSolverFinished(int exitCode, QProcess::ExitStatus status);
    void runTerminalCommand(const QString& command);

    QString solverExecutable() const;
    bool parallelCommand(const QString& configPath,
                         QString* program,
                         QStringList* prefixArguments) const;
    QString preparePreviewCase(QString* error);
    QString latestPreviewOutput() const;
    static bool copyDirectory(const QString& source,
                              const QString& destination,
                              QString* error);
    void log(const QString& message);
    QString initialDirectory() const;

    MainWindow* view_ = nullptr;
    QString sourceRoot_;
    CaseProject project_;
    QString importedSfmPath_;
    QString previewDirectory_;
    QString generatedMeshName_;
    /// @brief 当前右侧 VTK 窗口实际显示的 result/<name>.pvd 对应任务名。
    QString displayedOutputTask_;
    SolverProcessPurpose solverPurpose_ = SolverProcessPurpose::Solve;
    QProcess solverProcess_;
    QProcess terminalProcess_;
};

} // namespace SF::GUI
