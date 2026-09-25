/// @file SF_controller.cpp
/// @brief GUI 控制器对模型、视图与用户操作的协调。

#include "SF_controller.h"

#include "SF_mainWindow.h"

#include <QCoreApplication>
#include <QDateTime>
#include <QDir>
#include <QDirIterator>
#include <QFile>
#include <QFileInfo>
#include <QRegularExpression>
#include <QSet>
#include <QStandardPaths>
#include <QTimer>
#include <QXmlStreamReader>

#include <algorithm>
#include <utility>

namespace SF::GUI {
namespace {

QString readTextFile(const QString& path, QString* error = nullptr)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return {};
    }
    QString text = QString::fromUtf8(file.readAll());
    text.replace(QStringLiteral("\r\n"), QStringLiteral("\n"));
    return text;
}

bool writeTextFile(const QString& path,
                   const QStringList& lines,
                   QString* error = nullptr)
{
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text
                   | QIODevice::Truncate)) {
        if (error) *error = file.errorString();
        return false;
    }
    QByteArray data = lines.join(QLatin1Char('\n')).toUtf8();
    data.append('\n');
    if (file.write(data) != data.size()) {
        if (error) *error = file.errorString();
        return false;
    }
    return true;
}

struct PvdSummary {
    int frameCount = 0;
    QString latestTimeText;
    double latestTime = 0.0;
    bool hasLatestTime = false;
};

QString caseDirectoryFromPath(const QString& casePath)
{
    const QFileInfo info(casePath);
    if (info.isDir()) return info.absoluteFilePath();
    if (info.fileName().compare(QStringLiteral("controlDict"),
                                Qt::CaseInsensitive) == 0) {
        return QDir(info.dir().filePath(QStringLiteral(".."))).absolutePath();
    }
    return info.absolutePath();
}

QString stripFoamComments(QString text)
{
    text.remove(QRegularExpression(QStringLiteral(R"(/\*[\s\S]*?\*/)")));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\n]*)")));
    return text;
}

QString foamValue(const QString& content, const QString& key)
{
    const QString clean = stripFoamComments(content);
    const QRegularExpression pattern(
        QStringLiteral(R"((?:^|[;\n])\s*%1\s+([^;{}]+);)")
            .arg(QRegularExpression::escape(key)),
        QRegularExpression::CaseInsensitiveOption);
    const QRegularExpressionMatch match = pattern.match(clean);
    return match.hasMatch() ? match.captured(1).trimmed() : QString();
}

QString foamBlock(const QString& content, const QString& key)
{
    const QString clean = stripFoamComments(content);
    const QRegularExpression startPattern(
        QStringLiteral(R"(\b%1\b\s*\{)")
            .arg(QRegularExpression::escape(key)),
        QRegularExpression::CaseInsensitiveOption);
    const QRegularExpressionMatch start = startPattern.match(clean);
    if (!start.hasMatch()) return {};

    int depth = 0;
    const int begin = start.capturedEnd(0) - 1;
    for (int i = begin; i < clean.size(); ++i) {
        if (clean.at(i) == QLatin1Char('{')) {
            if (depth == 0) {
                ++depth;
                continue;
            }
            ++depth;
        } else if (clean.at(i) == QLatin1Char('}')) {
            --depth;
            if (depth == 0) return clean.mid(begin + 1, i - begin - 1);
        }
    }
    return {};
}

QString foamBlockValue(const QString& content,
                       const QString& block,
                       const QString& key)
{
    return foamValue(foamBlock(content, block), key);
}

QStringList foamWordList(QString text)
{
    text = text.trimmed();
    if (text.startsWith(QLatin1Char('(')) && text.endsWith(QLatin1Char(')'))) {
        text = text.mid(1, text.size() - 2);
    }
    text.replace(QLatin1Char(';'), QLatin1Char(' '));
    text.replace(QLatin1Char(','), QLatin1Char(' '));
    return text.split(QRegularExpression(QStringLiteral("\\s+")),
                      Qt::SkipEmptyParts);
}

QString foamWord(QString text)
{
    text = text.trimmed();
    text.remove(QLatin1Char(';'));
    if ((text.startsWith(QLatin1Char('"')) && text.endsWith(QLatin1Char('"')))
        || (text.startsWith(QLatin1Char('\'')) && text.endsWith(QLatin1Char('\'')))) {
        text = text.mid(1, text.size() - 2);
    }
    const QStringList words = text.split(
        QRegularExpression(QStringLiteral("\\s+")), Qt::SkipEmptyParts);
    return words.isEmpty() ? QString() : words.first();
}

enum class CaseBackendKind {
    SonicSolver,
    OpenFOAMScript,
    Unsupported
};

struct CaseBackend {
    CaseBackendKind kind = CaseBackendKind::Unsupported;
    QString application;
    QString scriptPath;
    QString error;
};

bool isSonicApplication(const QString& application)
{
    const QString name = application.toLower();
    return name == QStringLiteral("sonicsolver")
        || name == QStringLiteral("densitybasesolver")
        || name == QStringLiteral("pressurebasesolver")
        || name == QStringLiteral("createmesh")
        || name == QStringLiteral("meshgenerator");
}

bool hasSonicDictionaries(const QString& caseDir)
{
    const QDir dir(caseDir);
    return QFileInfo(dir.filePath(QStringLiteral("case.yaml"))).isFile()
        || (QFileInfo(dir.filePath(QStringLiteral("system/outputDict"))).isFile()
            && QFileInfo(dir.filePath(QStringLiteral("system/sonicDict"))).isFile());
}

CaseBackend detectCaseBackend(const QString& casePath)
{
    CaseBackend backend;
    const QString caseDir = caseDirectoryFromPath(casePath);
    if (hasSonicDictionaries(caseDir)) {
        backend.kind = CaseBackendKind::SonicSolver;
        backend.application = QStringLiteral("sonicSolver");
        return backend;
    }
    const QString controlPath =
        QDir(caseDir).filePath(QStringLiteral("system/controlDict"));
    const QString content = readTextFile(controlPath);
    if (content.isEmpty()) {
        backend.error = QStringLiteral("未找到或无法读取 system/controlDict");
        return backend;
    }

    backend.application = foamWord(foamValue(content, QStringLiteral("application")));
    if (backend.application.isEmpty()) {
        backend.application = foamWord(foamValue(content, QStringLiteral("solver")));
    }
    if (backend.application.isEmpty()) {
        backend.error =
            QStringLiteral("system/controlDict 未设置 application 或 solver");
        return backend;
    }

    if (isSonicApplication(backend.application)) {
        backend.kind = CaseBackendKind::SonicSolver;
        return backend;
    }

    const QDir dir(caseDir);
    const QString allrun = dir.filePath(QStringLiteral("Allrun"));
    const QString run = dir.filePath(QStringLiteral("run"));
    if (QFileInfo(allrun).isFile()) {
        backend.kind = CaseBackendKind::OpenFOAMScript;
        backend.scriptPath = allrun;
        return backend;
    }
    if (QFileInfo(run).isFile()) {
        backend.kind = CaseBackendKind::OpenFOAMScript;
        backend.scriptPath = run;
        return backend;
    }

    backend.error =
        QStringLiteral("application=%1 不是 sonicSolver 后端，且 case 目录中没有 Allrun/run 脚本")
            .arg(backend.application);
    return backend;
}

QString openFoamPointerFile(const QString& caseDir, QString* error)
{
    const QDir dir(caseDir);
    if (!QFileInfo(dir.filePath(QStringLiteral("constant/polyMesh"))).isDir()) {
        return {};
    }

    const QFileInfoList existing = dir.entryInfoList(
        {QStringLiteral("*.foam"), QStringLiteral("*.FOAM")},
        QDir::Files, QDir::Name);
    for (const QFileInfo& info : existing) {
        if (!info.fileName().startsWith(QStringLiteral("._"))) {
            return info.absoluteFilePath();
        }
    }

    const QString pointer =
        dir.filePath(QStringLiteral(".sonicGui.foam"));
    QFile file(pointer);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text
                   | QIODevice::Truncate)) {
        if (error) *error = file.errorString();
        return {};
    }
    file.write("// sonicGui OpenFOAM reader pointer\n");
    return pointer;
}

QString meshPreviewCandidate(const QString& caseDir, QString* error = nullptr)
{
    const QDir dir(caseDir);
    if (hasSonicDictionaries(caseDir)) {
        const bool native = QFileInfo(
            dir.filePath(QStringLiteral("case.yaml"))).isFile();
        const QDir meshDir(dir.filePath(native ? QStringLiteral("mesh")
                                               : QStringLiteral("constant")));
        const QStringList preferred{
            QStringLiteral("mesh.sfm"), QStringLiteral("mesh.SFM")};
        for (const QString& name : preferred) {
            const QFileInfo mesh(meshDir.filePath(name));
            if (mesh.isFile()) return mesh.absoluteFilePath();
        }
        const QFileInfoList sfmFiles = meshDir.entryInfoList(
            {QStringLiteral("*.sfm"), QStringLiteral("*.SFM")},
            QDir::Files, QDir::Name);
        for (const QFileInfo& mesh : sfmFiles) {
            if (!mesh.fileName().startsWith(QStringLiteral("._"))) {
                return mesh.absoluteFilePath();
            }
        }
        if (error) {
            *error = QStringLiteral(
                "检测到 sonicSolver case，但 mesh 目录下没有 SFM 网格文件。");
        }
        return {};
    }

    // 没有 sonic 字典对时按 OpenFOAM case 处理，只从真实 constant/polyMesh 读取。
    return openFoamPointerFile(caseDir, error);
}

PvdSummary summarizePvd(const QString& path)
{
    PvdSummary summary;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return summary;

    QXmlStreamReader xml(&file);
    while (!xml.atEnd()) {
        xml.readNext();
        if (!xml.isStartElement()
            || xml.name().compare(QStringLiteral("DataSet"),
                                  Qt::CaseInsensitive) != 0) {
            continue;
        }
        ++summary.frameCount;
        const QString timeText =
            xml.attributes().value(QStringLiteral("timestep")).toString();
        bool ok = false;
        const double time = timeText.toDouble(&ok);
        if (ok) {
            summary.latestTimeText = timeText;
            summary.latestTime = time;
            summary.hasLatestTime = true;
        }
    }
    return summary;
}

} // namespace

GuiController::GuiController(MainWindow* view,
                             QString sourceRoot,
                             QObject* parent)
    : QObject(parent)
    , view_(view)
    , sourceRoot_(std::move(sourceRoot))
{
    connectView();
    connectProcesses();
}

void GuiController::initialize()
{
    const QString testWorkspace =
        qEnvironmentVariable("SONIC_GUI_WORKSPACE");
    if (!testWorkspace.isEmpty()) {
        openWorkspace(testWorkspace);
        const QStringList testModules =
            qEnvironmentVariable("SONIC_GUI_MODULES")
                .split(QLatin1Char(','), Qt::SkipEmptyParts);
        for (const QString& module : testModules) {
            if (module.compare(QStringLiteral("IBM"),
                               Qt::CaseInsensitive) == 0) {
                addOptionalModule(ProjectModuleKind::IBM);
            } else if (module.compare(QStringLiteral("MRF"),
                                      Qt::CaseInsensitive) == 0) {
                addOptionalModule(ProjectModuleKind::MRF);
            } else if (module.compare(QStringLiteral("Parallel"),
                                      Qt::CaseInsensitive) == 0) {
                addOptionalModule(ProjectModuleKind::Parallel);
            }
        }
        return;
    }
    log(QStringLiteral("请先从“文件”菜单设置工作文件夹。"));
    if (!view_->vtkAvailable()) {
        log(QStringLiteral("当前未检测到 VTK；安装后可显示 VTS/VTM/STL。"));
    }
}

void GuiController::connectView()
{
    connect(view_, &MainWindow::workingDirectoryRequested,
            this, &GuiController::requestWorkingDirectory);
    connect(view_, &MainWindow::saveAllRequested,
            this, &GuiController::saveAll);
    connect(view_, &MainWindow::importSfmRequested,
            this, &GuiController::requestImportSfm);
    connect(view_, &MainWindow::generateMeshRequested,
            this, &GuiController::generateMesh);
    connect(view_, &MainWindow::checkMeshRequested,
            this, &GuiController::checkMesh);
    connect(view_, &MainWindow::importVtkRequested,
            this, &GuiController::requestImportVtk);
    connect(view_, &MainWindow::importStlRequested,
            this, &GuiController::requestImportStl);
    connect(view_, &MainWindow::addModuleRequested,
            this, &GuiController::addOptionalModule);
    connect(view_, &MainWindow::geometrySaveRequested,
            this, &GuiController::saveGeometry);
    connect(view_, &MainWindow::runRequested,
            this, &GuiController::runSolver);
    connect(view_, &MainWindow::outputTaskCreateRequested,
            this, &GuiController::createOutputTask);
    connect(view_, &MainWindow::outputTaskSubmitRequested,
            this, &GuiController::submitOutputTask);
    connect(view_, &MainWindow::outputTaskSelected,
            this, &GuiController::selectOutputTask);
    connect(view_, &MainWindow::outputResidualRequested,
            this, &GuiController::showOutputResidual);
    connect(view_, &MainWindow::stopRequested,
            this, &GuiController::stopSolver);
    connect(view_, &MainWindow::configValueEdited,
            this, &GuiController::editConfigValue);
    connect(view_, &MainWindow::configAssignmentsRequested,
            this, &GuiController::appendConfigAssignments);
    connect(view_, &MainWindow::blockMeshPatchAddRequested,
            this, &GuiController::addBlockMeshPatch);
    connect(view_, &MainWindow::blockMeshPatchFaceAddRequested,
            this, &GuiController::addBlockMeshPatchFace);
    connect(view_, &MainWindow::blockMeshPatchRenameRequested,
            this, &GuiController::renameBlockMeshPatch);
    connect(view_, &MainWindow::blockMeshPatchRemoveRequested,
            this, &GuiController::removeBlockMeshPatch);
    connect(view_, &MainWindow::configEntryRemoveRequested,
            this, &GuiController::removeConfigEntry);
    connect(view_, &MainWindow::fieldAddRequested,
            this, &GuiController::addField);
    view_->setVisualizationLoadedHandler(
        [this](const QString& path) { handleVisualizationLoaded(path); });
    connect(view_, &MainWindow::terminalCommandRequested,
            this, &GuiController::runTerminalCommand);
}

void GuiController::connectProcesses()
{
    connect(&solverProcess_, &QProcess::readyReadStandardOutput,
            this, &GuiController::appendSolverProcessOutput);
    connect(&solverProcess_, &QProcess::readyReadStandardError,
            this, &GuiController::appendSolverProcessOutput);
    connect(&solverProcess_,
            qOverload<int, QProcess::ExitStatus>(&QProcess::finished),
            this, &GuiController::handleSolverFinished);
    connect(&solverProcess_, &QProcess::errorOccurred, this,
            [this](QProcess::ProcessError) {
                log(QStringLiteral("求解器进程错误：%1")
                        .arg(solverProcess_.errorString()));
            });

    connect(&terminalProcess_, &QProcess::readyReadStandardOutput, this, [this] {
        view_->appendTerminalOutput(
            QString::fromLocal8Bit(terminalProcess_.readAllStandardOutput()));
    });
    connect(&terminalProcess_, &QProcess::readyReadStandardError, this, [this] {
        view_->appendTerminalOutput(
            QString::fromLocal8Bit(terminalProcess_.readAllStandardError()));
    });
}

QString GuiController::initialDirectory() const
{
    return project_.directory().isEmpty()
        ? QDir(sourceRoot_).filePath(QStringLiteral("test"))
        : project_.directory();
}

void GuiController::requestWorkingDirectory()
{
    const QString directory = view_->chooseWorkingDirectory(initialDirectory());
    if (!directory.isEmpty()) openWorkspace(directory);
}

void GuiController::openWorkspace(const QString& directory)
{
    if (solverProcess_.state() != QProcess::NotRunning) {
        view_->showInformation(QStringLiteral("求解器正在运行"),
                               QStringLiteral("请先停止当前进程。"));
        return;
    }

    const QDir selectedDirectory(directory);
    const bool existingCase = QFileInfo(selectedDirectory.filePath(
            QStringLiteral("system/controlDict"))).isFile()
        || QFileInfo(selectedDirectory.filePath(
            QStringLiteral("system/outputDict"))).isFile()
        || QFileInfo(selectedDirectory.filePath(
            QStringLiteral("system/sonicDict"))).isFile()
        || QFileInfo(selectedDirectory.filePath(
            QStringLiteral("constant/polyMesh"))).isDir();
    QString error;
    if (!project_.open(directory, &error)) {
        view_->showWarning(QStringLiteral("无法打开工作文件夹"), error);
        return;
    }
    importedSfmPath_.clear();
    displayedOutputTask_.clear();
    refreshProjectView();
    log(QStringLiteral("工作文件夹：%1").arg(project_.directory()));
    if (existingCase) {
        // 先按 case 类型读取真实网格：sonic 是 constant/*.sfm，其他一律是 polyMesh。
        const QString preview = meshPreviewCandidate(project_.directory(), &error);
        if (!preview.isEmpty() && view_->loadVisualization(preview, &error)) {
            log(QStringLiteral("已读取网格预览：%1").arg(preview));
        } else if (!error.isEmpty()) {
            log(QStringLiteral("网格预览未加载：%1").arg(error));
        }
        if (project_.isSonicCase() && preview.isEmpty()) generateInitialOutput();
        const QString firstJobResult = QDir(project_.resultDirectory()).filePath(
            project_.outputTaskName() + QStringLiteral(".pvd"));
        if (QFileInfo(firstJobResult).isFile()) {
            error.clear();
            if (view_->loadVisualization(firstJobResult, &error)) {
                log(QStringLiteral("默认显示输出任务：%1")
                        .arg(project_.outputTaskName()));
            } else if (!error.isEmpty()) {
                log(error);
            }
        }
    } else {
        log(QStringLiteral("空工作文件夹已初始化；尚无输入网格，不生成第 0 步结果。"));
    }
    // case 打开后的导航状态始终落在网格；PVD 行用独立底色表达当前显示结果。
    view_->focusMeshModule();
}

void GuiController::refreshProjectView()
{
    const QFileInfo info(project_.directory());
    view_->setWorkspace(info.fileName(), info.absoluteFilePath());
    view_->setProject(project_.documents(), project_.modules());
    view_->watchResultDirectory(project_.resultDirectory());
    refreshOutputTasks();

    const QString geometryPath =
        project_.moduleFile(ProjectModuleKind::Geometry);
    if (!geometryPath.isEmpty()) {
        QString error;
        if (!view_->loadGeometry(geometryPath, &error)) {
            log(QStringLiteral("几何显示未更新：%1").arg(error));
        }
    }
}

void GuiController::generateInitialOutput()
{
    if (project_.casePath().isEmpty()
        || solverProcess_.state() != QProcess::NotRunning) {
        return;
    }

    const CaseBackend backend = detectCaseBackend(project_.casePath());
    if (backend.kind != CaseBackendKind::SonicSolver) {
        const QString reason = backend.error.isEmpty()
            ? QStringLiteral("当前 case 由 OpenFOAM 脚本后端运行")
            : backend.error;
        log(QStringLiteral("跳过第 0 步输出：%1").arg(reason));
        return;
    }

    QString unavailableReason;
    if (!canGenerateInitialOutput(&unavailableReason)) {
        log(QStringLiteral("跳过第 0 步输出：%1").arg(unavailableReason));
        return;
    }

    const QString executable = solverExecutable();
    if (!QFileInfo::exists(executable)) {
        log(QStringLiteral("无法生成第 0 步：找不到求解器 %1")
                .arg(executable));
        return;
    }

    QString program = executable;
    QStringList arguments;
    QStringList prefixArguments;
    QString parallelPath = parallelControlFile(project_.casePath());
    if (parallelPath.isEmpty()) {
        parallelPath = solverControlFile(project_.casePath());
    }
    if (parallelCommand(parallelPath, &program, &prefixArguments)) {
        arguments = prefixArguments;
        arguments << executable
                  << QStringLiteral("--initial-output")
                  << project_.casePath();
    } else {
        arguments << QStringLiteral("--initial-output")
                  << project_.casePath();
    }
    startSolverProcess(program, arguments, project_.directory(),
                       SolverProcessPurpose::InitialOutput);
}

bool GuiController::canGenerateInitialOutput(QString* reason) const
{
    const QString caseDir = caseDirectoryFromPath(project_.casePath());
    const QString sonicPath =
        QDir(caseDir).filePath(QStringLiteral("system/sonicDict"));
    const QString content = readTextFile(sonicPath);
    if (content.isEmpty()) {
        if (reason) *reason = QStringLiteral("system/sonicDict 为空");
        return false;
    }

    const QString createMeshEnabled =
        foamBlockValue(content, QStringLiteral("createMesh"),
                       QStringLiteral("enabled"));
    const bool caseCreateMeshEnabled =
        createMeshEnabled.compare(QStringLiteral("true"),
                                  Qt::CaseInsensitive) == 0;
    if (caseCreateMeshEnabled) {
        if (reason) *reason = QStringLiteral("当前是 createMesh 模式");
        return false;
    }

    QStringList meshFiles =
        foamWordList(foamValue(content, QStringLiteral("meshFiles")));
    if (meshFiles.isEmpty()) {
        meshFiles << QStringLiteral("constant/mesh.sfm");
    }
    if (meshFiles.isEmpty()) {
        if (reason) *reason = QStringLiteral("尚未配置网格文件");
        return false;
    }
    const QDir caseDirectory(caseDir);
    for (const QString& meshFile : meshFiles) {
        const QString path = QFileInfo(meshFile).isAbsolute()
            ? meshFile : caseDirectory.filePath(meshFile);
        if (!QFileInfo(path).isFile()) {
            if (reason) {
                *reason = QStringLiteral("网格文件不存在：%1").arg(path);
            }
            return false;
        }
    }

    const QString solverCreateMesh =
        foamBlockValue(content, QStringLiteral("solver"),
                       QStringLiteral("createMesh"));
    if (solverCreateMesh.compare(QStringLiteral("true"),
                                 Qt::CaseInsensitive) == 0) {
        if (reason) *reason = QStringLiteral("当前是 createMesh 模式");
        return false;
    }
    return true;
}

void GuiController::editConfigValue(int documentIndex,
                                    int entryIndex,
                                    const QString& value)
{
    QString error;
    if (!project_.setValue(documentIndex, entryIndex, value, &error)) {
        view_->showWarning(QStringLiteral("保存失败"), error);
        return;
    }
    refreshProjectView();
    log(QStringLiteral("配置已保存。"));
}

void GuiController::appendConfigAssignments(
    int documentIndex,
    const QString& section,
    const QStringList& assignments,
    bool newSection)
{
    QString error;
    if (!project_.appendAssignments(documentIndex, section, assignments,
                                    newSection, &error)) {
        view_->showWarning(QStringLiteral("新增配置失败"), error);
        return;
    }
    refreshProjectView();
    log(QStringLiteral("[%1] 已新增配置。").arg(section));
}

void GuiController::addBlockMeshPatch(int documentIndex, const QString& name,
                                      const QString& type)
{
    QString error;
    if (!project_.addBlockMeshPatch(documentIndex, name, type, &error)) {
        view_->showWarning(QStringLiteral("新增 patch 失败"), error);
        return;
    }
    refreshProjectView();
}

void GuiController::addBlockMeshPatchFace(int documentIndex,
                                          const QString& patchName,
                                          const QString& face)
{
    QString error;
    if (!project_.addBlockMeshPatchFace(documentIndex, patchName, face, &error)) {
        view_->showWarning(QStringLiteral("新增 patch 面失败"), error);
        return;
    }
    refreshProjectView();
}

void GuiController::renameBlockMeshPatch(int documentIndex,
                                         const QString& oldName,
                                         const QString& newName)
{
    if (oldName == newName) return;
    QString error;
    if (!project_.renameBlockMeshPatch(documentIndex, oldName, newName, &error)) {
        view_->showWarning(QStringLiteral("重命名 patch 失败"), error);
        return;
    }
    refreshProjectView();
}

void GuiController::removeBlockMeshPatch(int documentIndex, const QString& name)
{
    QString error;
    if (!project_.removeBlockMeshPatch(documentIndex, name, &error)) {
        view_->showWarning(QStringLiteral("删除 patch 失败"), error);
        return;
    }
    refreshProjectView();
}

void GuiController::removeConfigEntry(int documentIndex, int entryIndex)
{
    QString error;
    if (!project_.removeEntry(documentIndex, entryIndex, &error)) {
        view_->showWarning(QStringLiteral("删除配置失败"), error);
        return;
    }
    refreshProjectView();
}

void GuiController::saveAll()
{
    if (project_.directory().isEmpty()) return;
    QString error;
    if (!project_.saveAll(&error)) {
        view_->showWarning(QStringLiteral("保存失败"), error);
        return;
    }
    log(QStringLiteral("已保存全部项目配置。"));
}

void GuiController::addOptionalModule(ProjectModuleKind kind)
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    QString error;
    if (!project_.addOptionalModule(kind, &error)) {
        view_->showWarning(QStringLiteral("新增模块失败"), error);
        return;
    }
    refreshProjectView();
    log(QStringLiteral("可选模块已加入项目。"));
}

void GuiController::addField(ProjectModuleKind kind, const QString& fieldName)
{
    if (project_.directory().isEmpty()) return;
    QString error;
    if (!project_.addField(kind, fieldName, &error)) {
        view_->showWarning(QStringLiteral("新增字段失败"), error);
        return;
    }
    refreshProjectView();
    log(QStringLiteral("已新增字段 0/%1。").arg(fieldName));
}

void GuiController::generateMesh()
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    if (solverProcess_.state() != QProcess::NotRunning) {
        view_->showInformation(QStringLiteral("求解器正在运行"),
                               QStringLiteral("请先停止当前进程。"));
        return;
    }

    if (!project_.isSonicCase()) {
        const QString blockMeshDict = QDir(project_.directory()).filePath(
            QStringLiteral("system/blockMeshDict"));
        if (!QFileInfo(blockMeshDict).isFile()) {
            view_->showWarning(
                QStringLiteral("无法生成网格"),
                QStringLiteral("OpenFOAM 网格生成需要 system/blockMeshDict。"));
            return;
        }
        if (QStandardPaths::findExecutable(QStringLiteral("blockMesh")).isEmpty()) {
            view_->showWarning(
                QStringLiteral("找不到 blockMesh"),
                QStringLiteral("请在启动 sonicGui 前加载 OpenFOAM 环境。"));
            return;
        }
        saveAll();
        startSolverProcess(QStringLiteral("blockMesh"), {}, project_.directory(),
                           SolverProcessPurpose::MeshGeneration);
        return;
    }

    view_->showInformation(
        QStringLiteral("SFM 网格输入"),
        QStringLiteral(
            "sonicSolver case 使用已有的 constant/*.sfm 网格；请通过“文件 -> 导入 SFM 网格”导入，"
            "或先在专用 createMesh 工作流中生成。GUI 的 blockMesh 几何生成属于 OpenFOAM Mesh 工作域。"));
}

void GuiController::checkMesh()
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    if (project_.isSonicCase()) {
        view_->showInformation(
            QStringLiteral("SFM 网格"),
            QStringLiteral("sonicSolver case 使用 constant/*.sfm；请运行 createMesh 或直接检查 SFM 预览。"));
        return;
    }
    if (solverProcess_.state() != QProcess::NotRunning) {
        view_->showInformation(QStringLiteral("网格任务正在运行"),
                               QStringLiteral("请等待当前任务结束。"));
        return;
    }
    if (QStandardPaths::findExecutable(QStringLiteral("checkMesh")).isEmpty()) {
        view_->showWarning(
            QStringLiteral("找不到 checkMesh"),
            QStringLiteral("请在启动 sonicGui 前加载 OpenFOAM 环境。"));
        return;
    }
    startSolverProcess(QStringLiteral("checkMesh"), {}, project_.directory(),
                       SolverProcessPurpose::MeshCheck);
}

void GuiController::requestImportSfm()
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    if (!project_.isSonicCase()) {
        view_->showInformation(
            QStringLiteral("OpenFOAM 网格"),
            QStringLiteral("OpenFOAM case 从 constant/polyMesh 读取；不能导入 SFM。"));
        return;
    }
    const QString path = view_->chooseSfmFile(initialDirectory());
    if (path.isEmpty()) return;
    QString error;
    if (!project_.importMesh(path, &error)) {
        view_->showWarning(QStringLiteral("导入网格失败"), error);
        return;
    }
    importedSfmPath_ = QDir(project_.directory()).filePath(
        QStringLiteral("constant/mesh.sfm"));
    refreshProjectView();
    view_->showVisualizationView();
    if (!view_->loadVisualization(importedSfmPath_, &error)) {
        log(QStringLiteral("SFM 预览未更新：%1").arg(error));
    }
    log(QStringLiteral("已导入 SFM：%1；请根据 3D sets 预览设置 IC/BC 后再运行。")
            .arg(importedSfmPath_));
}

void GuiController::requestImportVtk()
{
    const QString path = view_->chooseVtkFile(initialDirectory());
    if (path.isEmpty()) return;
    QString error;
    if (!view_->loadVisualization(path, &error)) log(error);
}

void GuiController::requestImportStl()
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    const QString path = view_->chooseStlFile(initialDirectory());
    if (path.isEmpty()) return;
    QString error;
    if (!project_.importStl(path, &error)) {
        view_->showWarning(QStringLiteral("导入 STL 失败"), error);
        return;
    }
    refreshProjectView();
    const QString target = QDir(project_.directory()).filePath(
        QStringLiteral("constant/triSurface/%1")
            .arg(QFileInfo(path).fileName()));
    if (!view_->loadVisualization(target, &error)) log(error);
    log(QStringLiteral("IBM STL 已加入项目：%1").arg(target));
}

void GuiController::saveGeometry(const GeometrySketch& sketch)
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    QString error;
    if (!project_.writeGeometry(sketch, &error)) {
        view_->showWarning(QStringLiteral("几何参数未生成"), error);
        return;
    }
    refreshProjectView();
    log(QStringLiteral("几何草图已转换为单块挤出网格参数。"));
}

void GuiController::createOutputTask(const QString& taskName)
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    QString error;
    if (!project_.setOutputTaskName(taskName, &error)) {
        view_->showWarning(QStringLiteral("新建任务失败"), error);
        return;
    }
    refreshOutputTasks();
    log(QStringLiteral("当前输出任务：%1；结果将写入 result/%1.pvd。")
            .arg(project_.outputTaskName()));
}

void GuiController::submitOutputTask(const QString& taskName)
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    if (solverProcess_.state() != QProcess::NotRunning) {
        view_->showInformation(QStringLiteral("求解器正在运行"),
                               QStringLiteral("请先停止当前进程。"));
        return;
    }
    QString error;
    if (!project_.setOutputTaskName(taskName, &error)) {
        view_->showWarning(QStringLiteral("提交任务失败"), error);
        return;
    }
    refreshOutputTasks();
    runSolver(false);
}

void GuiController::selectOutputTask(const QString& taskName)
{
    if (project_.directory().isEmpty()) return;
    QString error;
    if (!project_.setOutputTaskName(taskName, &error)) {
        view_->showWarning(QStringLiteral("切换任务失败"), error);
        return;
    }
    refreshOutputTasks();

    const QString pvdPath =
        QDir(project_.resultDirectory())
            .filePath(project_.outputTaskName() + QStringLiteral(".pvd"));
    if (!QFileInfo(pvdPath).isFile()) {
        view_->showVisualizationView();
        log(QStringLiteral("任务 %1 尚无 PVD 结果。")
                .arg(project_.outputTaskName()));
        return;
    }
    if (!view_->loadVisualization(pvdPath, &error)) {
        log(error);
    } else {
        log(QStringLiteral("已切换输出任务：%1")
                .arg(project_.outputTaskName()));
    }
}

void GuiController::handleVisualizationLoaded(const QString& path)
{
    displayedOutputTask_.clear();
    const QFileInfo loaded(path);
    if (loaded.suffix().compare(QStringLiteral("pvd"),
                                Qt::CaseInsensitive) == 0
        && loaded.absolutePath()
               == QFileInfo(project_.resultDirectory()).absoluteFilePath()) {
        displayedOutputTask_ = loaded.completeBaseName();
    }
    refreshOutputTasks();
}

void GuiController::showOutputResidual(const QString& taskName)
{
    Q_UNUSED(taskName)
    log(QStringLiteral("残差显示窗口预留：后续接入求解器残差输出后绘制曲线。"));
}

void GuiController::refreshOutputTasks()
{
    view_->setOutputTasks(collectOutputTasks());
}

QVector<OutputTaskInfo> GuiController::collectOutputTasks() const
{
    QVector<OutputTaskInfo> tasks;
    if (project_.directory().isEmpty()) return tasks;

    const QString activeTask = project_.outputTaskName();
    QStringList taskNames = project_.outputTaskNames();
    QDir resultDir(project_.resultDirectory());
    const QFileInfoList pvdFiles = resultDir.entryInfoList(
        {QStringLiteral("*.pvd"), QStringLiteral("*.PVD")},
        QDir::Files,
        QDir::Name);
    for (const QFileInfo& pvd : pvdFiles) {
        if (pvd.fileName().startsWith(QStringLiteral("._"))) continue;
        const QString name = pvd.completeBaseName();
        if (!taskNames.contains(name, Qt::CaseInsensitive)) {
            taskNames.push_back(name);
        }
    }

    double startTime = 0.0;
    double endTime = 0.0;
    bool hasEndTime = false;
    const QString scPath = solverControlFile(project_.casePath());
    const QString scContent = readTextFile(scPath);
    bool ok = false;
    const QString startText = foamValue(scContent, QStringLiteral("startTime"));
    const double parsedStart = startText.toDouble(&ok);
    if (ok) startTime = parsedStart;
    const QString endText = foamValue(scContent, QStringLiteral("endTime"));
    const double parsedEnd = endText.toDouble(&ok);
    if (ok) {
        endTime = parsedEnd;
        hasEndTime = endTime > startTime;
    }

    for (const QString& name : taskNames) {
        OutputTaskInfo task;
        task.name = name;
        task.active = name.compare(activeTask, Qt::CaseInsensitive) == 0;
        task.displayed =
            name.compare(displayedOutputTask_, Qt::CaseInsensitive) == 0;
        const QString pvdPath =
            resultDir.filePath(name + QStringLiteral(".pvd"));
        if (QFileInfo(pvdPath).isFile()) {
            task.pvdPath = QFileInfo(pvdPath).absoluteFilePath();
            const PvdSummary summary = summarizePvd(task.pvdPath);
            task.frameCount = summary.frameCount;
            task.latestTime = summary.latestTimeText;
            task.latestTimeValue = summary.latestTime;
            if (hasEndTime && summary.hasLatestTime) {
                task.progress = std::clamp(
                    (summary.latestTime - startTime)
                        / (endTime - startTime),
                    0.0, 1.0);
            }
        }
        tasks.push_back(task);
    }
    return tasks;
}

QString GuiController::solverExecutable() const
{
    QString path = QDir(QCoreApplication::applicationDirPath())
                       .filePath(QStringLiteral("sonicSolver"));
    if (QFileInfo::exists(path)) return path;
    return QDir(sourceRoot_).filePath(QStringLiteral("build/sonicSolver"));
}

void GuiController::runSolver(bool previewOneStep)
{
    if (project_.casePath().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
        return;
    }
    if (solverProcess_.state() != QProcess::NotRunning) {
        view_->showInformation(QStringLiteral("求解器正在运行"),
                               QStringLiteral("请先停止当前进程。"));
        return;
    }

    saveAll();
    const CaseBackend backend = detectCaseBackend(project_.casePath());
    if (backend.kind == CaseBackendKind::Unsupported) {
        view_->showWarning(QStringLiteral("无法运行当前算例"), backend.error);
        return;
    }

    if (backend.kind == CaseBackendKind::OpenFOAMScript) {
        if (previewOneStep) {
            view_->showInformation(
                QStringLiteral("OpenFOAM 脚本算例不支持单步预览"),
                QStringLiteral("请使用“运行当前算例”调用 case 目录中的 Allrun/run 脚本。"));
            return;
        }
        startSolverProcess(QStringLiteral("/bin/bash"),
                           {backend.scriptPath},
                           project_.directory(),
                           SolverProcessPurpose::Solve);
        return;
    }

    QString runCase = project_.casePath();
    QString workingDirectory = project_.directory();
    if (previewOneStep) {
        QString error;
        runCase = preparePreviewCase(&error);
        if (runCase.isEmpty()) {
            view_->showWarning(QStringLiteral("预览准备失败"), error);
            return;
        }
        workingDirectory = QFileInfo(runCase).absolutePath();
    }

    const QString executable = solverExecutable();
    if (!QFileInfo::exists(executable)) {
        view_->showWarning(
            QStringLiteral("找不到求解器"),
            QStringLiteral("请启用 BUILD_CLI 并构建 sonicSolver：\n%1")
                .arg(executable));
        return;
    }

    QString program = executable;
    QStringList arguments;
    QStringList prefixArguments;
    QString parallelPath = parallelControlFile(runCase);
    if (parallelPath.isEmpty()) parallelPath = solverControlFile(runCase);
    if (parallelCommand(parallelPath, &program, &prefixArguments)) {
        arguments = prefixArguments;
        arguments << executable << runCase;
    } else {
        arguments << runCase;
    }
    startSolverProcess(
        program, arguments, workingDirectory,
        previewOneStep ? SolverProcessPurpose::PreviewStep
                       : SolverProcessPurpose::Solve);
}

void GuiController::startSolverProcess(const QString& program,
                                       const QStringList& arguments,
                                       const QString& workingDirectory,
                                       SolverProcessPurpose purpose)
{
    solverPurpose_ = purpose;
    view_->clearSolverOutput();
    solverProcess_.setWorkingDirectory(workingDirectory);
    solverProcess_.setProgram(program);
    solverProcess_.setArguments(arguments);
    solverProcess_.setProcessChannelMode(QProcess::SeparateChannels);
    solverProcess_.start();
    view_->setSolverRunning(true);
    log(QStringLiteral("启动：%1 %2")
            .arg(program, arguments.join(QLatin1Char(' '))));
}

void GuiController::appendSolverProcessOutput()
{
    const QString out =
        QString::fromLocal8Bit(solverProcess_.readAllStandardOutput());
    const QString err =
        QString::fromLocal8Bit(solverProcess_.readAllStandardError());
    view_->appendSolverOutput(out);
    view_->appendSolverOutput(err);
    const bool pvdUpdated = out.contains(QStringLiteral("PVD file updated:"))
        || err.contains(QStringLiteral("PVD file updated:"));
    if (pvdUpdated) {
        QTimer::singleShot(120, this, [this] { refreshOutputTasks(); });
    }
    if (solverPurpose_ == SolverProcessPurpose::PreviewStep
        && pvdUpdated) {
        QTimer::singleShot(250, this, [this] {
            if (solverProcess_.state() != QProcess::NotRunning) {
                solverProcess_.terminate();
            }
        });
    }
}

void GuiController::handleSolverFinished(int exitCode,
                                         QProcess::ExitStatus status)
{
    appendSolverProcessOutput();
    view_->setSolverRunning(false);
    const SolverProcessPurpose completedPurpose = solverPurpose_;
    solverPurpose_ = SolverProcessPurpose::Solve;
    log(QStringLiteral("求解器结束：exit=%1, status=%2")
            .arg(exitCode)
            .arg(status == QProcess::NormalExit
                     ? QStringLiteral("normal")
                     : QStringLiteral("crashed")));

    if (completedPurpose == SolverProcessPurpose::InitialOutput) {
        refreshOutputTasks();
        view_->refreshVisualization();
        if (exitCode == 0 && status == QProcess::NormalExit) {
            view_->showVisualizationView();
            log(QStringLiteral("已生成并加载未推进时间步的第 0 步结果。"));
        }
        return;
    }

    if (completedPurpose == SolverProcessPurpose::MeshGeneration) {
        if (exitCode != 0 || status != QProcess::NormalExit) {
            log(QStringLiteral("网格生成失败，请查看求解器输出。"));
            return;
        }
        QString error;
        if (!project_.isSonicCase()) {
            const QString polyMesh = openFoamPointerFile(project_.directory(),
                                                          &error);
            if (!polyMesh.isEmpty()) {
                if (!view_->loadVisualization(polyMesh, &error)) log(error);
            } else if (!error.isEmpty()) {
                log(error);
            }
            view_->showVisualizationView();
            log(QStringLiteral("blockMesh 已完成；已从 constant/polyMesh 刷新真实网格。"));
            return;
        }
        if (!project_.setMeshFile(generatedMeshName_, &error)) {
            view_->showWarning(QStringLiteral("网格已生成但清单更新失败"),
                               error);
            return;
        }
        refreshProjectView();
        const QDir caseDir(project_.directory());
        const QString sfmPath = caseDir.filePath(generatedMeshName_);
        const QString vtmPath =
            caseDir.filePath(QStringLiteral("constant/mesh.vtm"));
        const QString vtsPath =
            caseDir.filePath(QStringLiteral("constant/mesh.vts"));
        if (QFileInfo::exists(sfmPath)) {
            if (!view_->loadVisualization(sfmPath, &error)) log(error);
        } else if (QFileInfo::exists(vtmPath)) {
            if (!view_->loadVisualization(vtmPath, &error)) log(error);
        } else if (QFileInfo::exists(vtsPath)) {
            if (!view_->loadVisualization(vtsPath, &error)) log(error);
        }
        view_->showVisualizationView();
        log(QStringLiteral("网格已生成：%1")
                .arg(caseDir.filePath(generatedMeshName_)));
        return;
    }

    if (completedPurpose == SolverProcessPurpose::MeshCheck) {
        log(exitCode == 0 && status == QProcess::NormalExit
                ? QStringLiteral("checkMesh 完成。")
                : QStringLiteral("checkMesh 失败，请查看求解器输出。"));
        return;
    }

    if (completedPurpose == SolverProcessPurpose::Solve) {
        refreshOutputTasks();
        const QString pvdPath =
            QDir(project_.resultDirectory())
                .filePath(project_.outputTaskName() + QStringLiteral(".pvd"));
        if (exitCode == 0 && status == QProcess::NormalExit
            && QFileInfo(pvdPath).isFile()) {
            QString error;
            if (!view_->loadVisualization(pvdPath, &error)) log(error);
        }
        return;
    }

    if (completedPurpose != SolverProcessPurpose::PreviewStep) return;
    const QString output = latestPreviewOutput();
    if (output.isEmpty()) {
        log(QStringLiteral("单步预览未找到 VTS/VTM 输出。"));
        return;
    }
    QString error;
    if (!view_->loadVisualization(output, &error)) log(error);
    else log(QStringLiteral("单步预览输出：%1").arg(output));
}

void GuiController::stopSolver()
{
    if (solverProcess_.state() == QProcess::NotRunning) return;
    solverProcess_.terminate();
    QTimer::singleShot(1800, this, [this] {
        if (solverProcess_.state() != QProcess::NotRunning) {
            solverProcess_.kill();
        }
    });
}

void GuiController::runTerminalCommand(const QString& command)
{
    if (command.isEmpty()
        || terminalProcess_.state() != QProcess::NotRunning) {
        return;
    }
    terminalProcess_.setWorkingDirectory(
        project_.directory().isEmpty() ? sourceRoot_ : project_.directory());
    terminalProcess_.start(QStringLiteral("/bin/zsh"),
                           {QStringLiteral("-lc"), command});
}

QString GuiController::manifestFile(const QString& casePath,
                                    const QString& section) const
{
    const QDir caseDir(caseDirectoryFromPath(casePath));
    if (section.compare(QStringLiteral("solverControl"),
                        Qt::CaseInsensitive) == 0) {
        return caseDir.filePath(QStringLiteral("system/controlDict"));
    }
    if (section.compare(QStringLiteral("parallel"),
                        Qt::CaseInsensitive) == 0) {
        return caseDir.filePath(QStringLiteral("system/sonicDict"));
    }
    return {};
}

QString GuiController::solverControlFile(const QString& casePath) const
{
    return manifestFile(casePath, QStringLiteral("solverControl"));
}

QString GuiController::parallelControlFile(const QString& casePath) const
{
    return manifestFile(casePath, QStringLiteral("parallel"));
}

bool GuiController::parallelCommand(const QString& configPath,
                                    QString* program,
                                    QStringList* prefixArguments) const
{
    if (configPath.isEmpty()) return false;
    const QString content = readTextFile(configPath);
    const QString enabled = foamBlockValue(content, QStringLiteral("parallel"),
                                           QStringLiteral("enabled"));
    if (enabled.compare(QStringLiteral("true"),
                        Qt::CaseInsensitive) != 0) return false;

    int nProcs = 1;
    const QStringList split =
        foamWordList(foamBlockValue(content, QStringLiteral("parallel"),
                                    QStringLiteral("split")));
    if (split.size() == 3) {
        bool ok0 = false;
        bool ok1 = false;
        bool ok2 = false;
        const int sx = split.at(0).toInt(&ok0);
        const int sy = split.at(1).toInt(&ok1);
        const int sz = split.at(2).toInt(&ok2);
        if (ok0 && ok1 && ok2) nProcs = sx * sy * sz;
    } else {
        bool ok = false;
        const int parsed =
            foamBlockValue(content, QStringLiteral("parallel"),
                           QStringLiteral("nProcs")).toInt(&ok);
        if (ok) nProcs = parsed;
    }
    if (nProcs <= 1) return false;
    *program = QStringLiteral("mpirun");
    *prefixArguments = {QStringLiteral("-np"), QString::number(nProcs)};
    return true;
}

QString GuiController::preparePreviewCase(QString* error)
{
    const QString base =
        QStandardPaths::writableLocation(QStandardPaths::TempLocation);
    previewDirectory_ = QDir(base).filePath(
        QStringLiteral("sonicGuiPreview/%1_%2")
            .arg(QFileInfo(project_.directory()).fileName(),
                 QDateTime::currentDateTime().toString(
                     QStringLiteral("yyyyMMdd_hhmmss_zzz"))));
    if (!copyDirectory(project_.directory(), previewDirectory_, error)) return {};

    const QString previewCase = previewDirectory_;
    const QString scPath = solverControlFile(previewCase);
    if (scPath.isEmpty()) {
        if (error) *error = QStringLiteral("未找到 system/controlDict");
        return {};
    }
    if (!patchPreviewControl(scPath, error)) return {};
    if (!importedSfmPath_.isEmpty()
        && !patchPreviewMesh(previewCase, error)) {
        return {};
    }
    return previewCase;
}

bool GuiController::copyDirectory(const QString& source,
                                  const QString& destination,
                                  QString* error)
{
    QDir sourceDir(source);
    if (!sourceDir.exists() || !QDir().mkpath(destination)) {
        if (error) *error = QStringLiteral("无法创建预览目录");
        return false;
    }
    const QFileInfoList entries = sourceDir.entryInfoList(
        QDir::NoDotAndDotDot | QDir::Files | QDir::Dirs);
    for (const QFileInfo& entry : entries) {
        if (entry.fileName() == QStringLiteral("result")) continue;
        const QString target = QDir(destination).filePath(entry.fileName());
        if (entry.isDir()) {
            if (!copyDirectory(entry.absoluteFilePath(), target, error)) {
                return false;
            }
        } else if (!QFile::copy(entry.absoluteFilePath(), target)) {
            if (error) *error = QStringLiteral("无法复制：%1")
                .arg(entry.absoluteFilePath());
            return false;
        }
    }
    return true;
}

bool GuiController::replaceSectionValue(QStringList* lines,
                                        const QString& section,
                                        const QString& key,
                                        const QString& value)
{
    if (section.isEmpty()) {
        const QRegularExpression foamKeyPattern(
            QStringLiteral(R"(^(\s*)%1\s+.*;\s*(?://.*)?$)")
                .arg(QRegularExpression::escape(key)),
            QRegularExpression::CaseInsensitiveOption);
        for (int i = 0; i < lines->size(); ++i) {
            const QRegularExpressionMatch match =
                foamKeyPattern.match(lines->at(i));
            if (!match.hasMatch()) continue;
            (*lines)[i] = match.captured(1) + key + QLatin1Char(' ')
                        + value + QLatin1Char(';');
            return true;
        }
        lines->append(key + QLatin1Char(' ') + value + QLatin1Char(';'));
        return true;
    }

    const QRegularExpression sectionPattern(
        QStringLiteral(R"(^\s*\[([^\]]+)\])"));
    const QRegularExpression keyPattern(
        QStringLiteral(R"(^(\s*)%1\s*=.*$)")
            .arg(QRegularExpression::escape(key)),
        QRegularExpression::CaseInsensitiveOption);
    bool inSection = false;
    int insertionLine = -1;
    for (int i = 0; i < lines->size(); ++i) {
        const QRegularExpressionMatch sectionMatch =
            sectionPattern.match(lines->at(i));
        if (sectionMatch.hasMatch()) {
            inSection = sectionMatch.captured(1).trimmed().compare(
                            section, Qt::CaseInsensitive) == 0;
            insertionLine = inSection ? i + 1 : -1;
            continue;
        }
        if (!inSection) continue;
        const QRegularExpressionMatch keyMatch =
            keyPattern.match(lines->at(i));
        if (keyMatch.hasMatch()) {
            (*lines)[i] = keyMatch.captured(1) + key
                        + QStringLiteral(" = ") + value;
            return true;
        }
        insertionLine = i + 1;
    }
    if (insertionLine >= 0) {
        lines->insert(insertionLine,
                      QStringLiteral("    %1 = %2").arg(key, value));
        return true;
    }
    lines->append(QStringLiteral("[%1]").arg(section));
    lines->append(QStringLiteral("    %1 = %2").arg(key, value));
    return true;
}

bool GuiController::patchPreviewControl(const QString& scPath,
                                        QString* error) const
{
    const QString content = readTextFile(scPath, error);
    if (content.isNull()) return false;
    QStringList lines =
        content.split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    replaceSectionValue(&lines, QString(), QStringLiteral("writeControl"),
                        QStringLiteral("timeStep"));
    replaceSectionValue(&lines, QString(), QStringLiteral("writeInterval"),
                        QStringLiteral("1"));
    return writeTextFile(scPath, lines, error);
}

bool GuiController::patchPreviewMesh(const QString& casePath,
                                     QString* error) const
{
    const QString sonicPath =
        QDir(caseDirectoryFromPath(casePath)).filePath(
            QStringLiteral("system/sonicDict"));
    const QString content = readTextFile(sonicPath, error);
    if (content.isNull()) return false;
    QStringList lines =
        content.split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    for (int i = 0; i < lines.size(); ++i) {
        const QString trimmed = lines.at(i).trimmed();
        if (trimmed.startsWith(QStringLiteral("meshFiles"),
                               Qt::CaseInsensitive)) {
            lines[i] = QStringLiteral("meshFiles (constant/mesh.sfm);");
            return writeTextFile(sonicPath, lines, error);
        }
    }
    if (error) *error = QStringLiteral("system/sonicDict 没有 meshFiles 配置");
    return false;
}

QString GuiController::latestPreviewOutput() const
{
    if (previewDirectory_.isEmpty()) return {};
    QFileInfoList candidates;
    QDirIterator iterator(
        previewDirectory_,
        {QStringLiteral("*.vts"), QStringLiteral("*.vtm"),
         QStringLiteral("*.vtu")},
        QDir::Files,
        QDirIterator::Subdirectories);
    while (iterator.hasNext()) candidates << QFileInfo(iterator.next());
    if (candidates.isEmpty()) return {};
    std::sort(candidates.begin(), candidates.end(),
              [](const QFileInfo& a, const QFileInfo& b) {
        return a.lastModified() > b.lastModified();
    });
    return candidates.first().absoluteFilePath();
}

void GuiController::log(const QString& message)
{
    view_->appendLog(
        QStringLiteral("[%1] %2")
            .arg(QDateTime::currentDateTime().toString(
                     QStringLiteral("HH:mm:ss")),
                 message));
}

} // namespace SF::GUI
