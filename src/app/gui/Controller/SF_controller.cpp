#include <QDebug>
#include "app/application/model/SF_model.h"
#include "infrastructure/io/case/SF_case.h"
#include "infrastructure/io/serialization/SF_serialization.h"
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

struct PvdSummary {
    int frameCount = 0;
    QString latestTimeText;
    double latestTime = 0.0;
    bool hasLatestTime = false;
};

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
        const auto result=qEnvironmentVariable("SONIC_GUI_TEST_RESULT");
        if(!result.isEmpty()) QTimer::singleShot(950,this,[this,result]{QString error;view_->showVisualizationView();if(!view_->loadVisualization(result,&error))log(error);});
        const auto action=qEnvironmentVariable("SONIC_GUI_TEST_ACTION");
        if(action=="preview" || action=="run" || action=="mesh") QTimer::singleShot(900,this,[this,action]{if(action=="mesh")generateMesh();else runSolver(action=="preview");});
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
                if(qEnvironmentVariableIsSet("SONIC_GUI_TEST_ACTION")) QCoreApplication::exit(3);
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

void GuiController::openWorkspace(const QString& directory) {
    if(solverProcess_.state()!=QProcess::NotRunning){view_->showWarning("Busy","Stop the active process first");return;}
    QString error;if(!project_.open(directory,&error)){view_->showWarning("Native project",error);return;}
    importedSfmPath_.clear();displayedOutputTask_.clear();refreshProjectView();
    log("Native case: "+project_.directory());
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

void GuiController::generateInitialOutput() {
    if(solverProcess_.state()!=QProcess::NotRunning)return;
    QString error;if(!canGenerateInitialOutput(&error)){log(error);return;}
    QString program=solverExecutable();QStringList args,prefix;
    try {if(parallelCommand(project_.casePath(),&program,&prefix))args=prefix+QStringList{solverExecutable()};
        args<<"initialOutput"<<project_.casePath();startSolverProcess(program,args,project_.directory(),SolverProcessPurpose::InitialOutput);
    }catch(const std::exception& e){log(QString::fromUtf8(e.what()));}
}

bool GuiController::canGenerateInitialOutput(QString* reason) const {
    try {auto c=Application::ModelLoader::read(project_.casePath().toStdString());
        if(c.createMesh)throw std::runtime_error("Generate the registered native mesh first");
        for(const auto& file:c.meshFiles)if(!QFileInfo(QDir(project_.directory()).filePath(QString::fromStdString(file))).isFile())throw std::runtime_error("Registered mesh payload is missing");return true;
    }catch(const std::exception& e){if(reason)*reason=QString::fromUtf8(e.what());return false;}
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



void GuiController::generateMesh() {
    if(project_.directory().isEmpty())return;
    try {auto model=CaseIO::read(project_.directory().toStdString());if(model.meshGeneration.empty())throw std::runtime_error("Register structuredMesh generation first");
        // A distinct native copy controls mesh generation; the user's runtime is not patched.
        QString error;const QString copy=preparePreviewCase(&error);if(copy.isEmpty())throw std::runtime_error(error.toStdString());
        model=CaseIO::read(copy.toStdString());model.runtime["createMesh"]=true;CaseIO::write(model,copy.toStdString());
        startSolverProcess(solverExecutable(),{"run",copy},copy,SolverProcessPurpose::MeshGeneration);
    }catch(const std::exception& e){view_->showWarning("Native mesh",QString::fromUtf8(e.what()));}
}

void GuiController::checkMesh() {
    if(!project_.directory().isEmpty())startSolverProcess(solverExecutable(),{"check",project_.casePath()},project_.directory(),SolverProcessPurpose::MeshCheck);
}

void GuiController::requestImportSfm()
{
    if (project_.directory().isEmpty()) {
        view_->showInformation(QStringLiteral("尚未设置工作文件夹"),
                               QStringLiteral("请先设置工作文件夹。"));
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
        QStringLiteral("mesh/")+QFileInfo(path).fileName());
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
        QStringLiteral("mesh/geometry/%1")
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

    const auto config=Application::ModelLoader::read(project_.casePath().toStdString());
    const double startTime=config.time.startTime,endTime=config.time.endTime;
    const bool hasEndTime=endTime>startTime;

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

void GuiController::runSolver(bool preview) {
    if(project_.casePath().isEmpty() || solverProcess_.state()!=QProcess::NotRunning)return;
    QString error;if(!project_.saveAll(&error)){view_->showWarning("Native validation",error);return;}
    try {
        const QString path=preview?preparePreviewCase(&error):project_.casePath();if(path.isEmpty())throw std::runtime_error(error.toStdString());
        QString program=solverExecutable();QStringList prefix,args;
        if(parallelCommand(path,&program,&prefix))args=prefix+QStringList{solverExecutable()};
        if(preview)args<<"steps"<<"1"<<path;else args<<"run"<<path;
        startSolverProcess(program,args,path,preview?SolverProcessPurpose::PreviewStep:SolverProcessPurpose::Solve);
    }catch(const std::exception& e){view_->showWarning("Native run",QString::fromUtf8(e.what()));}
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
    if(qEnvironmentVariableIsSet("SONIC_GUI_TEST_ACTION"))qInfo().noquote()<<program<<arguments;
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
    if(!qEnvironmentVariable("SONIC_GUI_TEST_ACTION").isEmpty())qInfo().noquote()<<out<<err;
    view_->appendSolverOutput(out);
    view_->appendSolverOutput(err);
    const bool pvdUpdated = out.contains(QStringLiteral("PVD file updated:"))
        || err.contains(QStringLiteral("PVD file updated:"));
    if (pvdUpdated) {
        QTimer::singleShot(120, this, [this] { refreshOutputTasks(); });
    }

}

void GuiController::handleSolverFinished(int code,QProcess::ExitStatus status) {
    appendSolverProcessOutput();view_->setSolverRunning(false);refreshOutputTasks();
    if(code!=0 || status!=QProcess::NormalExit){log("Process failed: "+QString::number(code));if(!qEnvironmentVariable("SONIC_GUI_TEST_ACTION").isEmpty())QCoreApplication::exit(code?code:1);return;}
    QString error;
    if(solverPurpose_==SolverProcessPurpose::MeshGeneration) {
        const auto input=CaseIO::read(previewDirectory_.toStdString());
        for(const auto& file:input.mesh.at("files")) {const auto relative=QString::fromStdString(file.get<std::string>());const auto source=QDir(previewDirectory_).filePath(relative);if(!project_.importMesh(source,&error,true)){log(error);return;}}
        refreshProjectView();
    } else if(solverPurpose_==SolverProcessPurpose::PreviewStep) {
        const auto output=latestPreviewOutput();if(!output.isEmpty()&&!view_->loadVisualization(output,&error))log(error);
    } else view_->refreshVisualization();
    log("Native process completed normally");
    if(!qEnvironmentVariable("SONIC_GUI_TEST_ACTION").isEmpty())QCoreApplication::exit(0);
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







bool GuiController::parallelCommand(const QString& path,QString* program,QStringList* prefix) const {
    const auto config=Application::ModelLoader::read(path.toStdString());
    if(!config.parallel.enabled)return false;
    if(config.parallel.processCount<1)throw std::runtime_error("Invalid native parallel rank count");
    *program="mpirun";*prefix={"-np",QString::number(config.parallel.processCount)};return true;
}

QString GuiController::preparePreviewCase(QString* error) {
    previewDirectory_=QDir(sourceRoot_).filePath("test/guiNativePreview/"+QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss_zzz"));
    if(!copyDirectory(project_.directory(),previewDirectory_,error))return {};
    try {
        auto model=CaseIO::read(previewDirectory_.toStdString());model.runtime["writeControl"]="timeStep";model.runtime["writeInterval"]=1;
        model.runtime["writeInitial"]=true;model.runtime["output"]["outputDir"]="result";
        CaseIO::write(model,previewDirectory_.toStdString());(void)Application::ModelLoader::read(previewDirectory_.toStdString());return previewDirectory_;
    }catch(const std::exception& e){if(error)*error=QString::fromUtf8(e.what());return {};}
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
        if (entry.fileName() == QStringLiteral("result") || entry.fileName().startsWith(".")) continue;
        if(entry.isSymLink()){if(error)*error="Preview refuses symlink payloads";return false;}
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
