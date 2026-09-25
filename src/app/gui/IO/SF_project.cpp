/// @file SF_project.cpp
/// @brief GUI 工程与中立配置文档的读写实现。

#include "SF_project.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QRegularExpression>
#include <QSet>

#include <algorithm>

namespace SF::GUI {
namespace {

QString projectName(const QString& directory)
{
    QString name = QFileInfo(directory).fileName();
    name.replace(QRegularExpression(QStringLiteral("[^A-Za-z0-9_-]")),
                 QStringLiteral("_"));
    return name.isEmpty() ? QStringLiteral("sonicCase") : name;
}

QString controlDictTemplate(const QString& name)
{
    Q_UNUSED(name);
    return QStringLiteral(
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object controlDict;\n"
        "}\n"
        "application sonicSolver;\n"
        "createMesh false;\n"
        "startTime 0;\n"
        "endTime 1;\n"
        "CFL 0.5;\n"
        "writeControl runTime;\n"
        "writeInterval 0.1;\n");
}

QString outputDictTemplate(const QString& name)
{
    return QStringLiteral(
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object outputDict;\n"
        "}\n"
        "jobName %1;\n"
        "jobs (%1);\n"
        "outputDir result;\n").arg(name);
}

QString sonicDictTemplate(const QString& name)
{
    Q_UNUSED(name);
    return QStringLiteral(
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object sonicDict;\n"
        "}\n");
}

QStringList sfmCompanionSetCandidates(const QFileInfo& source)
{
    QStringList candidates;
    Q_UNUSED(source);
    return candidates;
}

const char* geometryTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object blockMeshDict;\n"
        "}\n\n"
        "convertToMeters 1;\n\n"
        "vertices\n(\n);\n\n"
        "blocks\n(\n);\n\n"
        "edges\n(\n);\n\n"
        "boundary\n(\n);\n\n"
        "mergePatchPairs\n(\n);\n";
}

const char* solverPropertiesTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object solverProperties;\n"
        "}\n"
        "type densityBase;\n";
}

const char* turbulenceTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object turbulenceProperties;\n"
        "}\n"
        "simulationType laminar;\n\n"
        "RAS\n"
        "{\n"
        "    RASModel kOmegaSST;\n"
        "    turbulence on;\n"
        "    printCoeffs on;\n"
        "}\n";
}

const char* multiphaseTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object phaseProperties;\n"
        "}\n"
        "phases ();\n";
}

const char* phaseChangeTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object phaseChange;\n"
        "}\n"
        "enabled false;\n"
        "model none;\n";
}

const char* initialConditionTemplate()
{
    return "";
}

const char* boundaryConditionTemplate()
{
    return "";
}

const char* ibmTemplate()
{
    return
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object IBMProperties;\n"
        "}\n"
        "IBM\n"
        "{\n"
        "    type ghost;\n"
        "    geometryFiles ();\n"
        "}\n";
}

const char* mrfTemplate()
{
    return "";
}

const char* parallelTemplate()
{
    return "";
}

QString stripFoamComments(QString text)
{
    text.remove(QRegularExpression(QStringLiteral(R"(/\*[\s\S]*?\*/)")));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\n]*)")));
    return text;
}

QString readFile(const QString& path)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return {};
    return QString::fromUtf8(file.readAll());
}

QStringList foamBoundaryNames(const QString& path)
{
    const QStringList lines = stripFoamComments(readFile(path))
                                  .split(QLatin1Char('\n'));
    QStringList names;
    const QRegularExpression namePattern(
        QStringLiteral(R"(^\s*([A-Za-z_][A-Za-z0-9_.-]*)\s*$)"));
    for (int i = 0; i + 1 < lines.size(); ++i) {
        const QRegularExpressionMatch match = namePattern.match(lines.at(i));
        if (!match.hasMatch()) continue;
        int next = i + 1;
        while (next < lines.size() && lines.at(next).trimmed().isEmpty()) ++next;
        if (next < lines.size() && lines.at(next).trimmed().startsWith(QLatin1Char('{'))) {
            const QString name = match.captured(1);
            if (name.compare(QStringLiteral("FoamFile"), Qt::CaseInsensitive) != 0
                && !names.contains(name)) names.push_back(name);
        }
    }
    return names;
}

QStringList sfmSetNames(const QString& path)
{
    QStringList names;
    const QStringList lines = readFile(path).split(QLatin1Char('\n'));
    const QSet<QString> metadata{QStringLiteral("information"),
                                 QStringLiteral("point"),
                                 QStringLiteral("cell"),
                                 QStringLiteral("face"),
                                 QStringLiteral("patch"),
                                 QStringLiteral("pointset")};
    for (const QString& line : lines) {
        const QString trimmed = line.trimmed();
        if (!trimmed.startsWith(QLatin1Char('#'))) continue;
        const QString name = trimmed.mid(1).trimmed();
        if (name.isEmpty() || metadata.contains(name.toLower())) continue;
        if (!names.contains(name)) names.push_back(name);
    }
    return names;
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

bool replaceFoamValueInFile(const QString& path,
                            const QString& key,
                            const QString& value,
                            QString* error)
{
    auto writeLines = [&](const QStringList& lines) {
        QFile output(path);
        if (!output.open(QIODevice::WriteOnly | QIODevice::Text
                         | QIODevice::Truncate)) {
            if (error) *error = output.errorString();
            return false;
        }
        QByteArray data = lines.join(QLatin1Char('\n')).toUtf8();
        data.append('\n');
        if (output.write(data) != data.size()) {
            if (error) *error = output.errorString();
            return false;
        }
        return true;
    };

    QFile input(path);
    if (!input.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = input.errorString();
        return false;
    }
    QStringList lines = QString::fromUtf8(input.readAll())
                            .split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    const QRegularExpression keyPattern(
        QStringLiteral(R"(^(\s*)%1\s+.*;\s*(?://.*)?$)")
            .arg(QRegularExpression::escape(key)),
        QRegularExpression::CaseInsensitiveOption);
    for (int i = 0; i < lines.size(); ++i) {
        const QRegularExpressionMatch match = keyPattern.match(lines.at(i));
        if (!match.hasMatch()) continue;
        lines[i] = match.captured(1) + key + QLatin1Char(' ') + value
                 + QLatin1Char(';');
        return writeLines(lines);
    }
    if (!lines.isEmpty() && !lines.last().trimmed().isEmpty()) {
        lines.push_back(QString());
    }
    lines.push_back(key + QLatin1Char(' ') + value + QLatin1Char(';'));
    return writeLines(lines);
}

bool replaceFoamBlockValueInFile(const QString& path,
                                 const QString& block,
                                 const QString& key,
                                 const QString& value,
                                 QString* error)
{
    auto writeLines = [&](const QStringList& outputLines) {
        QFile output(path);
        if (!output.open(QIODevice::WriteOnly | QIODevice::Text
                         | QIODevice::Truncate)) {
            if (error) *error = output.errorString();
            return false;
        }
        QByteArray data = outputLines.join(QLatin1Char('\n')).toUtf8();
        data.append('\n');
        if (output.write(data) != data.size()) {
            if (error) *error = output.errorString();
            return false;
        }
        return true;
    };

    QFile input(path);
    if (!input.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = input.errorString();
        return false;
    }

    QStringList lines = QString::fromUtf8(input.readAll())
                            .split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    if (!lines.isEmpty() && lines.last().isEmpty()) lines.removeLast();

    const QRegularExpression keyPattern(
        QStringLiteral(R"(^(\s*)%1\s+.*;\s*(?://.*)?$)")
            .arg(QRegularExpression::escape(key)),
        QRegularExpression::CaseInsensitiveOption);
    bool waitingBrace = false;
    bool inBlock = false;
    int insertion = -1;
    for (int i = 0; i < lines.size(); ++i) {
        const QString trimmed = lines.at(i).trimmed();
        if (!inBlock) {
            if (waitingBrace && trimmed == QStringLiteral("{")) {
                inBlock = true;
                waitingBrace = false;
                insertion = i + 1;
                continue;
            }
            if (trimmed.compare(block, Qt::CaseInsensitive) == 0) {
                waitingBrace = true;
                continue;
            }
            if (trimmed.startsWith(block + QStringLiteral("{"),
                                   Qt::CaseInsensitive)) {
                inBlock = true;
                insertion = i + 1;
                continue;
            }
            continue;
        }
        if (trimmed.startsWith(QLatin1Char('}'))) break;
        const QRegularExpressionMatch match = keyPattern.match(lines.at(i));
        if (match.hasMatch()) {
            lines[i] = match.captured(1) + key + QLatin1Char(' ')
                     + value + QLatin1Char(';');
            return writeLines(lines);
        }
        if (!trimmed.isEmpty()) insertion = i + 1;
    }

    if (insertion < 0) {
        if (!lines.isEmpty() && !lines.last().trimmed().isEmpty()) {
            lines.push_back(QString());
        }
        lines.push_back(block);
        lines.push_back(QStringLiteral("{"));
        lines.push_back(QStringLiteral("    %1 %2;").arg(key, value));
        lines.push_back(QStringLiteral("}"));
    } else {
        lines.insert(insertion,
                     QStringLiteral("    %1 %2;").arg(key, value));
    }
    return writeLines(lines);
}

bool copyFileToPath(const QString& sourcePath,
                    const QString& targetPath,
                    QString* error)
{
    const QFileInfo source(sourcePath);
    if (!source.isFile()) {
        if (error) *error = QStringLiteral("文件不存在：%1").arg(sourcePath);
        return false;
    }
    const QFileInfo targetInfo(targetPath);
    if (!QDir().mkpath(targetInfo.absolutePath())) {
        if (error) {
            *error = QStringLiteral("无法创建目录：%1")
                .arg(targetInfo.absolutePath());
        }
        return false;
    }
    if (source.absoluteFilePath() == targetInfo.absoluteFilePath()) return true;
    if (targetInfo.exists() && !QFile::remove(targetInfo.absoluteFilePath())) {
        if (error) {
            *error = QStringLiteral("无法替换：%1")
                .arg(targetInfo.absoluteFilePath());
        }
        return false;
    }
    if (!QFile::copy(source.absoluteFilePath(), targetInfo.absoluteFilePath())) {
        if (error) {
            *error = QStringLiteral("无法复制到：%1")
                .arg(targetInfo.absoluteFilePath());
        }
        return false;
    }
    return true;
}

QString templateFor(ProjectModuleKind kind)
{
    switch (kind) {
        case ProjectModuleKind::Geometry: return QString::fromUtf8(geometryTemplate());
        case ProjectModuleKind::Solver:
            return controlDictTemplate(QStringLiteral("sonicCase"));
        case ProjectModuleKind::Turbulence: return QString::fromUtf8(turbulenceTemplate());
        case ProjectModuleKind::Multiphase: return QString::fromUtf8(multiphaseTemplate());
        case ProjectModuleKind::PhaseChange: return QString::fromUtf8(phaseChangeTemplate());
        case ProjectModuleKind::InitialCondition:
            return QString::fromUtf8(initialConditionTemplate());
        case ProjectModuleKind::BoundaryCondition:
            return QString::fromUtf8(boundaryConditionTemplate());
        case ProjectModuleKind::IBM: return QString::fromUtf8(ibmTemplate());
        case ProjectModuleKind::MRF: return QString::fromUtf8(mrfTemplate());
        case ProjectModuleKind::Parallel: return QString::fromUtf8(parallelTemplate());
        default: return {};
    }
}

} // namespace

bool CaseProject::open(const QString& directory, QString* error)
{
    const QFileInfo info(directory);
    directory_ = info.absoluteFilePath();
    if (!QDir().mkpath(directory_)) {
        if (error) *error = QStringLiteral("无法创建工作文件夹：%1").arg(directory_);
        return false;
    }
    casePath_ = directory_;
    const QDir caseDir(directory_);
    const bool hasCaseFiles = QFileInfo(caseDir.filePath(
            QStringLiteral("case.yaml"))).isFile()
        || QFileInfo(caseDir.filePath(
            QStringLiteral("system/controlDict"))).isFile()
        || QFileInfo(caseDir.filePath(
            QStringLiteral("system/outputDict"))).isFile()
        || QFileInfo(caseDir.filePath(
            QStringLiteral("system/sonicDict"))).isFile()
        || QFileInfo(caseDir.filePath(
            QStringLiteral("system/blockMeshDict"))).isFile()
        || QFileInfo(caseDir.filePath(
            QStringLiteral("constant/polyMesh"))).isDir();
    if (!hasCaseFiles && !createDefaultProject(error)) {
        return false;
    }
    if (!reload(error)) return false;
    if (!QDir().mkpath(resultDirectory())) {
        if (error) {
            *error = QStringLiteral("无法创建结果目录：%1")
                .arg(resultDirectory());
        }
        return false;
    }
    return true;
}

bool CaseProject::createDefaultProject(QString* error)
{
    const QDir dir(directory_);
    if (!QDir().mkpath(dir.filePath(QStringLiteral("system")))
        || !QDir().mkpath(dir.filePath(QStringLiteral("constant")))
        || !QDir().mkpath(dir.filePath(QStringLiteral("constant/triSurface")))
        || !QDir().mkpath(dir.filePath(QStringLiteral("0")))
        || !QDir().mkpath(dir.filePath(QStringLiteral("result")))) {
        if (error) *error = QStringLiteral("无法创建 OpenFOAM-like case 目录。");
        return false;
    }

    const QString name = projectName(directory_);
    const QString controlPath = dir.filePath(QStringLiteral("system/controlDict"));
    if (!QFileInfo::exists(controlPath)
        && !writeFile(controlPath, controlDictTemplate(name), error)) {
        return false;
    }
    const QString outputPath = dir.filePath(QStringLiteral("system/outputDict"));
    if (!QFileInfo::exists(outputPath)
        && !writeFile(outputPath, outputDictTemplate(name), error)) {
        return false;
    }
    const QString sonicPath = dir.filePath(QStringLiteral("system/sonicDict"));
    if (!QFileInfo::exists(sonicPath)
        && !writeFile(sonicPath, sonicDictTemplate(name), error)) {
        return false;
    }
    const QString solverPath =
        dir.filePath(QStringLiteral("system/solverProperties"));
    if (!QFileInfo::exists(solverPath)
        && !writeFile(solverPath,
                      QString::fromUtf8(solverPropertiesTemplate()), error)) {
        return false;
    }
    const QString turbulencePath =
        dir.filePath(QStringLiteral("constant/turbulenceProperties"));
    if (!QFileInfo::exists(turbulencePath)
        && !writeFile(turbulencePath, templateFor(ProjectModuleKind::Turbulence),
                      error)) {
        return false;
    }
    const QString geometryPath = dir.filePath(defaultFileName(ProjectModuleKind::Geometry));
    if (!QFileInfo::exists(geometryPath)
        && !writeFile(geometryPath, templateFor(ProjectModuleKind::Geometry),
                      error)) {
        return false;
    }
    return true;
}

bool CaseProject::reload(QString* error)
{
    documents_.clear();
    modules_.clear();

    const QDir caseDir(directory_);
    const auto addEmptyModule = [this](ProjectModuleKind kind) {
        ProjectModule module;
        module.kind = kind;
        module.title = moduleTitle(kind);
        module.documentIndex = -1;
        modules_.push_back(module);
    };
    const auto loadExisting = [this, error](const QString& path,
                                            ProjectModuleKind kind) {
        return loadDocument(path, kind, moduleTitle(kind), false, error);
    };

    // 原生 registry case 由 solver 的严格 YAML reader 管理。GUI 先提供模块导航、
    // 网格预览和运行入口，不用旧模板补写 system/constant/0 目录。
    if (QFileInfo(caseDir.filePath(QStringLiteral("case.yaml"))).isFile()) {
        for (const auto kind : {ProjectModuleKind::Geometry,
                 ProjectModuleKind::Mesh, ProjectModuleKind::Solver,
                 ProjectModuleKind::Turbulence, ProjectModuleKind::Multiphase,
                 ProjectModuleKind::PhaseChange,
                 ProjectModuleKind::InitialCondition,
                 ProjectModuleKind::BoundaryCondition, ProjectModuleKind::Output,
                 ProjectModuleKind::IBM, ProjectModuleKind::MRF,
                 ProjectModuleKind::Parallel}) {
            addEmptyModule(kind);
        }
        return true;
    }

    const QString blockMeshDict =
        caseDir.filePath(QStringLiteral("system/blockMeshDict"));
    if (QFileInfo(blockMeshDict).isFile()) {
        if (!loadExisting(blockMeshDict, ProjectModuleKind::Geometry)) return false;
    } else {
        addEmptyModule(ProjectModuleKind::Geometry);
    }

    const bool sonicCase = isSonicCase();
    const QString sonicDict = caseDir.filePath(QStringLiteral("system/sonicDict"));
    const QString snappyDict = caseDir.filePath(
        QStringLiteral("system/snappyHexMeshDict"));
    if (sonicCase) {
        if (!loadExisting(sonicDict, ProjectModuleKind::Mesh)) return false;
    } else if (QFileInfo(snappyDict).isFile()) {
        if (!loadExisting(snappyDict, ProjectModuleKind::Mesh)) return false;
    } else {
        addEmptyModule(ProjectModuleKind::Mesh);
    }

    const QString controlDict = caseDir.filePath(QStringLiteral("system/controlDict"));
    if (QFileInfo(controlDict).isFile()) {
        if (!loadExisting(controlDict, ProjectModuleKind::Solver)) return false;
    } else {
        addEmptyModule(ProjectModuleKind::Solver);
    }

    const QString turbulence = caseDir.filePath(
        QStringLiteral("constant/turbulenceProperties"));
    if (QFileInfo(turbulence).isFile()
        && !loadExisting(turbulence, ProjectModuleKind::Turbulence)) {
        return false;
    }

    const QString phaseProperties = caseDir.filePath(
        QStringLiteral("constant/phaseProperties"));
    if (QFileInfo(phaseProperties).isFile()) {
        if (!loadExisting(phaseProperties, ProjectModuleKind::Multiphase)) return false;
    } else {
        addEmptyModule(ProjectModuleKind::Multiphase);
    }

    const QString phaseChange = caseDir.filePath(
        QStringLiteral("constant/phaseChange"));
    if (QFileInfo(phaseChange).isFile()
        && !loadDocument(phaseChange, ProjectModuleKind::PhaseChange,
                         moduleTitle(ProjectModuleKind::PhaseChange), true, error)) {
        return false;
    }

    const QDir zeroDir(QDir(directory_).filePath(QStringLiteral("0")));
    ProjectModule initialModule;
    initialModule.kind = ProjectModuleKind::InitialCondition;
    initialModule.title = moduleTitle(ProjectModuleKind::InitialCondition);
    initialModule.densityBased = isDensityBased();
    ProjectModule boundaryModule;
    boundaryModule.kind = ProjectModuleKind::BoundaryCondition;
    boundaryModule.title = moduleTitle(ProjectModuleKind::BoundaryCondition);
    boundaryModule.densityBased = isDensityBased();
    boundaryModule.targetNames = boundaryTargetNames();
    const QFileInfoList fieldFiles = zeroDir.entryInfoList(
        QDir::Files | QDir::NoDotAndDotDot, QDir::Name);
    for (const QFileInfo& fieldFile : fieldFiles) {
        if (fieldFile.fileName().startsWith(QStringLiteral("._"))) continue;
        ConfigDocument document;
        if (!document.load(fieldFile.absoluteFilePath(), error)) return false;
        const int documentIndex = documents_.size();
        documents_.push_back(document);
        for (int entryIndex = 0;
             entryIndex < document.entries().size(); ++entryIndex) {
            const ConfigEntry& entry = document.entries().at(entryIndex);
            if (entry.key.compare(QStringLiteral("internalField"),
                                  Qt::CaseInsensitive) == 0) {
                initialModule.entries.push_back({documentIndex, entryIndex});
            } else if (entry.section.startsWith(QStringLiteral("boundaryField/"),
                                                Qt::CaseInsensitive)) {
                boundaryModule.entries.push_back({documentIndex, entryIndex});
            }
        }
    }
    modules_.push_back(initialModule);
    modules_.push_back(boundaryModule);

    ProjectModule outputModule;
    outputModule.kind = ProjectModuleKind::Output;
    outputModule.title = moduleTitle(ProjectModuleKind::Output);
    modules_.push_back(outputModule);

    const struct {
        const char* path;
        ProjectModuleKind kind;
    } optionalModules[] = {
        {"constant/IBMProperties", ProjectModuleKind::IBM},
        {"constant/MRFProperties", ProjectModuleKind::MRF}
    };
    for (const auto& optional : optionalModules) {
        const QString path = QDir(directory_).filePath(QString::fromUtf8(optional.path));
        if (!QFileInfo(path).isFile()) continue;
        if (!loadDocument(path, optional.kind, moduleTitle(optional.kind),
                          true, error)) {
            return false;
        }
    }
    return true;
}

bool CaseProject::loadDocument(const QString& path,
                               ProjectModuleKind kind,
                               const QString& title,
                               bool optional,
                               QString* error)
{
    if (!QFileInfo::exists(path)) {
        const QString content = templateFor(kind);
        if (!writeFile(path, content, error)) return false;
    }

    ConfigDocument document;
    if (!document.load(path, error)) return false;
    const int documentIndex = documents_.size();
    documents_.push_back(document);

    ProjectModule module;
    module.kind = kind;
    module.title = title;
    module.documentIndex = documentIndex;
    module.optional = optional;
    for (int i = 0; i < document.entries().size(); ++i) {
        module.entries.push_back({documentIndex, i});
    }
    modules_.push_back(module);
    return true;
}

bool CaseProject::setValue(int documentIndex,
                           int entryIndex,
                           const QString& value,
                           QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].setValue(entryIndex, value, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::appendAssignments(int documentIndex,
                                    const QString& section,
                                    const QStringList& assignments,
                                    bool newSection,
                                    QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].appendAssignments(
            section, assignments, newSection, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::addBlockMeshPatch(int documentIndex,
                                    const QString& name,
                                    const QString& type,
                                    QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].addBlockMeshPatch(name, type, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::addBlockMeshPatchFace(int documentIndex,
                                        const QString& patchName,
                                        const QString& face,
                                        QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].addBlockMeshPatchFace(patchName, face, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::renameBlockMeshPatch(int documentIndex,
                                       const QString& oldName,
                                       const QString& newName,
                                       QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].renameBlockMeshPatch(oldName, newName, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::removeBlockMeshPatch(int documentIndex, const QString& name,
                                       QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].removeBlockMeshPatch(name, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::removeEntry(int documentIndex, int entryIndex, QString* error)
{
    if (documentIndex < 0 || documentIndex >= documents_.size()) {
        if (error) *error = QStringLiteral("无效的业务配置文档索引");
        return false;
    }
    if (!documents_[documentIndex].removeEntry(entryIndex, error)
        || !documents_[documentIndex].save(error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::saveAll(QString* error) const
{
    for (const ConfigDocument& document : documents_) {
        if (!document.save(error)) return false;
    }
    return true;
}

QString CaseProject::moduleFile(ProjectModuleKind kind) const
{
    for (const ProjectModule& module : modules_) {
        if (module.kind != kind) continue;
        const int documentIndex = module.documentIndex;
        if (documentIndex < 0 || documentIndex >= documents_.size()) {
            return {};
        }
        return documents_.at(documentIndex).path();
    }
    return {};
}

bool CaseProject::isSonicCase() const
{
    const QDir caseDir(directory_);
    return QFileInfo(caseDir.filePath(QStringLiteral("case.yaml"))).isFile()
        || (QFileInfo(caseDir.filePath(QStringLiteral("system/outputDict"))).isFile()
            && QFileInfo(caseDir.filePath(QStringLiteral("system/sonicDict"))).isFile());
}

bool CaseProject::isDensityBased() const
{
    const QString native = readFile(
        QDir(directory_).filePath(QStringLiteral("solver.yaml")));
    if (!native.isEmpty()) {
        return native.contains(
            QRegularExpression(QStringLiteral(R"(["']?type["']?\s*:\s*["']?densityBase["']?)"),
                               QRegularExpression::CaseInsensitiveOption));
    }
    const QString properties = readFile(
        QDir(directory_).filePath(QStringLiteral("system/solverProperties")));
    const QString type = foamValue(properties, QStringLiteral("type"));
    return type.compare(QStringLiteral("densityBase"),
                        Qt::CaseInsensitive) == 0;
}

QStringList CaseProject::boundaryTargetNames() const
{
    const QDir caseDir(directory_);
    const QStringList foamCandidates{
        caseDir.filePath(QStringLiteral("constant/polyMesh/boundary")),
        caseDir.filePath(QStringLiteral("constant/boundary"))};
    for (const QString& boundary : foamCandidates) {
        if (QFileInfo(boundary).isFile()) return foamBoundaryNames(boundary);
    }

    const QStringList candidates{
        caseDir.filePath(QStringLiteral("mesh/mesh_sets.sfm")),
        caseDir.filePath(QStringLiteral("mesh/mesh.sfm")),
        caseDir.filePath(QStringLiteral("constant/mesh_sets.sfm")),
        caseDir.filePath(QStringLiteral("constant/mesh.sfm"))};
    for (const QString& candidate : candidates) {
        if (!QFileInfo(candidate).isFile()) continue;
        const QStringList names = sfmSetNames(candidate);
        if (!names.isEmpty()) return names;
    }
    return {};
}

bool CaseProject::addField(ProjectModuleKind kind,
                           const QString& fieldName,
                           QString* error)
{
    if (kind != ProjectModuleKind::InitialCondition
        && kind != ProjectModuleKind::BoundaryCondition) {
        if (error) *error = QStringLiteral("只能从初始条件或边界条件新增字段");
        return false;
    }
    const QString name = fieldName.trimmed();
    if (!QRegularExpression(QStringLiteral(R"(^[A-Za-z_][A-Za-z0-9_.-]*$)"))
             .match(name).hasMatch()) {
        if (error) *error = QStringLiteral("字段名无效");
        return false;
    }
    const QString path = QDir(directory_).filePath(QStringLiteral("0/%1").arg(name));
    if (QFileInfo(path).isFile()) return reload(error);
    const bool vector = name == QStringLiteral("U");
    QString content = QStringLiteral(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
        "    class %1;\n    object %2;\n}\n\n"
        "dimensions [0 0 0 0 0 0 0];\n"
        "internalField uniform %3;\n"
        "boundaryField\n{\n")
        .arg(vector ? QStringLiteral("volVectorField")
                    : QStringLiteral("volScalarField"),
             name, vector ? QStringLiteral("(0 0 0)") : QStringLiteral("0"));
    for (const QString& patch : boundaryTargetNames()) {
        content += QStringLiteral(
            "    %1\n    {\n        type zeroGradient;\n    }\n").arg(patch);
    }
    content += QStringLiteral("}\n");
    return writeFile(path, content, error) && reload(error);
}

QString CaseProject::resultDirectory() const
{
    return QDir(directory_).filePath(QStringLiteral("result"));
}

QString CaseProject::outputTaskName() const
{
    const QDir dir(directory_);
    QString value = foamValue(readFile(dir.filePath(QStringLiteral("system/outputDict"))),
                              QStringLiteral("jobName"));
    value = value.trimmed();
    if (value.isEmpty()) value = projectName(directory_);
    return value;
}

QStringList CaseProject::outputTaskNames() const
{
    const QString content = readFile(
        QDir(directory_).filePath(QStringLiteral("system/outputDict")));
    QStringList names;
    const QString jobs = foamValue(content, QStringLiteral("jobs"));
    const QRegularExpression namePattern(QStringLiteral(R"([A-Za-z_][A-Za-z0-9_.-]*)"));
    for (auto match = namePattern.globalMatch(jobs); match.hasNext();) {
        const QString name = match.next().captured();
        if (!names.contains(name, Qt::CaseInsensitive)) names.push_back(name);
    }
    const QString active = outputTaskName();
    if (!names.contains(active, Qt::CaseInsensitive)) names.prepend(active);
    return names;
}

bool CaseProject::setOutputTaskName(const QString& taskName, QString* error)
{
    QString name = taskName.trimmed();
    name.replace(QRegularExpression(QStringLiteral("\\s+")),
                 QStringLiteral("_"));
    name.replace(QRegularExpression(QStringLiteral(R"([/\\:])")),
                 QStringLiteral("_"));
    while (name.startsWith(QLatin1Char('.'))) name.remove(0, 1);
    if (name.isEmpty()) {
        if (error) *error = QStringLiteral("任务名称不能为空");
        return false;
    }
    QStringList jobs = outputTaskNames();
    if (!jobs.contains(name, Qt::CaseInsensitive)) jobs.push_back(name);
    const QDir dir(directory_);
    const QString outputPath = dir.filePath(QStringLiteral("system/outputDict"));
    return replaceFoamValueInFile(outputPath, QStringLiteral("jobName"), name,
                                  error)
        && replaceFoamValueInFile(outputPath, QStringLiteral("outputDir"),
                                  QStringLiteral("result"), error)
        && replaceFoamValueInFile(outputPath, QStringLiteral("jobs"),
                                  QStringLiteral("(%1)").arg(
                                      jobs.join(QLatin1Char(' '))), error)
        && reload(error);
}

bool CaseProject::addOptionalModule(ProjectModuleKind kind, QString* error)
{
    QString solverFlag;
    QString path;
    switch (kind) {
        case ProjectModuleKind::IBM:
            solverFlag = QStringLiteral("IBM");
            path = QDir(directory_).filePath(QStringLiteral("constant/IBMProperties"));
            break;
        case ProjectModuleKind::MRF:
            solverFlag = QStringLiteral("MRF");
            path = QDir(directory_).filePath(QStringLiteral("constant/MRFProperties"));
            break;
        case ProjectModuleKind::Parallel:
            return updateSolverFlag(QStringLiteral("parallel.enabled"), true, error)
                && reload(error);
            break;
        case ProjectModuleKind::Multiphase:
            path = QDir(directory_).filePath(QStringLiteral("constant/phaseProperties"));
            break;
        case ProjectModuleKind::PhaseChange:
            path = QDir(directory_).filePath(QStringLiteral("constant/phaseChange"));
            break;
        default:
            if (error) *error = QStringLiteral("该模块不是可选模块");
            return false;
    }

    if (!QFileInfo::exists(path)
        && !writeFile(path, templateFor(kind), error)) {
        return false;
    }
    if (!solverFlag.isEmpty()
        && !updateSolverFlag(solverFlag, true, error)) {
        return false;
    }
    return reload(error);
}

bool CaseProject::setMeshFile(const QString& fileName, QString* error)
{
    if (!isSonicCase()) {
        if (error) {
            *error = QStringLiteral(
                "OpenFOAM case 不使用 sonicDict 的 meshFiles；请运行 blockMesh 或 snappyHexMesh。");
        }
        return false;
    }
    QString files = fileName.trimmed();
    if (files.isEmpty()
        || files == QStringLiteral("mesh.sfm")
        || files == QStringLiteral("constant/mesh.sfm")) {
        files = QStringLiteral("constant/mesh.sfm");
    }
    const QString value = QStringLiteral("(%1)").arg(files);
    return replaceFoamValueInFile(
               QDir(directory_).filePath(QStringLiteral("system/sonicDict")),
               QStringLiteral("meshFiles"), value, error)
        && reload(error);
}

bool CaseProject::importMesh(const QString& sourcePath, QString* error)
{
    const QDir caseDir(directory_);
    const QString meshTarget =
        caseDir.filePath(QStringLiteral("constant/mesh.sfm"));
    if (!copyFileToPath(sourcePath, meshTarget, error)) {
        return false;
    }

    QStringList meshEntries{QStringLiteral("constant/mesh.sfm")};
    return setMeshFile(meshEntries.join(QLatin1Char(' ')), error);
}

bool CaseProject::importStl(const QString& sourcePath, QString* error)
{
    if (!addOptionalModule(ProjectModuleKind::IBM, error)) return false;
    const QString targetName = QFileInfo(sourcePath).fileName();
    const QString targetPath = QDir(directory_).filePath(
        QStringLiteral("constant/triSurface/%1").arg(targetName));
    if (!copyFileToPath(sourcePath, targetPath, error)) {
        return false;
    }

    const QString ibmPath = QDir(directory_).filePath(
        QStringLiteral("constant/IBMProperties"));
    return replaceFoamBlockValueInFile(
               ibmPath, QStringLiteral("IBM"),
               QStringLiteral("geometryFiles"),
               QStringLiteral("(%1)").arg(targetName), error)
        && reload(error);
}

bool CaseProject::writeGeometry(const GeometrySketch& sketch, QString* error)
{
    if (sketch.points.size() != 4) {
        if (error) {
            *error = QStringLiteral(
                "当前几何生成器只接受四个不同顶点，并将其挤出为一个结构块。");
        }
        return false;
    }

    QString text = QStringLiteral(
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class dictionary;\n"
        "    object blockMeshDict;\n"
        "}\n\n"
        "convertToMeters 1;\n\n"
        "vertices\n"
        "(\n");
    for (int i = 0; i < sketch.points.size(); ++i) {
        const QPointF point = sketch.points.at(i);
        text += QStringLiteral("    (%1 %2 0)\n")
                    .arg(point.x(), 0, 'g', 12)
                    .arg(-point.y(), 0, 'g', 12);
    }
    for (int i = 0; i < sketch.points.size(); ++i) {
        const QPointF point = sketch.points.at(i);
        text += QStringLiteral("    (%1 %2 1)\n")
                    .arg(point.x(), 0, 'g', 12)
                    .arg(-point.y(), 0, 'g', 12);
    }
    text += QStringLiteral(
        ");\n\n"
        "blocks\n"
        "(\n"
        "    hex (0 1 2 3 4 5 6 7) (40 20 1) simpleGrading (1 1 1)\n"
        ");\n\n"
        "edges\n"
        "(\n");

    for (const GeometryArc& arc : sketch.arcs) {
        if (arc.start < 0 || arc.start >= 4
            || arc.end < 0 || arc.end >= 4) {
            continue;
        }
        const QPointF middle = arc.middle;
        text += QStringLiteral("    arc %1 %2 (%3 %4 0)\n")
                    .arg(arc.start)
                    .arg(arc.end)
                    .arg(middle.x(), 0, 'g', 12)
                    .arg(-middle.y(), 0, 'g', 12);
        text += QStringLiteral("    arc %1 %2 (%3 %4 1)\n")
                    .arg(arc.start + 4)
                    .arg(arc.end + 4)
                    .arg(middle.x(), 0, 'g', 12)
                    .arg(-middle.y(), 0, 'g', 12);
    }

    text += QStringLiteral(
        ");\n\n"
        "boundary\n"
        "(\n"
        "    left\n"
        "    {\n"
        "        type patch;\n"
        "        faces ((0 3 7 4));\n"
        "    }\n"
        "    right\n"
        "    {\n"
        "        type patch;\n"
        "        faces ((1 2 6 5));\n"
        "    }\n"
        "    walls\n"
        "    {\n"
        "        type wall;\n"
        "        faces ((0 1 5 4) (3 2 6 7));\n"
        "    }\n"
        "    frontBack\n"
        "    {\n"
        "        type empty;\n"
        "        faces ((0 1 2 3) (4 5 6 7));\n"
        "    }\n"
        ");\n\n"
        "mergePatchPairs\n"
        "(\n"
        ");\n");

    const QString geometryPath = QDir(directory_).filePath(
        defaultFileName(ProjectModuleKind::Geometry));
    if (!writeFile(geometryPath, text, error)) return false;
    return reload(error);
}

bool CaseProject::updateSolverFlag(const QString& key,
                                   bool enabled,
                                   QString* error)
{
    QString block = QStringLiteral("solver");
    QString entryKey = key;
    const int dot = key.indexOf(QLatin1Char('.'));
    if (dot > 0) {
        block = key.left(dot);
        entryKey = key.mid(dot + 1);
    }

    const QString dictPath =
        QDir(directory_).filePath(QStringLiteral("system/sonicDict"));
    QFile file(dictPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = file.errorString();
        return false;
    }
    QStringList lines = QString::fromUtf8(file.readAll())
                            .split(QLatin1Char('\n'), Qt::KeepEmptyParts);
    if (!lines.isEmpty() && lines.last().isEmpty()) lines.removeLast();
    const QRegularExpression keyPattern(
        QStringLiteral(R"(^(\s*)%1\s+.*;\s*(?://.*)?$)")
            .arg(QRegularExpression::escape(entryKey)),
        QRegularExpression::CaseInsensitiveOption);
    int insertion = -1;
    bool waitingBrace = false;
    bool inBlock = false;
    for (int i = 0; i < lines.size(); ++i) {
        const QString trimmed = lines.at(i).trimmed();
        if (!inBlock) {
            if (waitingBrace && trimmed == QStringLiteral("{")) {
                inBlock = true;
                waitingBrace = false;
                insertion = i + 1;
                continue;
            }
            if (trimmed.compare(block, Qt::CaseInsensitive) == 0) {
                waitingBrace = true;
                continue;
            }
            if (trimmed.startsWith(block + QStringLiteral("{"),
                                   Qt::CaseInsensitive)) {
                inBlock = true;
                insertion = i + 1;
                continue;
            }
            continue;
        }
        if (trimmed.startsWith(QLatin1Char('}'))) break;
        const QRegularExpressionMatch keyMatch = keyPattern.match(lines.at(i));
        if (keyMatch.hasMatch()) {
            lines[i] = keyMatch.captured(1) + entryKey + QLatin1Char(' ')
                     + (enabled ? QStringLiteral("true")
                                : QStringLiteral("false"))
                     + QLatin1Char(';');
            return writeFile(dictPath, lines.join(QLatin1Char('\n')), error);
        }
        if (!trimmed.isEmpty()) insertion = i + 1;
    }

    if (insertion < 0) {
        if (!lines.isEmpty() && !lines.last().trimmed().isEmpty()) {
            lines.push_back(QString());
        }
        lines.push_back(block);
        lines.push_back(QStringLiteral("{"));
        lines.push_back(QStringLiteral("    %1 %2;")
                            .arg(entryKey,
                                 enabled ? QStringLiteral("true")
                                         : QStringLiteral("false")));
        lines.push_back(QStringLiteral("}"));
        return writeFile(dictPath, lines.join(QLatin1Char('\n')), error);
    }

    lines.insert(insertion,
                 QStringLiteral("    %1 %2;")
                     .arg(entryKey, enabled ? QStringLiteral("true")
                                            : QStringLiteral("false")));
    return writeFile(dictPath, lines.join(QLatin1Char('\n')), error);
}

QString CaseProject::moduleTitle(ProjectModuleKind kind)
{
    switch (kind) {
        case ProjectModuleKind::Geometry: return QStringLiteral("几何");
        case ProjectModuleKind::Mesh: return QStringLiteral("网格");
        case ProjectModuleKind::Solver: return QStringLiteral("求解器");
        case ProjectModuleKind::Turbulence: return QStringLiteral("湍流模型");
        case ProjectModuleKind::Multiphase: return QStringLiteral("多相流");
        case ProjectModuleKind::PhaseChange: return QStringLiteral("相变");
        case ProjectModuleKind::InitialCondition: return QStringLiteral("初始条件");
        case ProjectModuleKind::BoundaryCondition: return QStringLiteral("边界条件");
        case ProjectModuleKind::Output: return QStringLiteral("输出");
        case ProjectModuleKind::IBM: return QStringLiteral("IBM");
        case ProjectModuleKind::MRF: return QStringLiteral("MRF");
        case ProjectModuleKind::Parallel: return QStringLiteral("Parallel");
    }
    return {};
}

QString CaseProject::defaultFileName(ProjectModuleKind kind)
{
    switch (kind) {
        case ProjectModuleKind::Geometry: return QStringLiteral("system/blockMeshDict");
        case ProjectModuleKind::Solver: return QStringLiteral("system/controlDict");
        case ProjectModuleKind::Turbulence:
            return QStringLiteral("constant/turbulenceProperties");
        case ProjectModuleKind::Multiphase:
            return QStringLiteral("constant/phaseProperties");
        case ProjectModuleKind::PhaseChange:
            return QStringLiteral("constant/phaseChange");
        case ProjectModuleKind::InitialCondition: return QStringLiteral("0/rho");
        case ProjectModuleKind::BoundaryCondition: return QStringLiteral("0/U");
        case ProjectModuleKind::IBM: return QStringLiteral("constant/IBMProperties");
        case ProjectModuleKind::MRF: return QStringLiteral("constant/MRFProperties");
        case ProjectModuleKind::Parallel: return QStringLiteral("system/sonicDict");
        default: return {};
    }
}

bool CaseProject::writeFile(const QString& path,
                            const QString& content,
                            QString* error)
{
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text
                   | QIODevice::Truncate)) {
        if (error) *error = file.errorString();
        return false;
    }
    const QByteArray data = content.toUtf8();
    if (file.write(data) != data.size()) {
        if (error) *error = file.errorString();
        return false;
    }
    return true;
}

} // namespace SF::GUI
