/// @file SF_project.h
/// @brief GUI 工程与中立配置文档的读写实现。

#pragma once

#include "SF_configDocument.h"

#include <QMetaType>
#include <QPointF>
#include <QString>
#include <QVector>

namespace SF::GUI {

/// @brief GUI 左侧业务模块类型。
enum class ProjectModuleKind {
    Geometry,
    Mesh,
    Solver,
    Turbulence,
    Multiphase,
    PhaseChange,
    InitialCondition,
    BoundaryCondition,
    Output,
    IBM,
    MRF,
    Parallel
};

/// @brief 指向某个配置文档条目的稳定引用。
struct ProjectEntryRef {
    int documentIndex = -1;
    int entryIndex = -1;
};

/// @brief 左侧树中的一个业务模块。
struct ProjectModule {
    ProjectModuleKind kind = ProjectModuleKind::Geometry;
    QString title;
    int documentIndex = -1;
    QVector<ProjectEntryRef> entries;
    /// @brief 边界条件模块可选择的 polyMesh patch 或 SFM set 名。
    QStringList targetNames;
    bool densityBased = false;
    bool optional = false;
};

/// @brief 输出任务列表中显示的一条 PVD 任务。
struct OutputTaskInfo {
    QString name;
    QString pvdPath;
    QString latestTime;
    double latestTimeValue = 0.0;
    double progress = 0.0;
    int frameCount = 0;
    /// @brief 该任务的 PVD 当前正在右侧 VTK 窗口显示。
    bool displayed = false;
    /// @brief 该任务是下一次求解写入的 jobName。
    bool active = false;
};

/// @brief 几何草图中的圆弧。
struct GeometryArc {
    int start = -1;
    int end = -1;
    QPointF middle;
};

/// @brief 可转换成单块挤出网格参数的二维几何草图。
struct GeometrySketch {
    QVector<QPointF> points;
    QVector<GeometryArc> arcs;
};

/// @brief 同时支持 sonicSolver 与 OpenFOAM 的 case 工作目录。
class CaseProject {
public:
    bool open(const QString& directory, QString* error = nullptr);
    bool setValue(int documentIndex,
                  int entryIndex,
                  const QString& value,
                  QString* error = nullptr);
    /// @brief 向某个配置文档 section 追加配置行。
    bool appendAssignments(int documentIndex,
                           const QString& section,
                           const QStringList& assignments,
                           bool newSection,
                           QString* error = nullptr);
    bool addBlockMeshPatch(int documentIndex,
                           const QString& name,
                           const QString& type,
                           QString* error = nullptr);
    bool addBlockMeshPatchFace(int documentIndex,
                               const QString& patchName,
                               const QString& face,
                               QString* error = nullptr);
    bool renameBlockMeshPatch(int documentIndex,
                              const QString& oldName,
                              const QString& newName,
                              QString* error = nullptr);
    bool removeBlockMeshPatch(int documentIndex, const QString& name,
                              QString* error = nullptr);
    bool removeEntry(int documentIndex, int entryIndex,
                     QString* error = nullptr);
    bool saveAll(QString* error = nullptr) const;
    bool addOptionalModule(ProjectModuleKind kind,
                           QString* error = nullptr);
    /// @brief 在 0/ 目录新增一个 OpenFOAM 字段文件。
    bool addField(ProjectModuleKind kind, const QString& fieldName,
                  QString* error = nullptr);
    /// @brief 设置 system/sonicDict 中的 meshFiles 项，通常在 GUI 生成 mesh.sfm 后调用。
    bool setMeshFile(const QString& fileName, QString* error = nullptr);
    bool importMesh(const QString& sourcePath, QString* error = nullptr);
    bool importStl(const QString& sourcePath, QString* error = nullptr);
    bool writeGeometry(const GeometrySketch& sketch,
                       QString* error = nullptr);
    /// @brief 设置当前输出任务名；求解结果将写为 result/<name>.pvd。
    bool setOutputTaskName(const QString& taskName, QString* error = nullptr);

    const QString& directory() const { return directory_; }
    const QString& casePath() const { return casePath_; }
    const QVector<ConfigDocument>& documents() const { return documents_; }
    const QVector<ProjectModule>& modules() const { return modules_; }
    /// @brief 存在 case.yaml 或旧 sonic 字典对时使用 sonicSolver/SFM 工作流。
    bool isSonicCase() const;
    bool isDensityBased() const;

    /// @brief 返回业务模块对应的配置文件路径。
    /// @param kind 业务模块类型。
    /// @return 配置文件绝对路径；模块不存在时返回空字符串。
    QString moduleFile(ProjectModuleKind kind) const;

    /// @brief 返回 case 配置的结果输出目录。
    /// @return 结果目录绝对路径。
    QString resultDirectory() const;
    /// @brief 返回当前任务名，也就是 controlDict/sonicDict 的 jobName。
    QString outputTaskName() const;
    /// @brief 返回 outputDict 中持久化的全部 job 名。
    QStringList outputTaskNames() const;

private:
    bool createDefaultProject(QString* error);
    bool reload(QString* error);
    bool loadDocument(const QString& path,
                      ProjectModuleKind kind,
                      const QString& title,
                      bool optional,
                      QString* error);
    bool updateSolverFlag(const QString& key,
                          bool enabled,
                          QString* error);
    static QString moduleTitle(ProjectModuleKind kind);
    static QString defaultFileName(ProjectModuleKind kind);
    static bool writeFile(const QString& path,
                          const QString& content,
                          QString* error);
    QStringList boundaryTargetNames() const;
    QString directory_;
    QString casePath_;
    QVector<ConfigDocument> documents_;
    QVector<ProjectModule> modules_;
};

} // namespace SF::GUI

Q_DECLARE_METATYPE(SF::GUI::GeometrySketch)
