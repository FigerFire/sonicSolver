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
    /// @brief 边界条件模块可选择的 SFM set 名。
    QStringList targetNames;
    bool optional = false;
    bool enabled = true;
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

/// @brief 使用严格 CaseIO 的原生 case 工作目录。
class CaseProject {
public:
    bool open(const QString& directory, QString* error = nullptr);
    bool setValue(int documentIndex,
                  int entryIndex,
                  const QString& value,
                  QString* error = nullptr);
    bool saveAll(QString* error = nullptr) const;
    bool addOptionalModule(ProjectModuleKind kind,
                           QString* error = nullptr);
    /// @brief 显式更新 mesh registry，不扫描未注册输入。
    bool setMeshFile(const QString& fileName, QString* error = nullptr);
    bool importMesh(const QString& sourcePath, QString* error = nullptr, bool replaceGenerated = false);
    bool importStl(const QString& sourcePath, QString* error = nullptr);
    bool writeGeometry(const GeometrySketch& sketch,
                       QString* error = nullptr);
    /// @brief 设置当前输出任务名；求解结果将写为 result/<name>.pvd。
    bool setOutputTaskName(const QString& taskName, QString* error = nullptr);

    const QString& directory() const { return directory_; }
    const QString& casePath() const { return casePath_; }
    const QVector<ConfigDocument>& documents() const { return documents_; }
    const QVector<ProjectModule>& modules() const { return modules_; }
    /// @brief 返回业务模块对应的配置文件路径。
    /// @param kind 业务模块类型。
    /// @return 配置文件绝对路径；模块不存在时返回空字符串。
    QString moduleFile(ProjectModuleKind kind) const;

    /// @brief 返回 case 配置的结果输出目录。
    /// @return 结果目录绝对路径。
    QString resultDirectory() const;
    /// @brief 返回当前任务名，也就是 native runtime 的 jobName。
    QString outputTaskName() const;
    /// @brief 返回 native runtime 的当前 job 名（历史结果由结果视图读取）。
    QStringList outputTaskNames() const;

private:
    bool reload(QString* error);
    QString directory_,casePath_;
    QVector<ConfigDocument> documents_;
    QVector<ProjectModule> modules_;
};

} // namespace SF::GUI
Q_DECLARE_METATYPE(SF::GUI::GeometrySketch)
