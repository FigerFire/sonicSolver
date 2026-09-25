/// @file SF_configDocument.h
/// @brief GUI 工程与中立配置文档的读写实现。

#pragma once

#include <QString>
#include <QStringList>
#include <QVector>

namespace SF::GUI {

/// @brief GUI 属性树中的一个可编辑配置项。
struct ConfigEntry {
    QString section;
    QString sectionLabel;
    QString key;
    QString value;
    int startLine = -1;
    int endLine = -1;
    /// @brief 该条目使用 OpenFOAM 的 `key value;` 语法。
    bool foamAssignment = false;
    /// @brief blockMesh 顶点坐标；树中以 `x y z` 显示，写回时恢复圆括号。
    bool coordinateTuple = false;
    bool caseSection = false;

    bool isBoolean() const;
};

/// @brief GUI 使用的轻量配置文档模型，支持 `.sf` 与 OpenFOAM 字典。
class ConfigDocument {
public:
    bool load(const QString& path, QString* error = nullptr);
    bool save(QString* error = nullptr) const;
    bool setValue(int entryIndex, const QString& value, QString* error = nullptr);
    /// @brief 向已有 section 追加若干配置行。
    /// @param section section 名称，不包含方括号。
    /// @param assignments 完整赋值行，例如 `v0 = [0,0,0]`。
    /// @param newSection 是否创建一个同名的新 section。
    /// @param error 失败诊断。
    /// @return 成功写入内存文档返回 true。
    bool appendAssignments(const QString& section,
                           const QStringList& assignments,
                           bool newSection,
                           QString* error = nullptr);
    /// @brief 在 blockMeshDict 的 boundary 列表中新增一个 patch。
    bool addBlockMeshPatch(const QString& name,
                           const QString& type,
                           QString* error = nullptr);
    /// @brief 向指定 blockMesh patch 的 faces 列表新增一个四边形面。
    bool addBlockMeshPatchFace(const QString& patchName,
                               const QString& face,
                               QString* error = nullptr);
    /// @brief 重命名 blockMeshDict 中的 patch。
    bool renameBlockMeshPatch(const QString& oldName,
                              const QString& newName,
                              QString* error = nullptr);
    /// @brief 从 blockMeshDict 的 boundary 列表删除一个完整 patch（含全部 faces）。
    bool removeBlockMeshPatch(const QString& name, QString* error = nullptr);
    /// @brief 删除一个可编辑的 blockMesh 行条目，例如 patch face。
    bool removeEntry(int entryIndex, QString* error = nullptr);

    const QString& path() const { return path_; }
    QString fileName() const;
    const QVector<ConfigEntry>& entries() const { return entries_; }
    bool dirty() const { return dirty_; }

private:
    void parse();
    void parseFoamDictionary();
    void parseBlockMeshDict();
    static int squareBracketBalance(const QString& text);
    static QString valueWithoutComment(const QString& text);
    static QString inlineComment(const QString& text);

    QString path_;
    QStringList lines_;
    QVector<ConfigEntry> entries_;
    bool dirty_ = false;
};

} // namespace SF::GUI
