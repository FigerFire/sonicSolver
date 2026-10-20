#pragma once
#include <QString>
#include <QVector>
namespace SF::GUI {
/// @brief 原生文档文本项；复杂 WHAT/HOW 由同一 YAML 文档编辑，不建立 GUI 数学副本。
struct ConfigEntry {
    QString section,sectionLabel,key,value;
};
class ConfigDocument {
public:
    bool load(const QString& path,QString* error=nullptr);
    bool setValue(int entry,const QString& text,QString* error=nullptr);
    const QString& path() const {return path_;}
    QString fileName() const;
    const QVector<ConfigEntry>& entries() const {return entries_;}
    bool dirty() const {return dirty_;}
private:
    QString path_;QVector<ConfigEntry> entries_;bool dirty_=false;
};
}
