#include "SF_configDocument.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <QFile>
#include <QFileInfo>
#include <yaml-cpp/yaml.h>
namespace SF::GUI {
bool ConfigDocument::load(const QString& path,QString* error) {
    try {
        const auto doc=Serialization::readSonicFile(path.toStdString());
        QFile f(path);if(!f.open(QIODevice::ReadOnly))throw std::runtime_error(f.errorString().toStdString());
        path_=path;entries_={{QString::fromStdString(doc.object),QString::fromStdString(doc.type),fileName(),QString::fromUtf8(f.readAll())}};dirty_=false;return true;
    }catch(const std::exception& e){if(error)*error=QString::fromUtf8(e.what());return false;}
}
bool ConfigDocument::setValue(int index,const QString& text,QString* error) {
    try {
        if(index!=0 || entries_.size()!=1)throw std::runtime_error("Invalid native document entry");
        const auto parsed=YAML::Load(text.toStdString());
        if(!parsed.IsMap() || !parsed["SonicFile"].IsMap())throw std::runtime_error("Native document requires SonicFile header");
        entries_[0].value=text;dirty_=true;return true;
    }catch(const std::exception& e){if(error)*error=QString::fromUtf8(e.what());return false;}
}
QString ConfigDocument::fileName() const{return QFileInfo(path_).fileName();}
}
