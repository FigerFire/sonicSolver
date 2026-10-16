#include "app/application/SF_inspection.h"
#include "SF_project.h"
#include "app/application/model/SF_caseTemplate.h"
#include "app/application/model/SF_model.h"
#include "infrastructure/io/case/SF_case.h"
#include "infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QSaveFile>
#include <QTemporaryDir>
#include <filesystem>
#include <map>
namespace SF::GUI {
namespace {
bool reject(QString* error,const std::exception& e){if(error)*error=QString::fromUtf8(e.what());return false;}
void put(const QString& path,const QByteArray& data){QDir().mkpath(QFileInfo(path).absolutePath());QSaveFile f(path);if(!f.open(QIODevice::WriteOnly)||f.write(data)!=data.size()||!f.commit())throw std::runtime_error(("Cannot save "+path).toStdString());}
Model::Description description(const QString& dir){return CaseIO::read(dir.toStdString());}
void commit(const Model::Description& model,const QString& directory) {
    QTemporaryDir stage(QDir(directory).filePath(".native-composition-XXXXXX"));
    if(!stage.isValid())throw std::runtime_error("Cannot stage native composition");
    CaseIO::write(model,stage.path().toStdString());auto candidate=CaseIO::read(stage.path().toStdString());
    candidate.caseDir=directory.toStdString();(void)Application::inspectCase(Application::ModelLoader::build(candidate));
    std::map<QString,QByteArray> previous;
    std::vector<QString> changed;
    try {
        for(const auto& document:candidate.comments) {
            const QString relative=QString::fromStdString(document.first);
            const QString target=QDir(directory).filePath(relative);
            QFile old(target);if(old.exists()) {
                if(!old.open(QIODevice::ReadOnly))throw std::runtime_error("Cannot read document before replacement");
                previous[target]=old.readAll();
            }
            QFile source(QDir(stage.path()).filePath(relative));
            if(!source.open(QIODevice::ReadOnly))throw std::runtime_error("Cannot read validated staged document");
            put(target,source.readAll());changed.push_back(target);
        }
    } catch(...) {
        for(const auto& target:changed)if(previous.count(target))put(target,previous.at(target));
        else QFile::remove(target);
        throw;
    }
}

ProjectModuleKind kind(const ConfigEntry& header) {
    if(header.section=="mesh")return header.sectionLabel=="generation" ? ProjectModuleKind::Geometry : ProjectModuleKind::Mesh;
    if(header.section=="fields") {
        if(header.sectionLabel=="internalField")return ProjectModuleKind::InitialCondition;
        if(header.sectionLabel=="boundaries")return ProjectModuleKind::BoundaryCondition;
    }
    return ProjectModuleKind::Solver;
}
}
bool CaseProject::open(const QString& directory,QString* error) {
    try {
        const QString absolute=QFileInfo(directory).absoluteFilePath();
        if(!QFileInfo(QDir(absolute).filePath("case.yaml")).exists()) {
            Application::createNativeCase(absolute.toStdString());
        }
        const auto old=directory_;directory_=absolute;casePath_=absolute;
        if(!reload(error)){directory_=old;return false;}return true;
    }catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::reload(QString* error) {
    try {
        const auto m=description(directory_);(void)Application::ModelLoader::build(m);
        QVector<ConfigDocument> documents;QVector<ProjectModule> modules;
        for(const auto& entry:m.comments) {
            ConfigDocument doc;const QString path=QDir(directory_).filePath(QString::fromStdString(entry.first));
            if(!doc.load(path,error))return false;
            ProjectModule module;module.kind=kind(doc.entries().front());module.title=QString::fromStdString(entry.first);
            module.documentIndex=documents.size();module.entries={{static_cast<int>(documents.size()),0}};
            for(const auto& object:m.objects)if(object.source==path.toStdString()) {
                module.optional=object.category=="models";
                if(object.type=="turbulence")module.kind=ProjectModuleKind::Turbulence;
                else if(object.type=="IBM")module.kind=ProjectModuleKind::IBM;
                else if(object.type=="multiPhase")module.kind=ProjectModuleKind::Multiphase;
                const auto& params=object.parameters;
                const bool disabled=params.value("enabled",true)==false || (params.contains(object.name) && params.at(object.name).is_object() && !params.at(object.name).value("enabled",true));
                module.enabled=!disabled;
                if(disabled)module.title+=" [configured, inactive]";
            }
            documents.push_back(std::move(doc));modules.push_back(std::move(module));
        }
        ProjectModule output;output.kind=ProjectModuleKind::Output;output.title="Results";modules.push_back(output);
        documents_=std::move(documents);modules_=std::move(modules);return true;
    }catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::setValue(int di,int ei,const QString& value,QString* error) {
    if(di<0||di>=documents_.size()){if(error)*error="Invalid native document";return false;}
    const auto old=documents_[di];if(!documents_[di].setValue(ei,value,error))return false;
    if(!saveAll(error)){documents_[di]=old;return false;}return reload(error);
}
bool CaseProject::saveAll(QString* error) const {
    std::map<QString,QByteArray> originals;std::vector<QString> committed;
    try {
        QTemporaryDir stage(QDir(directory_).filePath(".native-edit-XXXXXX"));
        if(!stage.isValid())throw std::runtime_error("Cannot stage native documents");
        for(const auto& doc:documents_) {
            QFile input(doc.path());if(!input.open(QIODevice::ReadOnly))throw std::runtime_error("Cannot read existing document");
            originals[doc.path()]=input.readAll();
            const QString relative=QDir(directory_).relativeFilePath(doc.path());
            put(QDir(stage.path()).filePath(relative),doc.entries().front().value.toUtf8());
        }
        // The production reader performs identity, duplicate-key, dependency and canonical-path validation.
        const auto validated=description(stage.path());
        // Preserve the actual case root for mesh/geometry references in typed configuration validation.
        auto input=validated;input.caseDir=directory_.toStdString();(void)Application::inspectCase(Application::ModelLoader::build(input));
        for(const auto& doc:documents_)if(doc.dirty()) {put(doc.path(),doc.entries().front().value.toUtf8());committed.push_back(doc.path());}
        return true;
    }catch(const std::exception& e){for(const auto& path:committed)put(path,originals.at(path));return reject(error,e);}
}
bool CaseProject::addOptionalModule(ProjectModuleKind selected,QString* error) {
    try {
        auto m=description(directory_);std::string name;
        if(selected==ProjectModuleKind::Turbulence)name="turbulence";
        else if(selected==ProjectModuleKind::IBM)name="IBM";
        else if(selected==ProjectModuleKind::MRF)name="MRF";
        else throw std::runtime_error("Unsupported module template; edit the native registry/document explicitly");
        for(const auto& object:m.objects)if(object.name==name)throw std::runtime_error("Module already registered");
        Model::ObjectDescriptor o;o.name=o.type=name;o.category="models";
        o.parameters={{"enabled",false}};
        if(name=="turbulence")o.parameters["simulationType"]="laminar";
        if(name=="IBM")o.parameters={{"IBM",{{"enabled",false},{"type","ghost"}}}};
        m.objects.push_back(o);(void)Application::ModelLoader::build(m);commit(m,directory_);return reload(error);
    }catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::setMeshFile(const QString& fileName,QString* error) {
    try {auto m=description(directory_);m.mesh["files"]={fileName.toStdString()};(void)Application::ModelLoader::build(m);commit(m,directory_);return reload(error);}catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::importMesh(const QString& source,QString* error,bool replaceGenerated) {
    const auto target=QDir(directory_).filePath("mesh/"+QFileInfo(source).fileName());
    QByteArray previous;bool replaced=false,existed=QFileInfo(target).exists();
    try {
        const auto inspection=Application::inspectCase(directory_.toStdString());
        MultiBlockMesh mesh;
        if(!mesh.loadFiles(QFileInfo(source).absolutePath().toStdString(),
            {QFileInfo(source).fileName().toStdString()},inspection.mesh))
            throw std::runtime_error("Invalid/unsupported SFM payload: production mesh reload failed");
        if(source!=target) {
            if(existed && !replaceGenerated)throw std::runtime_error("Existing SFM is not overwritten by import");
            if(existed){QFile old(target);if(!old.open(QIODevice::ReadOnly))throw std::runtime_error("Cannot preserve previous generated mesh");previous=old.readAll();}
            QFile input(source);if(!input.open(QIODevice::ReadOnly))throw std::runtime_error("Cannot read SFM payload");
            put(target,input.readAll());replaced=true;
        }
        if(!setMeshFile(QDir(directory_).relativeFilePath(target),error))throw std::runtime_error(error?error->toStdString():"Cannot register SFM");
        return true;
    }catch(const std::exception& e){if(replaced){if(existed)put(target,previous);else QFile::remove(target);}return reject(error,e);}
}
bool CaseProject::importStl(const QString& source,QString* error) {
    try {auto m=description(directory_);const auto name=QFileInfo(source).completeBaseName().toStdString();for(const auto& o:m.objects)if(o.name==name)throw std::runtime_error("Geometry already registered");
        QDir().mkpath(QDir(directory_).filePath("mesh/geometry"));const auto target=QDir(directory_).filePath("mesh/geometry/"+QFileInfo(source).fileName());if(!QFile::copy(source,target))throw std::runtime_error("Cannot copy STL");
        Model::ObjectDescriptor o;o.name=name;o.category="models";o.type="geometry";o.parameters={{"file",QDir(directory_).relativeFilePath(target).toStdString()}};m.objects.push_back(o);commit(m,directory_);return reload(error);
    }catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::writeGeometry(const GeometrySketch& sketch,QString* error) {
    try {
        if(sketch.points.size()!=4)throw std::runtime_error("Unsupported sketch topology: native extrusion requires four ordered corners; edit structuredMesh for multiblock geometry");
        auto m=description(directory_);using P=Model::Parameters;P vertices=P::array();for(double z:{0.,.1})for(const auto& p:sketch.points)vertices.push_back({p.x(),p.y(),z});
        P edges=P::array();for(const auto& arc:sketch.arcs)for(int offset:{0,4})edges.push_back({{"type","arc"},{"vertices",{arc.start+offset,arc.end+offset}},{"point",{arc.middle.x(),arc.middle.y(),offset?.1:0.}}});
        auto g=m.meshGeneration;
        if(g.empty())throw std::runtime_error("Register structuredMesh generation before editing geometry");
        g["vertices"]=vertices;g["edges"]=edges;
        m.meshGeneration=g;commit(m,directory_);return reload(error);
    }catch(const std::exception& e){return reject(error,e);}
}
bool CaseProject::setOutputTaskName(const QString& task,QString* error) {
    try {if(task.isEmpty()||task.contains('/')||task.contains('\\')||task=="..")throw std::runtime_error("Invalid output task name");auto m=description(directory_);m.runtime["output"]["jobName"]=task.toStdString();commit(m,directory_);return reload(error);}catch(const std::exception& e){return reject(error,e);}
}
QString CaseProject::moduleFile(ProjectModuleKind requested) const {for(const auto& m:modules_)if(m.kind==requested && m.documentIndex>=0)return documents_[m.documentIndex].path();return {};}
QString CaseProject::resultDirectory() const {const auto m=description(directory_);return QDir(directory_).filePath(QString::fromStdString(m.runtime.value("output",Model::Parameters::object()).value("outputDir",std::string("result"))));}
QString CaseProject::outputTaskName() const {const auto m=description(directory_);return QString::fromStdString(m.runtime.value("output",Model::Parameters::object()).value("jobName",m.name));}
QStringList CaseProject::outputTaskNames() const{return {outputTaskName()};}
}
