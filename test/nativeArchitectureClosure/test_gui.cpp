#include "app/gui/IO/SF_project.h"
#include "infrastructure/io/case/SF_case.h"
#include "app/application/model/SF_model.h"
#include <QCoreApplication>
#include <QFile>
#include <QDir>
#include <fstream>
#include <filesystem>
#include <iostream>
using namespace SF;
void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}
int main(int argc,char** argv) {
    QCoreApplication app(argc,argv);namespace fs=std::filesystem;
    const auto work=fs::weakly_canonical(fs::current_path());
    const auto test=fs::weakly_canonical(fs::path(SF_TEST_SOURCE_DIR)/"test");
    require(work==test || work.string().rfind(test.string()+"/",0)==0,"GUI tests require a repository test/ working directory");
    const auto root=work/"guiContract";require(!fs::exists(root),"Use a fresh test directory");
    GUI::CaseProject p;QString error;require(p.open(QString::fromStdString(root.string()),&error),qPrintable(error));
    require(fs::exists(root/"case.yaml")&&!fs::exists(root/"system"),"GUI did not create a native-only project");
    auto m=CaseIO::read(root.string());require(!m.meshGeneration.empty(),"Generator was not registered");
    int runtime=-1;for(int i=0;i<p.documents().size();++i)if(p.documents()[i].fileName()=="runtime.yaml")runtime=i;
    require(runtime>=0,"Registered runtime missing from GUI navigation");
    auto text=p.documents()[runtime].entries()[0].value;const auto original=text;
    text.replace("type: runtime","type: invalidIdentity");require(!p.setValue(runtime,0,text,&error),"Invalid identity was saved");
    QFile file(QString::fromStdString((root/"solvers/runtime.yaml").string()));require(file.open(QIODevice::ReadOnly),"Missing original document");require(QString::fromUtf8(file.readAll())==original,"Rejected edit corrupted original");
    require(p.setValue(runtime,0,original+"\n# preserved native comment\n",&error),qPrintable(error));
    require(p.addOptionalModule(GUI::ProjectModuleKind::Turbulence,&error),qPrintable(error));
    require(p.importMesh(QString::fromStdString((root/"mesh/mesh.sfm").string()),&error,true),qPrintable(error));
    auto model=CaseIO::read(root.string());require(model.objects.back().type=="turbulence","Module registration missing");
    require(!model.objects.back().parameters.value("enabled",true),"New optional module was silently enabled");
    const auto rewritten=root.parent_path()/"guiRoundTrip";CaseIO::write(model,rewritten.string());
    const auto b=CaseIO::read(rewritten.string());require(b.meshGeneration==model.meshGeneration,"Mesh generation lost in round trip");
    QFile unregistered(QString::fromStdString((root/"models/ignored.yaml").string()));require(unregistered.open(QIODevice::WriteOnly),"Cannot create ignored test input");unregistered.write("bad unregistered document\n");unregistered.close();require(p.open(QString::fromStdString(root.string()),&error),"Unregistered input was loaded");
    const auto old=root.parent_path()/"oldOnly";fs::create_directories(old/"system");std::ofstream(old/"system/controlDict")<<"application sonicSolver;";GUI::CaseProject q;require(!q.open(QString::fromStdString(old.string()),&error),"Old-only case accepted");require(!fs::exists(old/"case.yaml"),"Old-only case mutated");
    std::cout<<"Native GUI creation, transactional save, optional registration and round-trip passed\n";
}
