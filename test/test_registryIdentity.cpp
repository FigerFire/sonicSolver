/// @brief 注册身份、内部角色与 object 目录的双重核验；所有证据位于 test/。
#include "infrastructure/io/case/SF_case.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs=std::filesystem;
using namespace SF;
using P=Model::Parameters;

void require(bool valid,const std::string& message) {
    if(!valid)throw std::runtime_error(message);
}
int main(int argc,char** argv) {try {
    if(argc==2) {
        const auto read=CaseIO::read(argv[1]);
        std::cout<<"Read registered case "<<argv[1]<<": "<<read.objects.size()<<" objects, "<<read.fields.size()<<" fields\n";
        return 0;
    }
    const auto source=fs::path(SF_TEST_SOURCE_DIR);
    const auto seed=CaseIO::read((source/"test/Sod/sodCase").string());
    const auto directory=source/"test/registryIdentity/cases/identity";
    auto reset=[&] {fs::remove_all(directory);CaseIO::write(seed,directory.string());};
    auto doc=[&](const std::string& path){return Serialization::readSonicFile((directory/path).string());};
    auto save=[&](const std::string& path,const Serialization::SonicDocument& value){Serialization::writeSonicFile(value,(directory/path).string());};
    auto reject=[&](const std::string& label,const std::string& token) {
        try{(void)CaseIO::read(directory.string());}
        catch(const std::exception& e){require(std::string(e.what()).find(token)!=std::string::npos,label+": "+e.what());std::cout<<"rejected "<<label<<": "<<e.what()<<'\n';return;}
        throw std::runtime_error("Accepted invalid "+label);
    };
    reset();
    require(doc("solvers/state.yaml").type=="state","STATE title is not registry identity");
    require(doc("solvers/execution.yaml").type=="execution","HOW title is not registry identity");
    require(doc("solvers/solvers.yaml").body.at("equations").at("type")=="equationSystem","Default flow identity was not renamed");
    require(doc("fields/boundaries.yaml").type=="boundaries","Boundary title is not registry identity");
    const auto loaded=CaseIO::read(directory.string());
    require(std::any_of(loaded.objects.begin(),loaded.objects.end(),[](const auto& o){return o.name=="state" && o.category=="solver" && o.type=="stateRegistry";}),"Registry role lost during typed dispatch");
    require(loaded.fields.size()==seed.fields.size(),"Registered field metadata lost");

    // Unregistered files may even be malformed: their existence is not a loading request.
    for(const auto* relative:{"solvers/unregistered.yaml","state/state.yaml","equations/equations.yaml","algorithms/algorithms.yaml"}) {
        const auto file=directory/relative;fs::create_directories(file.parent_path());std::ofstream(file)<<"not: [valid YAML\n";
    }
    require(CaseIO::read(directory.string()).objects.size()==loaded.objects.size(),"Unregistered document was loaded");
    auto registry=doc("solvers/solvers.yaml");registry.body["unregistered"]={{"type","custom"},{"file","solvers/unregistered.yaml"}};
    save("solvers/solvers.yaml",registry);reject("registered malformed file","sequence");

    reset();auto state=doc("solvers/state.yaml");state.type="stateRegistry";save("solvers/state.yaml",state);reject("role used as file identity","expected object solver type state");
    reset();state=doc("solvers/state.yaml");state.object="state";save("solvers/state.yaml",state);reject("wrong object","expected object solver type state");
    reset();std::ofstream(directory/"solvers/state.yaml")<<"object: solver\ntype: state\nuse: [rho, rhoU, rhoE]\n";reject("missing SonicFile header","missing explicit SonicFile header");
    reset();registry=doc("solvers/solvers.yaml");registry.body["state"]["type"]="custom";save("solvers/solvers.yaml",registry);reject("reserved key wrong role","reserved registry role");
    reset();registry=doc("solvers/solvers.yaml");registry.body["otherState"]=registry.body.at("state");registry.body.erase("state");save("solvers/solvers.yaml",registry);reject("reserved role renamed","requires registry key state");
    reset();registry=doc("solvers/solvers.yaml");registry.body["State"]=registry.body.at("state");registry.body.erase("state");save("solvers/solvers.yaml",registry);reject("reserved key case mismatch","exact identity state");
    reset();registry=doc("solvers/solvers.yaml");registry.body["state"]["type"]="StateRegistry";save("solvers/solvers.yaml",registry);reject("reserved role case mismatch","reserved registry role");
    reset();registry=doc("solvers/solvers.yaml");registry.body["nested"]={{"type","registry"},{"file","solvers/nested.yaml"}};save("solvers/solvers.yaml",registry);reject("second folder registry","reserved registry identity");
    reset();registry=doc("solvers/solvers.yaml");registry.body["zAlias"]={{"type","custom"},{"file","solvers/./state.yaml"}};save("solvers/solvers.yaml",registry);reject("duplicate canonical path","duplicate referenced file");
    reset();registry=doc("solvers/solvers.yaml");registry.body["state"]["file"]="models/state.yaml";save("solvers/solvers.yaml",registry);reject("cross-object reference","must stay in solvers/");
    reset();fs::create_symlink(directory/"models/models.yaml",directory/"solvers/link.yaml");registry=doc("solvers/solvers.yaml");registry.body["state"]["file"]="solvers/link.yaml";save("solvers/solvers.yaml",registry);reject("symlink object escape","leaves its object directory");
    reset();fs::create_symlink(directory/"solvers/state.yaml",directory/"solvers/link.yaml");registry=doc("solvers/solvers.yaml");registry.body["zAlias"]={{"type","custom"},{"file","solvers/link.yaml"}};save("solvers/solvers.yaml",registry);reject("duplicate symlink","duplicate referenced file");
    reset();registry=doc("solvers/solvers.yaml");registry.body["state"]["file"]="solvers/missing.yaml";save("solvers/solvers.yaml",registry);reject("registered missing file","Cannot read");
    reset();registry=doc("solvers/solvers.yaml");registry.body.erase("state");save("solvers/solvers.yaml",registry);
    const auto withoutState=CaseIO::read(directory.string());
    require(std::none_of(withoutState.objects.begin(),withoutState.objects.end(),[](const auto& o){return o.type=="stateRegistry";}),"STATE was discovered by filename");

    reset();auto fields=doc("fields/fields.yaml");fields.body.erase("internalField");fields.body.erase("boundaries");save("fields/fields.yaml",fields);
    std::ofstream(directory/"fields/internalField.yaml")<<"invalid: [\n";
    std::ofstream(directory/"fields/boundaries.yaml")<<"invalid: [\n";
    const auto withoutData=CaseIO::read(directory.string());
    require(std::all_of(withoutData.fields.begin(),withoutData.fields.end(),[](const auto& f){return f.initial.is_null() && f.boundaries.empty();}),"Unregistered field data was loaded");

    reset();auto models=doc("models/models.yaml");models.body["surface"]={{"type","geometry"},{"file","models/surface.yaml"}};save("models/models.yaml",models);
    Serialization::SonicDocument surface;surface.object="models";surface.type="surface";surface.body={{"file","mesh/body.STL"}};save("models/surface.yaml",surface);
    const auto named=CaseIO::read(directory.string());
    require(std::any_of(named.objects.begin(),named.objects.end(),[](const auto& o){return o.name=="surface" && o.type=="geometry";}),"Named model identity lost its role");
    surface.type="geometry";save("models/surface.yaml",surface);reject("model role used as identity","expected object models type surface");
    reset();models=doc("models/models.yaml");models.body["otherState"]={{"type","stateRegistry"},{"file","models/otherState.yaml"}};save("models/models.yaml",models);reject("reserved role wrong object","belongs to object solver");
    reset();auto boundary=doc("fields/boundaries.yaml");boundary.type="boundary";save("fields/boundaries.yaml",boundary);reject("boundary role used as identity","expected object fields type boundaries");
    std::cout<<"Registry-driven identity/role IO checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
