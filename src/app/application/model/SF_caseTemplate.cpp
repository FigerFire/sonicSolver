#include "app/application/SF_inspection.h"
#include "infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.h"
#include "SF_caseTemplate.h"
#include "infrastructure/io/case/SF_case.h"
#include "infrastructure/io/mesh/SF_meshParameters.h"
#include "SF_meshGen.h"
#include <filesystem>
#include <stdexcept>
namespace SF::Application {
void createNativeCase(const std::string& directory,const std::string& recipe) {
    namespace fs=std::filesystem;using P=Model::Parameters;
    if(recipe!="compressibleSingleFluid")throw std::runtime_error("Unsupported native case recipe '"+recipe+"': no distinct validated module template; configure registered modules explicitly");
    const fs::path path=fs::absolute(directory);
    if(fs::exists(path) && !fs::is_empty(path))throw std::runtime_error("Refusing to initialize non-empty directory: "+path.string());
    Model::Description m;m.caseDir=path.string();m.name=path.filename().string();
    m.runtime={{"startTime",0},{"endTime",.01},{"endStep",0},{"CFL",.3},{"maxDeltaT",.001},
        {"writeControl","timeStep"},{"writeInterval",1},{"writeInitial",true},{"createMesh",false},
        {"output",{{"jobName",m.name},{"outputDir","result"}}},
        {"parallel",{{"enabled",false},{"nProcs",1},{"split",{1,1,1}}}}};
    m.numerics={{"time",{{"default","forwardEuler"}}},{"terms",{{"convection","teno5Steger"}}}};
    auto object=[&](std::string name,std::string category,std::string type,P parameters){Model::ObjectDescriptor o;o.name=name;o.category=category;o.type=type;o.parameters=std::move(parameters);m.objects.push_back(std::move(o));};
    object("equations","solver","equationSystem",{{"preset","compressible"}});
    object("state","solver","stateRegistry",{{"use",{"rho","rhoU","rhoE"}}});
    object("equationRegistry","solver","equationRegistry",{{"use",{"Continuity","Momentum","Energy"}}});
    object("execution","solver","algorithmRegistry",{{"Explicit",P::object()}});
    object("thermoDynamics","models","thermoDynamics",{
        {"thermoDynamics",{{"equationOfState","perfectGas"},{"thermo","hConst"},{"transport","const"}}},
        {"properties",{{"equationOfState",{{"gamma",1.4},{"R",287.05}}},{"transport",{{"mu",0.0},{"Pr",.72}}}}}});
    for(const auto& name:std::vector<std::string>{"rho","U","p"}) {
        Model::FieldDescriptor f;f.name=f.symbol=f.outputName=name;
        f.type=name=="U"?Model::ValueType::Vector:Model::ValueType::Scalar;
        f.initial=name=="U"?P::array({1.,0.,0.}):P(name=="rho"?1.:100000.);
        f.boundaries={{"Left",{{"type","zeroGradient"}}},{"Right",{{"type","zeroGradient"}}},
            {"sideWalls",{{"type","zeroGradient"}}},{"frontBack",{{"type","empty"}}}};
        m.fields.push_back(std::move(f));
    }
    m.mesh={{"files",{"mesh/mesh.sfm"}}};
    m.meshGeneration={{"scale",1.},{"nGhost",3},
        {"vertices",{{0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,.1},{2,0,.1},{2,1,.1},{0,1,.1}}},
        {"blocks",{{{"vertices",{0,1,2,3,4,5,6,7}},{"cells",{20,10,1}},{"grading",{1,1,1}}}}},
        {"edges",P::array()},{"patches",{
            {{"name","Left"},{"faces",{{0,3,7,4}}}},{{"name","Right"},{"faces",{{1,5,6,2}}}},
            {{"name","sideWalls"},{"faces",{{0,4,5,1},{3,2,6,7}}}},
            {{"name","frontBack"},{"faces",{{0,1,2,3},{4,7,6,5}}}}}}};
    MeshParameters parameters;MeshIO::decodeParameters(m.meshGeneration,parameters);
    CaseIO::write(m,path.string());
    std::vector<double>x,y,z;std::vector<int> sizes;
    MeshGen::generateStructuredMesh(parameters,x,y,z,sizes);
    if(!MeshGen::writeCombinedMesh((path/"mesh/mesh").string(),parameters,x,y,z,sizes,parameters.nGhost))
        throw std::runtime_error("Native mesh generation failed");
    const auto inspection=inspectCase(path.string());
    MultiBlockMesh mesh;
    if(!mesh.loadFiles(path.string(),inspection.config.meshFiles,inspection.mesh))throw std::runtime_error("Generated native SFM failed production reload");
}
}
