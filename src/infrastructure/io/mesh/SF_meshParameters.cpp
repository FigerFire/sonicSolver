/// @brief Native structuredMesh 参数验证；几何生成内核只接收解析后的 MeshParameters。
#include "SF_meshParameters.h"
#include "SF_meshGen.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <cmath>
#include <stdexcept>
#include <limits>
namespace SF::MeshIO {
void decodeParameters(const Model::Parameters& input,MeshParameters& out) {
    using P=Model::Parameters;
    Model::Schema{{{"scale","number"},{"nGhost","integer",true},{"vertices","array",true},
        {"blocks","array",true},{"edges","array"},{"patches","array",true}},false}.validate(input,"structuredMesh");
    MeshParameters value;value.scale=input.value("scale",1.0);value.nGhost=input.at("nGhost");
    auto finite=[](double n){if(!std::isfinite(n))throw std::runtime_error("structuredMesh requires finite coordinates/parameters");};
    auto vector=[&](const P& p){if(!p.is_array() || p.size()!=3)throw std::runtime_error("structuredMesh vector requires three components");
        Vector3 v(p.at(0).get<double>(),p.at(1).get<double>(),p.at(2).get<double>());finite(v.x);finite(v.y);finite(v.z);return v;};
    finite(value.scale);if(value.scale<=0 || value.nGhost<=0)throw std::runtime_error("structuredMesh scale/nGhost must be positive");
    for(const auto& v:input.at("vertices"))value.vertices.push_back(vector(v));
    const auto integer=[](const P& p){
        if(!p.is_number_integer() || p.get<double>()<std::numeric_limits<int>::min()
            || p.get<double>()>std::numeric_limits<int>::max())
            throw std::runtime_error("structuredMesh counts/indices require in-range integers");
        return p.get<int>();
    };
    const auto vertex=[&](int id){if(id<0 || id>=int(value.vertices.size()))throw std::runtime_error("structuredMesh invalid vertex reference");};
    for(const auto& b:input.at("blocks")) {
        Model::Schema{{{"vertices","array",true},{"cells","array",true},{"grading","array",true}},false}.validate(b,"structuredMesh block");
        if(b.at("vertices").size()!=8 || b.at("cells").size()!=3 || b.at("grading").size()!=3)throw std::runtime_error("structuredMesh hex requires 8 vertices and 3 cell/grading counts");
        MeshBlock block{};for(int i=0;i<8;++i){block.verts[i]=integer(b.at("vertices").at(i));vertex(block.verts[i]);}
        for(int i=0;i<3;++i){block.cells[i]=integer(b.at("cells").at(i));block.grade[i]=b.at("grading").at(i).get<double>();finite(block.grade[i]);if(block.cells[i]<=0 || block.grade[i]<=0)throw std::runtime_error("structuredMesh cells/grading must be positive");}
        value.blocks.push_back(block);
    }
    if(value.blocks.empty())throw std::runtime_error("structuredMesh requires at least one block");
    for(const auto& e:input.value("edges",P::array())) {
        Model::Schema{{{"type","string",true},{"vertices","array",true},{"point","array",true}},false}.validate(e,"structuredMesh edge");
        if(e.at("type")!="arc" || e.at("vertices").size()!=2)throw std::runtime_error("Unsupported structuredMesh edge; requires arc with two vertices and through point");
        MeshEdge edge;edge.type=EdgeType::ARC;edge.v0=integer(e.at("vertices").at(0));edge.v1=integer(e.at("vertices").at(1));vertex(edge.v0);vertex(edge.v1);edge.arcP=vector(e.at("point"));value.edges.push_back(edge);
    }
    for(const auto& p:input.at("patches")) {
        Model::Schema{{{"name","string",true},{"faces","array",true}},false}.validate(p,"structuredMesh patch");
        MeshFacePatch patch;patch.name=p.at("name");
        for(const auto& existing:value.facePatches)if(existing.name==patch.name)throw std::runtime_error("Duplicate structuredMesh patch");
        for(const auto& f:p.at("faces")){if(f.size()!=4)throw std::runtime_error("structuredMesh patch requires quadrilateral faces");std::array<int,4> face;for(int i=0;i<4;++i){face[i]=integer(f.at(i));vertex(face[i]);}patch.faces.push_back(face);}
        value.facePatches.push_back(patch);
    }
    out=std::move(value);
}
bool readParameters(const std::string& path,MeshParameters& out) {
    const auto doc=Serialization::readSonicFile(path);
    if(doc.object!="mesh" || doc.type!="generation")throw std::runtime_error("Expected native mesh/generation document");
    decodeParameters(doc.body,out);return true;
}
}
