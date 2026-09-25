#include "compatibility/SF_compatibility.h"
#include "infrastructure/io/mesh/SF_sfmMesh.h"
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <set>
namespace SF {
namespace {
using P=Model::Parameters;

P expandBoundaryDefault(P value,const std::string& caseDir,const std::vector<std::string>& meshFiles) {
    if(!value.is_object()||!value.contains("default"))return value;
    const auto fallback=value.at("default");value.erase("default");
    std::set<std::string> names;
    SFMMeshIO mesh(caseDir);
    for(const auto& file:meshFiles) {
        SFMMeshBlock block;
        if(mesh.readSingleBlock(file,block))for(const auto& [name,points]:block.pointSets){ (void)points; names.insert(name); }
        std::map<std::string,std::vector<int>> sets;
        if(mesh.readSetsOnly(file,sets))for(const auto& [name,points]:sets){ (void)points; names.insert(name); }
    }
    if(names.empty())throw std::runtime_error("boundary default requires readable mesh patch sets");
    for(const auto& name:names)if(!value.contains(name))value[name]=fallback;
    return value;
}
// Typed native sections are stored in CaseAdapter::sections_;
// native decoding lives in SF_nativeDecode.cpp.
}
CaseConfig CaseAdapter::build(const Model::Description& m) {
    native_=true;sections_={};
    EquationCompositionConfig composition;
    // 语义角色 → parameters。同一个角色只能被一个对象占用。
    std::map<std::string,P> roles;
    auto assignRole=[&](const std::string& role,const P& parameters) {
        if(!roles.emplace(role,parameters).second)
            throw std::runtime_error("Multiple model objects target the same built-in configuration: "+role);
    };
    Model::Schema{{{"startTime","number"},{"endTime","number"},{"endStep","integer"},
        {"CFL","number"},{"maxDeltaT","number"},{"writeControl","string",true},
        {"writeInterval","number",true},{"writeInitial","boolean"},{"createMesh","boolean"},
        {"parallel","object"},{"output","object",true}},false}.validate(m.runtime,"runtime.yaml");
    // 每个语义对象直接落到自己的 typed slot；不经过任何磁盘布局名字。
    sections_.runtime=m.runtime;
    sections_.output=m.runtime.at("output");sections_.hasOutput=true;
    if(m.runtime.contains("parallel")) { sections_.parallel=m.runtime.at("parallel");sections_.hasParallel=true; }
    sections_.algorithm=m.solver;sections_.hasAlgorithm=!m.solver.is_null();
    const bool hasTermRecipes=m.numerics.contains("terms");
    if(hasTermRecipes&&(m.numerics.contains("convection")||m.numerics.contains("diffusion")))
        throw std::runtime_error("numerics.yaml cannot mix terms with legacy convection/diffusion selection");
    for(auto it=m.numerics.begin();it!=m.numerics.end();++it) {
        const bool known=it.key()=="time"||it.key()=="convection"
            ||it.key()=="diffusion"||it.key()=="sources"
            ||it.key()=="transport"||it.key()=="pressureCorrection"
            ||it.key()=="terms";
        if(!known)throw std::runtime_error("numerics.yaml: unknown operator group "+it.key());
        if(it.key()=="terms") {
            if(!it.value().is_object())throw std::runtime_error("numerics.yaml: terms must be a mapping");
            for(auto term=it.value().begin();term!=it.value().end();++term)
                if(term.key()!="convection"&&term.key()!="diffusion")
                    throw std::runtime_error("numerics.yaml: unknown term recipe role "+term.key());
        }
    }
    sections_.numerics=m.numerics;
    Model::Schema{{{"files","array",true},{"generator","string"}},false}.validate(m.mesh,"mesh/mesh.yaml");
    meshFiles_=m.mesh.at("files").get<std::vector<std::string>>();
    if(m.runtime.value("createMesh",false)) {
        // Geometry format adapter: blockMeshDict is a mesh format, not the model case format.
        const std::string generator=m.mesh.at("generator");
        std::ifstream input((std::filesystem::path(caseDir_)/generator).string());
        if(!input)throw std::runtime_error("Missing mesh generator: "+generator);
        meshGenerator_=generator;
    }
    using Context=std::function<void(const std::string&,const P&)>;
    Model::FactoryRegistry<Context> factories;
    factories.add("equationRegistry",{{{"use","array",true},{"add","object"},
        {"extend","object"},{"replace","object"},{"disable","array"}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            composition.declared = true;
            composition.equations = object.parameters.at("use")
                .get<std::vector<std::string>>();
            for (const auto& equation : composition.equations) {
                if (equation == "Momentum" || equation == "Continuity"
                    || equation == "Energy") continue;
                if (equation.find('/') != std::string::npos) {
                    throw std::runtime_error(
                        "Custom equation '"+equation
                        +"' is registered but its Equation compiler is not implemented.");
                }
                throw std::runtime_error("Unknown built-in equation '"+equation+"'.");
            }
            if (object.parameters.contains("add") || object.parameters.contains("extend")
                || object.parameters.contains("replace") || object.parameters.contains("disable")) {
                throw std::runtime_error(
                    "Equation modifications are declared by the input contract, "
                    "but lowering custom modifications is not implemented.");
            }
        });
    factories.add("algorithmRegistry",{{{"Explicit","object"},{"SIMPLE","object"},
        {"PISO","object"},{"PIMPLE","object"}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            composition.declared = true;
            int selected = 0;
            const auto select = [&](const char* name) {
                if (!object.parameters.contains(name)) return;
                if (++selected != 1) throw std::runtime_error(
                    "algorithms.yaml must select exactly one algorithm preset.");
                composition.algorithm = name;
                const auto& parameters = object.parameters.at(name);
                if (!parameters.is_object()) throw std::runtime_error(
                    std::string("Algorithm '")+name+"' must be a mapping.");
                composition.outerCorrectors = parameters.value("outerCorrectors",1);
                composition.pressureCorrectors = parameters.value(
                    composition.algorithm == "PISO" ? "correctors" : "pressureCorrectors",1);
                composition.nonOrthogonalCorrectors = parameters.value(
                    "nonOrthogonalCorrectors",0);
                if (composition.outerCorrectors <= 0 || composition.pressureCorrectors <= 0
                    || composition.nonOrthogonalCorrectors < 0) throw std::runtime_error(
                    std::string("Algorithm '")+name+"' has invalid corrector counts.");
            };
            select("Explicit");
            select("SIMPLE");
            select("PISO");
            select("PIMPLE");
            if (selected == 0) throw std::runtime_error(
                "algorithms.yaml must select one registered algorithm preset.");
        });
    // Native object type → typed semantic slot. The sink name is the model's own
    // semantic role, never a disk layout path.
    static const std::map<std::string,std::string> typedRoles{
        {"Thermophysical","Thermophysical"},{"PhaseSystem","PhaseSystem"},
        {"PhaseChange","PhaseChange"},{"Turbulence","Turbulence"},
        {"ImmersedBoundary","ImmersedBoundary"},{"ILW","ILW"},
        {"Gravity","Gravity"},{"MRF","MRF"},{"WallHeat","WallHeat"},
        {"thermophysical","Thermophysical"},{"phaseChange","PhaseChange"},
        {"turbulence","Turbulence"},{"IBM","ImmersedBoundary"},
        {"gravity","Gravity"},{"wallHeat","WallHeat"}};
    for(const auto& role:typedRoles)
        factories.add(role.first,{{},true},[role](const Model::ObjectDescriptor& o,Context& emit){emit(role.second,o.parameters);});
    factories.add("thermoDynamics",{{{"thermoDynamics","object",true},
        {"properties","object",true}},false},[&](const Model::ObjectDescriptor& object,Context& emit) {
            const auto& selection = object.parameters.at("thermoDynamics");
            if (selection.contains("energy")) throw std::runtime_error(
                "thermoDynamics.yaml: energy selects an equation and is forbidden; use equations.yaml.");
            Model::Schema{{{"equationOfState","string",true},{"thermo","string"},
                {"transport","string"}},false}.validate(selection,"thermoDynamics");
            composition.declared = true;
            composition.thermoDynamics.equationOfState =
                selection.at("equationOfState").get<std::string>();
            if (selection.contains("thermo")) composition.thermoDynamics.thermo =
                selection.at("thermo").get<std::string>();
            if (selection.contains("transport")) composition.thermoDynamics.transport =
                selection.at("transport").get<std::string>();
            const auto& properties = object.parameters.at("properties");
            if (composition.thermoDynamics.equationOfState == "rhoConst") {
                if (!properties.contains("equationOfState")
                    || !properties.at("equationOfState").contains("rho")) throw std::runtime_error(
                    "rhoConst requires properties.equationOfState.rho.");
                composition.thermoDynamics.constantDensity = properties.at("equationOfState")
                    .at("rho").get<double>();
                if (!(composition.thermoDynamics.constantDensity > 0.0)) {
                    throw std::runtime_error("rhoConst density must be positive.");
                }
            }
            emit("Thermophysical",object.parameters);
        });
    factories.add("multiPhase",{{},true},[](const Model::ObjectDescriptor& o,Context& emit) {
        auto parameters=o.parameters;
        if(!o.selection.empty()&&!parameters.contains("phaseSystem")) {
            parameters["phaseSystem"]=o.selection=="eulerian-eulerian"?"eulerianEulerian":o.selection;
        }
        emit("PhaseSystem",parameters);
    });
    // Equation expressions are deserialized generically. This solver only executes its registered presets.
    factories.add("BuiltinEquationSystem",{{{"preset","string",true}},false},[](const Model::ObjectDescriptor& o,Context&) {
        const auto preset=o.parameters.at("preset").get<std::string>();
        if(preset!="compressible"&&preset!="eulerianEulerian")throw std::runtime_error("Unknown equation preset: "+preset);
        if(!o.expression.empty())throw std::runtime_error("Custom equation expressions require an Equation compiler; built-in presets cannot ignore an expression.");
    });
    factories.add("equationSystem",{{{"preset","string",true}},false},[](const Model::ObjectDescriptor& o,Context&) {
        const auto preset=o.parameters.at("preset").get<std::string>();
        if(preset!="compressible"&&preset!="eulerianEulerian")throw std::runtime_error("Unknown equation preset: "+preset);
    });
    factories.add("File",{{{"file","string",true}},false},[&](const Model::ObjectDescriptor& o,Context&) {
        if(!std::filesystem::exists(std::filesystem::path(caseDir_)/o.parameters.at("file").get<std::string>()))throw std::runtime_error("Missing geometry: "+o.name);
    });
    factories.add("geometry",{{{"file","string",true}},false},[&](const Model::ObjectDescriptor& o,Context&) {
        if(!std::filesystem::exists(std::filesystem::path(caseDir_)/o.parameters.at("file").get<std::string>()))throw std::runtime_error("Missing geometry: "+o.name);
    });
    Context emit=assignRole;
    for(const auto& o:m.objects)factories.create(o,emit);
    for(const auto& role:roles) {
        const std::string& name=role.first;
        const P& value=role.second;
        if(name=="Thermophysical"){sections_.thermoDynamics=value;sections_.hasThermoDynamics=true;}
        else if(name=="PhaseSystem"){sections_.phaseSystem=value;sections_.hasPhaseSystem=true;}
        else if(name=="PhaseChange"){sections_.phaseChange=value;sections_.hasPhaseChange=true;}
        else if(name=="Turbulence"){sections_.turbulence=value;sections_.hasTurbulence=true;}
        else if(name=="ImmersedBoundary"){sections_.ibm=value;sections_.hasIbm=true;}
        else if(name=="ILW"){sections_.ilw=value;sections_.hasIlw=true;}
        else if(name=="Gravity"){sections_.gravity=value;sections_.hasGravity=true;}
        else if(name=="MRF"){sections_.mrf=value;sections_.hasMrf=true;}
        else if(name=="WallHeat"){sections_.wallHeat=value;sections_.hasWallHeat=true;}
        else throw std::runtime_error("Unknown built-in configuration role: "+name);
    }
    // field 的语义名字就是它的解码入口；不再合成 <startTime>/<name> 伪文档。
    for(const auto& f:m.fields) {
        if(f.initial.is_null())continue;
        Model::FieldDescriptor field=f;
        field.boundaries=expandBoundaryDefault(f.boundaries,caseDir_,meshFiles_);
        sections_.fields.push_back(std::move(field));
    }
    auto result=decodeNativeCase();if(!result)throw std::runtime_error("Cannot bind model "+m.name);
    result->caseName=m.name;
    if(result->createMesh)result->meshParameterFile=meshGenerator_;
    result->modelDescription=std::make_shared<const Model::Description>(m);
    result->composition=std::move(composition);
    return *result;
}
}
