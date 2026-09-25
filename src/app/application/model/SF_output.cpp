#include "SF_output.h"
#include "core/field/SF_field.h"
#include <map>
#include <stdexcept>
namespace SF::Application::ModelLoader {
namespace {
using Array=ResultWriter::ScalarField;
std::vector<Array> fluidFields(const Field& field,double gamma) {
    Array density{"Density",[&field](int i,int j,int k) {
        return field.hasStateModel()?field.thermodynamicState(i,j,k).density:field(i,j,k,RHO);
    }};
    Array velocity;velocity.name="Velocity";velocity.components=3;
    velocity.componentAt=[&field](int i,int j,int k,int c) {
        if(field.hasStateModel())return field.thermodynamicState(i,j,k).velocity[c];
        return field(i,j,k,RU+c)/field(i,j,k,RHO);
    };
    Array pressure{"Pressure",[&field,gamma](int i,int j,int k) {
        if(field.hasStateModel())return field.thermodynamicState(i,j,k).pressure;
        const double r=field(i,j,k,RHO);
        const double kinetic=0.5*(field(i,j,k,RU)*field(i,j,k,RU)+field(i,j,k,RV)*field(i,j,k,RV)+field(i,j,k,RW)*field(i,j,k,RW))/r;
        return (gamma-1)*(field(i,j,k,E)-kinetic);
    }};
    return {density,velocity,pressure};
}
}
void configureOutput(ResultWriter& writer,const CaseConfig& config) {
    const auto description=config.modelDescription;
    const double gamma=config.solver.numerics.idealGasGamma;
    writer.setFieldProvider([description,gamma](const Field& field,const std::vector<Array>& extra) {
        auto arrays=fluidFields(field,gamma);
        arrays.insert(arrays.end(),extra.begin(),extra.end());
        if(!description)return arrays;
        std::map<std::string,Array> available;
        for(const auto& a:arrays)available[a.name]=a;
        const std::map<std::string,std::string> aliases{{"rho","Density"},{"U","Velocity"},{"p","Pressure"},
            {"T","Temperature"},{"k","TurbulenceK"},{"epsilon","TurbulenceEpsilon"},{"omega","TurbulenceOmega"}};
        for(const auto& f:description->fields) {
            auto it=available.find(f.name);
            const auto alias=aliases.find(f.name);
            if(!f.output) { it=available.end(); }
            else {
                if(it==available.end()&&alias!=aliases.end())it=available.find(alias->second);
                if(it==available.end())it=available.find(f.outputName);
            }
            Array selected;
            if(it!=available.end()) {selected=it->second;available.erase(it);}
            else if(f.output) {
                // Eulerian models attach a phase suffix (for example k.liquid).
                // An unqualified descriptor is accepted only when it resolves uniquely.
                std::vector<std::string> matches;
                const auto prefix=f.name+".";
                for(const auto& [name,array]:available)
                    if(name.rfind(prefix,0)==0)matches.push_back(name);
                if(matches.size()>1) {
                    // Built-in unqualified aliases such as T can be a boundary/initial
                    // descriptor while the runtime exposes only phase-qualified arrays.
                    // The qualified descriptors will emit those arrays below; a custom
                    // ambiguous descriptor remains an error.
                    if(alias!=aliases.end())continue;
                    throw std::runtime_error("Ambiguous phase field provider for output field: "+f.name);
                }
                if(matches.size()==1) {
                    selected=available.at(matches.front());
                    available.erase(matches.front());
                }
            }
            if(f.output && !selected.valueAt && !selected.componentAt
                    && f.type==Model::ValueType::Vector&&f.name.rfind("U.",0)==0) {
                const auto suffix=f.name.substr(1);
                std::vector<Array> components;
                for(const auto& axis:{"Ux","Uy","Uz"}) {
                    const auto component=available.find(std::string(axis)+suffix);
                    if(component==available.end())throw std::runtime_error("Missing phase velocity provider: "+f.name);
                    components.push_back(component->second);
                    available.erase(component);
                }
                selected.components=3;selected.componentAt=[components](int i,int j,int k,int c){return components.at(c).value(i,j,k,0);};
            } else if(!selected.valueAt && !selected.componentAt && f.output) {
                throw std::runtime_error(
                    "No runtime field provider registered for output field: "+f.name);
            }
            if(!f.output)continue;
            selected.name=f.outputName;
            if(selected.components!=Model::components(f.type))throw std::runtime_error("Output provider type mismatch: "+f.name);
            available[selected.name]=selected;
        }
        // Every array passes through the same registry. Models retain ownership of their storage.
        Model::FieldRegistry registry;
        for(const auto& [name,a]:available) {
            Model::FieldDescriptor d;d.name=name;d.outputName=name;
            d.type=a.components==1?Model::ValueType::Scalar:a.components==3?Model::ValueType::Vector:a.components==6?Model::ValueType::SymmetricTensor:Model::ValueType::Tensor;
            registry.add({d,static_cast<std::size_t>(field.TotalSize()),[a=Array(a),&field](std::size_t index,int c) {
                int i,j,k;field.getIJK(static_cast<int>(index),i,j,k);return a.value(i,j,k,c);
            },{}});
        }
        std::vector<Array> result;
        for(const auto& [name,view]:registry.fields()) {
            Array a;a.name=name;a.components=Model::components(view.descriptor.type);
            a.writesTemperatureBoundaryId=available.at(name).writesTemperatureBoundaryId;
            a.componentAt=[read=view.read,&field](int i,int j,int k,int c){return read(field.getIdx(i,j,k),c);};
            result.push_back(std::move(a));
        }
        return result;
    });
}
}
