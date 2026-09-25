#include "SF_serialization.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <sstream>
#include <regex>
#include <iomanip>
namespace SF::Serialization {
namespace {
Model::Parameters convert(const YAML::Node& node,const std::string& path) {
    using P=Model::Parameters;
    if(node.IsNull()) return nullptr;
    if(node.IsSequence()) { auto p=P::array();for(const auto& v:node)p.push_back(convert(v,path));return p; }
    if(node.IsMap()) { auto p=P::object();for(const auto& v:node) {
        auto k=v.first.as<std::string>();
        if(p.contains(k))throw std::runtime_error(path+":"+std::to_string(v.first.Mark().line+1)+": duplicate key "+k);
        p[k]=convert(v.second,path);
    }return p; }
    const auto s=node.Scalar();
    if(node.Tag()=="!" || node.Tag()=="tag:yaml.org,2002:str")return s;
    if(s=="true")return true;if(s=="false")return false;
    try { std::size_t n=0;auto v=std::stoll(s,&n);if(n==s.size())return v; }catch(...){}
    try { std::size_t n=0;auto v=std::stod(s,&n);if(n==s.size()) {
        if(!std::isfinite(v))throw std::runtime_error(path+": nonfinite numeric value");return v;
    }}catch(const std::invalid_argument&){}catch(const std::out_of_range&){}
    return s;
}
std::string readText(const std::string& path) {
    std::ifstream input(path);
    if(!input)throw std::runtime_error("Cannot read "+path);
    return std::string(std::istreambuf_iterator<char>(input),{});
}
/// SonicFile 在 YAML 上增加两项纯输入便利语法：
///   patch: fixedValue: [1 0 0]
///   [1 0 0]
/// 两者在送入 YAML 解析器前都被归一为标准 YAML，不改变底层通用数据模型。
std::string normalizeSonicSyntax(const std::string& text) {
    const std::regex compact(R"(^([ \t]*)([^:#\n][^:\n]*):[ \t]*([A-Za-z][A-Za-z0-9_-]*)[ \t]*:[ \t]*(.*)$)");
    const std::regex sequence(R"(\[([^\[\],#'\"]+)\])");
    std::istringstream input(text);
    std::ostringstream output;
    std::string line;
    while(std::getline(input,line)) {
        std::smatch match;
        if(std::regex_match(line,match,compact)) {
            line=match[1].str()+match[2].str()+":\n"+match[1].str()+"  "+match[3].str()+": "+match[4].str();
        }
        std::string normalized;
        std::size_t cursor=0;
        for(std::sregex_iterator it(line.begin(),line.end(),sequence),end;it!=end;++it) {
            const auto& found=*it;
            normalized.append(line,cursor,static_cast<std::size_t>(found.position())-cursor);
            const std::string contents=found[1].str();
            std::istringstream values(contents);
            std::vector<std::string> tokens;
            std::string token;
            while(values>>token)tokens.push_back(token);
            if(tokens.size()>1) {
                normalized+='[';
                for(std::size_t i=0;i<tokens.size();++i) {
                    if(i)normalized+=", ";
                    normalized+=tokens[i];
                }
                normalized+=']';
            } else normalized+=found.str();
            cursor=static_cast<std::size_t>(found.position()+found.length());
        }
        normalized.append(line,cursor,std::string::npos);
        output<<normalized<<'\n';
    }
    return output.str();
}
std::string yamlScalar(const Model::Parameters& value) {
    if(value.is_null())return "null";
    if(value.is_boolean())return value.get<bool>()?"true":"false";
    if(value.is_number_integer()||value.is_number_unsigned())return value.dump();
    if(value.is_number_float()) {
        std::ostringstream out;out<<std::setprecision(15)<<value.get<double>();return out.str();
    }
    const auto text=value.get<std::string>();
    static const std::regex plain(R"(^[A-Za-z_][A-Za-z0-9_./+-]*$)");
    return std::regex_match(text,plain)?text:value.dump();
}
std::string yamlInline(const Model::Parameters& value) {
    if(!value.is_array())return yamlScalar(value);
    std::ostringstream out;out<<'[';
    for(std::size_t i=0;i<value.size();++i) {
        if(i)out<<", ";
        if(value.at(i).is_array())out<<yamlInline(value.at(i));
        else if(value.at(i).is_object())out<<value.at(i).dump();
        else out<<yamlScalar(value.at(i));
    }
    return out.str()+']';
}
void emitMapping(std::ostream& out,const Model::Parameters& value,int indent,bool compactConditions) {
    if(!value.is_object())throw std::runtime_error("SonicFile body must be a mapping");
    const std::string pad(static_cast<std::size_t>(indent),' ');
    for(auto it=value.begin();it!=value.end();++it) {
        const auto& child=it.value();
        out<<pad<<it.key();
        if(child.is_object()) {
            if(compactConditions&&child.size()==1) {
                const auto condition=child.begin();
                if(!condition.value().is_object()) {
                    out<<": "<<condition.key();
                    if(!condition.value().is_null())out<<": "<<(condition.value().is_array()?yamlInline(condition.value()):yamlScalar(condition.value()));
                    out<<'\n';
                    continue;
                }
            }
            out<<":\n";emitMapping(out,child,indent+2,compactConditions);
        } else if(child.is_array())out<<": "<<yamlInline(child)<<'\n';
        else out<<": "<<yamlScalar(child)<<'\n';
    }
}
}
Document readDocument(const std::string& path) {
    Document document;
    const auto text=readText(path);
    {
        std::istringstream input(text);
        std::string line;
        while(std::getline(input,line)) {
            bool single=false,doubleQuote=false,escaped=false;
            for(std::size_t i=0;i<line.size();++i) {
                const char c=line[i];
                if(escaped) { escaped=false; continue; }
                if(doubleQuote&&c=='\\') { escaped=true; continue; }
                if(!doubleQuote&&c=='\'') { single=!single; continue; }
                if(!single&&c=='\"') { doubleQuote=!doubleQuote; continue; }
                if(!single&&!doubleQuote&&c=='#') {
                    document.comments.push_back(line.substr(i));
                    break;
                }
            }
        }
    }
    try { document.value=convert(YAML::Load(normalizeSonicSyntax(text)),path); return document; }
    catch(const YAML::Exception& e) { throw std::runtime_error(path+":"+std::to_string(e.mark.line+1)+": "+e.msg); }
}
Model::Parameters read(const std::string& path) { return readDocument(path).value; }
void write(const Model::Parameters& value,const std::string& path) {
    writeDocument({value,{}},path);
}
void writeDocument(const Document& document,const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write "+path);
    // JSON is a canonical, unambiguous subset of YAML 1.2.
    for(const auto& comment:document.comments)out<<comment<<'\n';
    out<<document.value.dump(2)<<'\n';if(!out)throw std::runtime_error("Failed writing "+path);
}
SonicDocument readSonicFile(const std::string& path) {
    const auto document=readDocument(path);
    if(!document.value.is_object())throw std::runtime_error(path+": SonicFile must be a mapping");
    SonicDocument result;result.source=path;result.comments=document.comments;
    auto root=document.value;
    Model::Parameters header;
    if(root.contains("SonicFile")) {
        header=root.at("SonicFile");root.erase("SonicFile");
    } else {
        header=root;
        root.erase("object");root.erase("type");
    }
    if(!header.is_object()||!header.contains("object")||!header.contains("type")
       ||!header.at("object").is_string()||!header.at("type").is_string())
        throw std::runtime_error(path+": SonicFile requires string object and type");
    result.object=header.at("object").get<std::string>();
    result.type=header.at("type").get<std::string>();
    result.body=std::move(root);
    return result;
}
void writeSonicFile(const SonicDocument& document,const std::string& path,bool compactConditions) {
    if(document.object.empty()||document.type.empty())throw std::runtime_error("SonicFile requires object and type");
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write "+path);
    for(const auto& comment:document.comments)out<<comment<<'\n';
    out<<"SonicFile:\n  object: "<<yamlScalar(document.object)<<"\n  type: "<<yamlScalar(document.type)<<"\n";
    emitMapping(out,document.body,0,compactConditions);
    if(!out)throw std::runtime_error("Failed writing "+path);
}
}
