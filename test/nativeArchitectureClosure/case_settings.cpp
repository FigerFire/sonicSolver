#include "infrastructure/io/case/SF_case.h"
#include <iostream>
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("Expected case directory");
 const auto m=SF::CaseIO::read(argv[1]);SF::Model::Parameters out={{"name",m.name},{"runtime",m.runtime},{"numerics",m.numerics},{"algorithm",m.solver},{"mesh",m.mesh},{"generation",m.meshGeneration}};
 out["fields"]=SF::Model::Parameters::object();
 for(const auto& f:m.fields)out["fields"][f.name]={{"symbol",f.symbol},{"shape",f.type==SF::Model::ValueType::Vector?"vector":"scalar"},{"initial",f.initial},{"sets",f.sets},{"boundaries",f.boundaries},{"output",f.output},{"outputName",f.outputName}};
 out["registeredObjects"]=SF::Model::Parameters::object();
 for(const auto& o:m.objects)out["registeredObjects"][o.category+"/"+o.name]={{"role",o.type},{"selection",o.selection},{"parameters",o.parameters}};
 std::cout<<out.dump(2)<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
