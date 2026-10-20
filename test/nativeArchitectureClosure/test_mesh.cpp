#include "infrastructure/io/mesh/SF_meshParameters.h"
#include "SF_meshGen.h"
#include <stdexcept>
#include <iostream>
using namespace SF;
int main(){
 using P=Model::Parameters;
 P p={{"nGhost",3},{"scale",1.},{"vertices",{{0,0,0},{2,0,0},{2,1,0},{0,1,0},{0,0,.1},{2,0,.1},{2,1,.1},{0,1,.1}}},
    {"blocks",{{{"vertices",{0,1,2,3,4,5,6,7}},{"cells",{20,10,1}},{"grading",{1,1,1}}}}},
    {"patches",{{{"name","left"},{"faces",{{0,3,7,4}}}}}}};
 MeshParameters decoded;MeshIO::decodeParameters(p,decoded);
 std::vector<double>x,y,z;std::vector<int> sizes;MeshGen::generateStructuredMesh(decoded,x,y,z,sizes);
 if(x.size()!=21*11*2 || sizes!=std::vector<int>{21,11,2})throw std::runtime_error("Native mesh cell/node counts changed");
 for(int k=0;k<2;++k)for(int j=0;j<11;++j)for(int i=0;i<21;++i){int n=(k*11+j)*21+i;
  if(std::abs(x[n]-.1*i)>1e-14 || std::abs(y[n]-.1*j)>1e-14 || std::abs(z[n]-.1*k)>1e-14)throw std::runtime_error("Native mesh coordinates changed");}
 for(int n=0;n<7;++n){P bad=p;
  if(n==0)bad["blocks"][0]["cells"][0]=20.5;
  if(n==1)bad["blocks"][0]["vertices"][0]=100;
  if(n==2)bad["blocks"][0]["grading"][0]=0;
  if(n==3)bad["nGhost"]=3.5;
  if(n==4)bad["unregistered"]=1;
  if(n==5)bad["patches"].push_back(bad["patches"][0]);
  if(n==6)bad.erase("patches");
  bool rejected=false;try{MeshIO::decodeParameters(bad,decoded);}catch(const std::exception&){rejected=true;}
  if(!rejected)throw std::runtime_error("Invalid native mesh parameter silently accepted");
 }
 std::cout<<"native mesh coordinates/counts and 7 typed parameter negatives passed\n";
}
