#pragma once

/// @file SF_signature.h
/// @brief 与存储和离散解耦的数学类型、算子身份及构造期规则。
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <sstream>
#include <cstdint>
#include <memory>
#include <iomanip>

namespace SF::Math {
/// @brief SI顺序：质量、长度、时间、温度、物质的量、电流、发光强度。未声明不是无量纲。
struct Dimensions {
    std::array<double,7> exponent{};
    bool known=false;
    static Dimensions specified(std::array<double,7> value) {
        for(double x:value)if(!std::isfinite(x))throw std::runtime_error("Non-finite dimension exponent");
        return {value,true};
    }
    static Dimensions dimensionless() {return specified({});}
};
inline void requireKnown(const Dimensions& value) {
    for(double x:value.exponent)if(!std::isfinite(x))throw std::runtime_error("Non-finite mathematical dimension");
    if(!value.known)throw std::runtime_error("Unspecified mathematical dimensions");
}
inline bool equivalent(const Dimensions& a,const Dimensions& b) {
    requireKnown(a);requireKnown(b);
    for(std::size_t i=0;i<7;++i)if(std::abs(a.exponent[i]-b.exponent[i])>1.e-12)return false;
    return true;
}
inline Dimensions scaled(Dimensions value,double factor) {
    requireKnown(value);
    for(double& x:value.exponent)x*=factor;
    return Dimensions::specified(value.exponent);
}
inline Dimensions combine(Dimensions a,const Dimensions& b,double sign=1) {
    requireKnown(a);requireKnown(b);
    for(std::size_t i=0;i<7;++i)a.exponent[i]+=sign*b.exponent[i];
    return Dimensions::specified(a.exponent);
}
/// @brief Full张量每个指标的extent相同；Symmetric仅用于二阶对称张量。
enum class TensorStructure { Full, Symmetric };
struct TensorSignature {
    int rank=0;
    Dimensions dimensions;
    int extent=3;
    TensorStructure structure=TensorStructure::Full;
};
/// @brief 只有异质product才使用tuple；它不冒充同量纲scalar/tensor。
struct Signature : TensorSignature {
    std::shared_ptr<const std::array<TensorSignature,4>> tuple;
    unsigned tupleSize=0;
    bool callable=false;
    Signature()=default;
    Signature(int r,Dimensions d):TensorSignature{r,d} {}
};
inline void requireKnown(const Signature& s) {
    if(s.tupleSize) {
        if(s.tupleSize>4 || !s.tuple)throw std::runtime_error("Oversized mathematical product");
        for(unsigned i=0;i<s.tupleSize;++i) {requireKnown((*s.tuple)[i].dimensions);if((*s.tuple)[i].rank<0)throw std::runtime_error("Invalid product member");}
        return;
    }
    requireKnown(s.dimensions);
    if(s.rank<0 || s.extent<=0 || (s.structure==TensorStructure::Symmetric && s.rank!=2))
        throw std::runtime_error("Invalid mathematical tensor signature");
}
inline bool compatible(const Signature& a,const Signature& b) {
    requireKnown(a);requireKnown(b);
    if(a.tupleSize || b.tupleSize) {
        if(a.tupleSize!=b.tupleSize || a.callable!=b.callable)return false;
        for(unsigned i=0;i<a.tupleSize;++i) {
            const auto& x=(*a.tuple)[i];const auto& y=(*b.tuple)[i];
            if(x.rank!=y.rank || x.extent!=y.extent || x.structure!=y.structure || !equivalent(x.dimensions,y.dimensions))return false;
        }return true;
    }
    return a.rank==b.rank && a.extent==b.extent && a.structure==b.structure && equivalent(a.dimensions,b.dimensions);
}
inline Signature signature(int rank,std::array<double,7> dimensions) {return {rank,Dimensions::specified(dimensions)};}
inline std::string text(const Signature& a) {
    std::ostringstream out;out<<std::setprecision(17);if(a.tupleSize){out<<(a.callable?"function(inputs..., result)(":"product(");for(unsigned i=0;i<a.tupleSize;++i){Signature member;static_cast<TensorSignature&>(member)=(*a.tuple)[i];out<<text(member)<<"; ";}return out.str()+")";}out<<"rank="<<a.rank<<" extent="<<a.extent<<" structure="<<static_cast<int>(a.structure)<<" dimensions=";
    if(!a.dimensions.known)out<<"Unspecified";else {out<<"[";for(auto x:a.dimensions.exponent)out<<x<<" ";out<<"]";}return out.str();
}
inline Signature product(std::initializer_list<Signature> members) {
    if(members.size()>4 || members.size()==0)throw std::runtime_error("Mathematical product requires 1..4 members");
    Signature result;result.rank=-1;std::array<TensorSignature,4> values{};
    for(const auto& value:members){requireKnown(value);if(value.tupleSize)throw std::runtime_error("Nested product unsupported");values[result.tupleSize++]=value;}
    result.tuple=std::make_shared<const std::array<TensorSignature,4>>(values);return result;
}
/// @brief 只声明数学函数输入和结果；复合元数据不可变，不对应Field storage。
inline Signature function(Signature result,std::initializer_list<Signature> inputs) {
    if(inputs.size()>3)throw std::runtime_error("This function signature supports at most 3 inputs");
    Signature value;value.rank=-1;value.callable=true;std::array<TensorSignature,4> types{};
    for(const auto& input:inputs){requireKnown(input);if(input.tupleSize)throw std::runtime_error("Nested function arguments unsupported");types[value.tupleSize++]=input;}
    requireKnown(result);if(result.tupleSize)throw std::runtime_error("Function result must be a tensor quantity");
    types[value.tupleSize++]=result;value.tuple=std::make_shared<const std::array<TensorSignature,4>>(types);return value;
}
inline void requireCompatible(const Signature& a,const Signature& b,std::string_view operation) {
    if(!compatible(a,b))throw std::runtime_error("Incompatible mathematical signatures for "+std::string(operation)+": operand 0 "+text(a)+"; operand 1 "+text(b));
}
/// @brief Identity不含occurrence或数值provider。扩展规则使用独立的显式名称和签名函数。
enum class Op { Add, Subtract, Multiply, Divide, Negate, Dot, Outer, Cross, DoubleDot,
    Grad, Div, Curl, Laplacian, Diffusion, Ddt, Partial, Trace, Pow, Sin, Source, Extension };
inline const char* name(Op op) {
    switch(op) {
#define SF_MATH_NAME(ID,TEXT) case Op::ID:return TEXT;
        SF_MATH_NAME(Add,"add") SF_MATH_NAME(Subtract,"subtract") SF_MATH_NAME(Multiply,"multiply")
        SF_MATH_NAME(Divide,"divide") SF_MATH_NAME(Negate,"negate") SF_MATH_NAME(Dot,"dot")
        SF_MATH_NAME(Outer,"outer") SF_MATH_NAME(Cross,"cross") SF_MATH_NAME(DoubleDot,"doubleDot") SF_MATH_NAME(Curl,"curl") SF_MATH_NAME(Grad,"grad")
        SF_MATH_NAME(Div,"div") SF_MATH_NAME(Laplacian,"laplacian") SF_MATH_NAME(Diffusion,"diffusion")
        SF_MATH_NAME(Ddt,"ddt") SF_MATH_NAME(Partial,"partial") SF_MATH_NAME(Trace,"trace")
        SF_MATH_NAME(Pow,"pow") SF_MATH_NAME(Sin,"sin") SF_MATH_NAME(Source,"source")
        SF_MATH_NAME(Extension,"extension")
#undef SF_MATH_NAME
    }
    throw std::runtime_error("Unknown mathematical operator identity");
}
inline Op resolveOperator(std::string_view text) {
    if(text=="gradient")return Op::Grad;
    for(auto op:{Op::Add,Op::Subtract,Op::Multiply,Op::Divide,Op::Negate,Op::Dot,Op::Outer,Op::Cross,Op::DoubleDot,
            Op::Grad,Op::Div,Op::Curl,Op::Laplacian,Op::Diffusion,Op::Ddt,Op::Partial,Op::Trace,Op::Pow,Op::Sin,Op::Source})
        if(text==name(op))return op;
    throw std::runtime_error("Unknown mathematical operator: "+std::string(text));
}
/// @brief 只进行数学推导；指数参数仅对常数幂有效，Numerics不参与此规则。
inline Signature resultSignature(Op op,const std::vector<Signature>& a,double exponent=1);
namespace detail {
inline Signature inferSignature(Op op,const std::vector<Signature>& a,double exponent) {
    if(op==Op::Source && !a.empty() && a[0].callable) {
        requireKnown(a[0]);
        if(a.size()!=a[0].tupleSize)throw std::runtime_error("Source function argument count mismatch");
        for(unsigned i=1;i<a.size();++i){Signature input;static_cast<TensorSignature&>(input)=(*a[0].tuple)[i-1];requireCompatible(a[i],input,"source function argument");}
        Signature result;static_cast<TensorSignature&>(result)=(*a[0].tuple)[a[0].tupleSize-1];return result;
    }
    if(op==Op::Pow && a.size()==2) {
        requireCompatible(a[1],signature(0,{}),"constant power exponent");
        return resultSignature(op,{a[0]},exponent);
    }
    if(op==Op::Div && a.size()==2) {
        return resultSignature(Op::Div,{resultSignature(Op::Multiply,a)});
    }
    const auto count=(op==Op::Add || op==Op::Subtract || op==Op::Multiply || op==Op::Divide
        || op==Op::Dot || op==Op::Outer || op==Op::Cross || op==Op::DoubleDot || op==Op::Diffusion)?2u:1u;
    if(a.size()!=count)throw std::runtime_error("Wrong operand count for "+std::string(name(op)));
    for(const auto& s:a){requireKnown(s);if(s.tupleSize)throw std::runtime_error("Tensor operator cannot act on heterogeneous product");}
    auto result=a[0];
    const auto length=Dimensions::specified({0,1,0,0,0,0,0});
    const auto time=Dimensions::specified({0,0,1,0,0,0,0});
    switch(op) {
    case Op::Add:case Op::Subtract:requireCompatible(a[0],a[1],name(op));return result;
    case Op::Negate:case Op::Source:return result;
    case Op::Multiply:
        if(a[0].rank && a[1].rank)throw std::runtime_error("Tensor multiplication requires explicit dot/outer/contraction");
        result=a[a[0].rank?0:1];result.dimensions=combine(a[0].dimensions,a[1].dimensions);break;
    case Op::Divide:
        if(a[1].rank)throw std::runtime_error("Mathematical division requires scalar denominator");
        result.dimensions=combine(a[0].dimensions,a[1].dimensions,-1);break;
    case Op::Outer:case Op::Dot:case Op::Cross:case Op::DoubleDot:
        if(a[0].extent!=a[1].extent)throw std::runtime_error("Incompatible contraction extents");
        if(op!=Op::Outer && (!a[0].rank || !a[1].rank))throw std::runtime_error("Contraction requires tensor operands");
        if(op==Op::Cross && (a[0].rank!=1 || a[1].rank!=1 || a[0].extent!=3))throw std::runtime_error("Cross requires two 3D vectors");
        if(op==Op::DoubleDot && (a[0].rank<2 || a[1].rank<2))throw std::runtime_error("Double contraction requires rank >= 2 operands");
        result.rank=op==Op::Cross?1:a[0].rank+a[1].rank-(op==Op::Dot?2:op==Op::DoubleDot?4:0);
        result.structure=TensorStructure::Full;result.dimensions=combine(a[0].dimensions,a[1].dimensions);break;
    case Op::Grad:case Op::Div:case Op::Partial:
        if(op==Op::Div && result.rank==0)throw std::runtime_error("Continuous div requires tensor rank >= 1");
        result.rank+=op==Op::Grad?1:op==Op::Div?-1:0;result.structure=TensorStructure::Full;
        result.dimensions=combine(result.dimensions,length,-1);break;
    case Op::Curl:
        if(result.rank!=1 || result.extent!=3)throw std::runtime_error("Curl requires a 3D vector");
        result.dimensions=combine(result.dimensions,length,-1);break;
    case Op::Laplacian:result.dimensions=combine(result.dimensions,scaled(length,2),-1);break;
    case Op::Diffusion:
        if(a[0].rank)throw std::runtime_error("Diffusion coefficient must be scalar");
        result=a[1];result.dimensions=combine(combine(a[0].dimensions,a[1].dimensions),scaled(length,2),-1);break;
    case Op::Ddt:result.dimensions=combine(result.dimensions,time,-1);break;
    case Op::Trace:
        if(result.rank!=2)throw std::runtime_error("Trace requires rank 2");
        result.rank=0;result.structure=TensorStructure::Full;break;
    case Op::Pow:
        if(result.rank || !std::isfinite(exponent))throw std::runtime_error("Power requires scalar and finite constant exponent");
        result.dimensions=scaled(result.dimensions,exponent);break;
    case Op::Sin:
        if(result.rank || !equivalent(result.dimensions,Dimensions::dimensionless()))throw std::runtime_error("Sin requires dimensionless scalar");break;
    default:throw std::runtime_error("Unknown mathematical operator identity");
    case Op::Extension:throw std::runtime_error("Extension requires its declared mathematical signature rule");
    }
    requireKnown(result);return result;
}
} // namespace detail
inline Signature resultSignature(Op op,const std::vector<Signature>& operands,double exponent) {
    try {return detail::inferSignature(op,operands,exponent);}
    catch(const std::runtime_error& error) {
        std::string label;try{label=name(op);}catch(const std::runtime_error&){label="<unknown identity>";}
        std::ostringstream out;out<<"Math operator "<<label<<": ";
        for(std::size_t i=0;i<operands.size();++i)out<<"operand "<<i<<" {"<<text(operands[i])<<"}; ";
        throw std::runtime_error(out.str()+"reason: "+error.what());
    }
}
/// @brief 扩展的稳定身份；名称只在声明边界映射，编译器比较此值。
constexpr std::uint64_t extensionIdentity(std::string_view value) {
    std::uint64_t hash=14695981039346656037ull;
    for(char c:value){hash^=static_cast<unsigned char>(c);hash*=1099511628211ull;}return hash;
}
/// @brief 模块拥有的静态数学规则；不登记数值算法，不使用全局可变注册状态。
struct OperatorRule {
    const char* id;
    Signature (*result)(const std::vector<Signature>&);
};
} // namespace SF::Math
