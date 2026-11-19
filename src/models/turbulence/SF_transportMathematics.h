#pragma once

/// @brief Single-fluid RAS mathematical contributions.
/// The stage recipe declares full conservative transport. The explicit legacy
/// selection retains its real no-convection, once-step arithmetic contract.
#include "core/system/SF_formula.h"
#include <stdexcept>

namespace SF::Turbulence {
/// @brief 原 stage SST 固定常量，供数学声明与未迁移闭合使用。
inline constexpr double stageSstBetaStar=0.09;
inline constexpr double stageSstSigmaK1=.85,stageSstSigmaK2=1.;
inline constexpr double stageSstSigmaW1=.5,stageSstSigmaW2=.856;
inline constexpr double stageSstBeta1=.075,stageSstBeta2=.0828;
inline Math::Signature mathematicalSignature(std::string name) {
    const auto scalar=[](std::array<double,7> d){return Math::signature(0,d);};
    const auto prefix=name.find('.');
    if(prefix!=std::string::npos && (name.substr(0,prefix)=="SST" || name.substr(0,prefix)=="kOmegaSST" || name.substr(0,prefix)=="kEpsilon")) {
        name=name.substr(prefix+1);name=name.substr(0,name.find('.',name.rfind("nu.",0)==0 || name.rfind("production.",0)==0 || name.rfind("sink.",0)==0 ? name.find('.')+1:0));
        if(name=="nu.k" || name=="nu.omega")return scalar({0,2,-1,0,0,0,0});
        if(name=="production.k" || name=="limitedProduction")return scalar({1,-1,-3,0,0,0,0});
        if(name=="production.omega" || name=="cross" || name=="crossSource" || name=="crossDiffusion")return scalar({1,-3,-2,0,0,0,0});
        if(name=="sink.k" || name=="sink.omega")return scalar({1,-3,-1,0,0,0,0});
        if(name=="omegaFloor")return scalar({0,0,-1,0,0,0,0});
        if(name=="maxRhoK")return scalar({1,-1,-2,0,0,0,0});
        if(name=="maxK")return scalar({0,2,-2,0,0,0,0});
        if(name=="maxMuT" || name=="laminarMu" || name=="eddyDynamicViscosity")return scalar({1,-1,-1,0,0,0,0});
        for(const auto* key:{"betaStar","beta","gamma","C1","C2","sigmaK","sigmaOmega","sigmaEpsilon","F1","sigmaK1","sigmaK2","sigmaW1","sigmaW2","beta1","beta2","gamma1","gamma2","a1"})if(name==key)return scalar({});
    } else {
        name=name.substr(0,prefix);
        if(name=="rho" || name=="phaseMass")return scalar({1,-3,0,0,0,0,0});
        if(name=="mu" || name=="mu_t")return scalar({1,-1,-1,0,0,0,0});
        if(name=="k")return scalar({0,2,-2,0,0,0,0});
        if(name=="omega")return scalar({0,0,-1,0,0,0,0});
        if(name=="epsilon")return scalar({0,2,-3,0,0,0,0});
        if(name=="alpha")return scalar({});
        if(name=="massFaceFlux")return scalar({1,0,-1,0,0,0,0});
        if(name=="U")return Math::signature(1,{0,1,-1,0,0,0,0});
        if(name=="rhoU")return Math::signature(1,{1,-2,-1,0,0,0,0});
        if(name=="Smagorinsky")return scalar({1,-1,-1,0,0,0,0});
    }
    throw std::runtime_error("Unknown declared turbulence mathematical symbol: "+name);
}
inline System::FormulaExpr mathematicalSymbol(const std::string& name) {return System::FormulaExpr::symbol(name,mathematicalSignature(name));}
/// Integrated mass rate divided by cell volume is distinct from continuous div.
inline System::FormulaExpr faceMassDivergence(System::FormulaExpr value) {
    static const Math::OperatorRule rule{"faceMassDivergence",[](const std::vector<Math::Signature>& a) {
        if(a.size()!=1 || a[0].rank)throw std::runtime_error("Face mass divergence requires scalar face rate");
        auto result=a[0];result.dimensions=Math::combine(result.dimensions,Math::Dimensions::specified({0,3,0,0,0,0,0}),-1);return result;
    }};return System::FormulaExpr::op(rule,{std::move(value)});
}
}

namespace SF::Turbulence {
/// @brief 原局部销毁率；执行频率和目标由 HOW/WHICH 独立绑定。
inline System::Equation stageKDestructionMathematics() {
    using E=System::FormulaExpr;
    return {"SST.sink.k",mathematicalSymbol("SST.sink.k"),
        E::multiply(E::multiply(mathematicalSymbol("kOmegaSST.betaStar"),mathematicalSymbol("rho")),mathematicalSymbol("omega")),
        {System::OriginKind::Model,"kOmegaSST.stageConservative"}};
}
/// @brief 原 stage 闭合的局部关系；非局部 F1/mu_t 由既有 prepare 发布。
inline std::vector<System::Equation> stageLocalMathematics(
        const System::FormulaExpr& viscosity=mathematicalSymbol("mu_t")) {
    using E=System::FormulaExpr;
    const auto s=[](const std::string& n){return mathematicalSymbol(n);};
    const auto blend=[&](const char* first,const char* second) {
        return E::add(E::multiply(s("SST.F1"),s(first)),
            E::multiply(E::subtract(E::constantValue(1),s("SST.F1")),s(second)));
    };
    const auto diffusion=[&](const char* target,const char* first,const char* second) {
        return System::Equation{target,s(target),E::divide(E::add(s("kOmegaSST.laminarMu"),
            E::multiply(blend(first,second),viscosity)),s("rho")),
            {System::OriginKind::Model,"kOmegaSST.stageConservative"}};
    };
    return {stageKDestructionMathematics(),
        diffusion("SST.nu.k","kOmegaSST.sigmaK1","kOmegaSST.sigmaK2"),
        diffusion("SST.nu.omega","kOmegaSST.sigmaW1","kOmegaSST.sigmaW2"),
        {"SST.sink.omega",s("SST.sink.omega"),
            E::multiply(E::multiply(blend("kOmegaSST.beta1","kOmegaSST.beta2"),s("rho")),s("omega")),
            {System::OriginKind::Model,"kOmegaSST.stageConservative"}}};
}
/// @brief 固定 SST 空间算子的真实输入/输出签名；并非任意 callback。
template<int Output> inline Math::Signature stageClosureSignature(const std::vector<Math::Signature>& args) {
    if(args.size()!=(Output==1?7:11))throw std::runtime_error("SST operator operand count mismatch");
    const std::array<Math::Signature,11> expected={mathematicalSignature("rho"),mathematicalSignature("U"),mathematicalSignature("k"),mathematicalSignature("omega"),mathematicalSignature("mu"),mathematicalSignature("omega"),Math::signature(0,{}),Math::signature(0,{}),Math::signature(0,{}),Math::signature(0,{}),Math::signature(0,{})};
    for(std::size_t i=0;i<(Output==1?7:10);++i)Math::requireCompatible(args[i],Output==1 && i>=5?Math::signature(0,{}):expected[i],"SST operand");
    if constexpr(Output!=1)Math::requireCompatible(args[10],mathematicalSignature("mu_t"),"published dynamic viscosity");
    const std::array<const char*,5> outputs={"SST.F1","mu_t","SST.production.k","SST.production.omega","SST.cross"};
    return mathematicalSignature(outputs[Output]);
}
inline const std::array<Math::OperatorRule,5>& sstRules() {
    static const std::array<Math::OperatorRule,5> rules={{{"sst.F1",stageClosureSignature<0>},{"sst.eddyViscosity",stageClosureSignature<1>},{"sst.kProduction",stageClosureSignature<2>},{"sst.omegaProduction",stageClosureSignature<3>},{"sst.crossDiffusion",stageClosureSignature<4>}}};
    return rules;
}
/// @brief 公式身份独立于模型活动性；实际字段仍须显式声明/由活动模型提供。
inline void declareClosureOperators(System::FormulaSymbols& symbols) {
    for(const auto& rule:sstRules())symbols.declareOperator(rule);
}
/// @brief 可复用独立黏度公式；输入顺序为 rho/U/k/omega/mu/betaStar/a1。
inline System::Equation sstViscosity(std::string id,System::FormulaExpr output,std::vector<System::FormulaExpr> inputs) {
    return {std::move(id),std::move(output),System::FormulaExpr::op(sstRules()[1],std::move(inputs)),{System::OriginKind::Model,"SST eddy dynamic viscosity"}};
}
/// @brief 原非局部 SST 数学组；梯度、壁距、配置与 regularization 属于明确所选实现。
inline std::vector<System::Equation> stageClosureMathematics(const std::string& density="rho",const std::string& velocity="",
        const System::FormulaExpr& viscosity=mathematicalSymbol("mu_t"),bool declareViscosity=true) {
    using E=System::FormulaExpr;
    const std::array<const char*,5> outputs={"SST.F1","mu_t","SST.production.k","SST.production.omega","SST.cross"};
    const auto& rules=sstRules();
    const auto rho=E::symbol(density,mathematicalSignature("rho"));
    const auto u=velocity.empty()?E::divide(mathematicalSymbol("rhoU"),rho):E::symbol(velocity,mathematicalSignature("U"));
    std::vector<System::Equation> result;
    for(std::size_t i=0;i<outputs.size();++i) {
        std::vector<E> inputs={rho,u,mathematicalSymbol("k"),mathematicalSymbol("omega"),mathematicalSymbol("kOmegaSST.laminarMu"),mathematicalSymbol("kOmegaSST.omegaFloor"),mathematicalSymbol("kOmegaSST.betaStar"),mathematicalSymbol("kOmegaSST.sigmaW2"),mathematicalSymbol("kOmegaSST.gamma1"),mathematicalSymbol("kOmegaSST.gamma2"),mathematicalSymbol("kOmegaSST.a1")};
        if(i==1) {
            if(declareViscosity)result.push_back(sstViscosity(viscosity.name(),viscosity,{inputs[0],inputs[1],inputs[2],inputs[3],inputs[4],inputs[6],inputs[10]}));continue;
        }
        inputs.pop_back(); // a1 belongs only to the independent viscosity relation.
        inputs.push_back(viscosity);
        result.emplace_back(outputs[i],mathematicalSymbol(outputs[i]),E::op(rules[i],std::move(inputs)),System::Provenance{System::OriginKind::Model,"kOmegaSST.nonlocal"});
    }
    return result;
}
/// @brief 动力黏度的局部公式，同一表达式可生成 WHAT AST 或求已绑定数值。
///        已有 floor 由选择该关系的调用入口提供。
template<class Value>
Value kEpsilonViscosity(const Value& cMu,const Value& rho,const Value& k,const Value& epsilon) {
    return cMu*rho*k*k/epsilon;
}
inline System::Equation stageTransportMathematics(const std::string& target) {
    using E=System::FormulaExpr;auto s=[](const std::string& n){return mathematicalSymbol(n);};
    if(target!="k" && target!="omega")throw std::runtime_error("Stage SST target unsupported.");
    auto q=E::multiply(s("rho"),s(target));
    auto rhs=E::add(E::op(Math::Op::Diffusion,{E::multiply(s("rho"),s("SST.nu."+target)),s(target)}),
        E::subtract(s("SST.production."+target),E::multiply(s("SST.sink."+target),s(target))));
    if(target=="omega")rhs=E::add(rhs,s("SST.cross"));
    return {target,E::add(E::op(Math::Op::Ddt,{q}),E::op(Math::Op::Div,{s("rhoU"),s(target)})),rhs,{System::OriginKind::Model,"kOmegaSST.stageConservative"}};
}
/// @brief 常密度规范化输运；WHAT保留div(U,a)，WHICH显式选择同一已发布phi。
inline System::Equation fixedVolumeMathematics(const std::string& target,const std::string& velocity) {
    using E=System::FormulaExpr;
    const auto a=mathematicalSymbol(target),rho=mathematicalSymbol("rho");
    const auto u=E::symbol(velocity,mathematicalSignature("U"));
    const auto production=target=="k"?mathematicalSymbol("SST.production.k")
        :E::add(mathematicalSymbol("SST.production.omega"),mathematicalSymbol("SST.cross"));
    return {target,E::op(Math::Op::Ddt,{a}),
        E::add(E::negate(E::op(Math::Op::Div,{u,a})),
            E::add(E::op(Math::Op::Diffusion,{mathematicalSymbol("SST.nu."+target),a},"diffusion."+target),
                E::subtract(E::divide(production,rho),E::multiply(E::divide(mathematicalSymbol("SST.sink."+target),rho),a)))),
        {System::OriginKind::Model,"kOmegaSST.fixedVolume"}};
}
inline System::Equation transportMathematics(const std::string& second,
                                            const std::string& target) {
    using E=System::FormulaExpr;
    const auto s=[](const std::string& name) { return mathematicalSymbol(name); };
    const auto mul=[](E a,E b) { return E::multiply(std::move(a),std::move(b)); };
    const auto div=[](E a,E b) { return E::divide(std::move(a),std::move(b)); };
    if ((second!="epsilon" && second!="omega") || (target!="k" && target!=second))
        throw std::runtime_error("Unsupported single-fluid RAS mathematical contract.");
    const bool epsilon=second=="epsilon";
    const std::string model=epsilon ? "kEpsilon" : "kOmegaSST";
    const auto coefficient=[&](const std::string& name) { return s(model+"."+name); };
    const auto production=coefficient("limitedProduction");
    E source;
    E diffusivity;
    if (target=="k") {
        const auto sink=epsilon ? s("epsilon")
            : mul(coefficient("betaStar"),mul(s("k"),s("omega")));
        source=E::subtract(div(production,s("rho")),sink);
    } else if (epsilon) {
        source=E::subtract(div(mul(coefficient("C1"),mul(s("epsilon"),production)),
                                  coefficient("maxRhoK")),
            div(mul(coefficient("C2"),mul(s("epsilon"),s("epsilon"))),coefficient("maxK")));
    } else {
        source=E::subtract(div(mul(coefficient("gamma"),production),coefficient("maxMuT")),
            mul(coefficient("beta"),mul(s("omega"),s("omega"))));
    }
    diffusivity=epsilon ? div(s("mu_t"),coefficient(target=="k" ? "sigmaK" : "sigmaEpsilon"))
        : E::add(coefficient("laminarMu"),mul(coefficient(target=="k" ? "sigmaK" : "sigmaOmega"),s("mu_t")));
    auto diffusion=E::op(Math::Op::Diffusion,{diffusivity,s(target)});
    if (!epsilon && target=="omega")
        diffusion=E::add(std::move(diffusion),coefficient("crossDiffusion"));
    return {target,E::op(Math::Op::Ddt,{s(target)}),
        E::add(std::move(source),div(std::move(diffusion),s("rho"))),
        {System::OriginKind::Model,model}};
}
} // namespace SF::Turbulence
