#pragma once
#include "solver/algorithm/SF_singleFluidBinding.h"
#include "solver/algorithm/time/SF_stageGroup.h"
#include "solver/algorithm/SF_stageFlux.h"
#include <set>
namespace SF::SolverAlgorithm {
/// @brief 普通标量/向量共用的 occurrence 资源查找；不创建存储或选择方法。
const FDM::ScalarInstanceData& equationData(const FDM::SolverServices&,const std::string&);
/// @brief 普通数值调用共享的物理区间戳；不选择 dt、不推进 clock。
struct EquationInterval {
    int step=-1;double time=0,dt=0;
    void open(const State::StateBundle&,const std::string&);
    void require(const State::StateBundle&,const std::string&) const;
};
/// @brief 校验并借用普通 Current/Working 场；不创建数组或运行准备。
State::DistributedFieldView& equationField(System::StateRealization&,const System::StateUse&,
    Field&,bool writable=false,bool unsliced=false);
/// @brief 单元素显式recipe组合；输入R为正时间导数，不包含时间或存储所有权。
inline double explicitValue(const System::CompiledTimeRecipe& recipe,int n,double dt,double old,double current,
        double rhs,const std::array<double,4>& increments) {
    if(recipe.backend()==System::ExplicitStageBackend::SSPRK3)
        return recipe.stage(n).baseWeight*old+recipe.stage(n).incrementWeight*(current+dt*rhs);
    if(n+1<recipe.stageCount() || recipe.stageCount()==1)return old+recipe.stage(n).incrementWeight*dt*rhs;
    double sum=0;const auto weights=recipe.finalWeights();
    for(int i=0;i<recipe.stageCount();++i)sum+=weights[i]*increments[i];
    return old+dt*sum/recipe.finalDivisor();
}

/// @brief 绑定读者时可借用已实例化的共享通量。函数本身保留实际 owner。
struct FluxReaders {
    std::function<StageMassFluxView(StageFluxIdentity)> mass;
    std::function<double(std::string_view,int,int)> volume;
};

/// @brief 物理区间的启动动作与 dt 提议；带 prepare 的调用先建立状态和边界读取。
struct IntervalBinding {
    std::function<void()> prepare;
    std::function<double()> stepSize;
    std::function<void()> snapshot,accept;
    bool publishesOldTime=false;
};

/// @brief 一次 occurrence 的启动接线。构造时已挂接输出布局，数值操作随后安装。
struct BoundEquation {
    std::function<void(const FluxReaders&)> connect;
    std::function<void(Run::OpRegistry&)> install;
    Time::TemporalCallbacks stage;
    IntervalBinding interval;
    /// @brief 仅在实际借用 SolverServices::scalarEquation 时记录原端口。
    const FDM::ScalarEquationData* singleScalarPort=nullptr;
};

/// @brief 调用绑定阶段借用实际执行生命周期的数据；不能被运行回调按引用捕获。
struct EquationBindingInput {
    const System::CompiledEquationCall& call;
    const System::CompiledNumericalSystem& numerical;
    State::StateBundle& state;
    const FDM::SolverServices& services;
    System::StateRealization& views;
    const double& limit;
    std::vector<const System::ResidualDefinition*> residuals;
};

BoundEquation bindScalar(const EquationBindingInput&);
BoundEquation bindAffine(const EquationBindingInput&,const std::shared_ptr<std::set<std::string>>&);
BoundEquation bindSst(const EquationBindingInput&);
BoundEquation bindBlock(const EquationBindingInput&);
void bindVariableCommits(Run::OpRegistry&,const System::ExecutionScope&,State::StateBundle&,System::StateRealization&);
}
