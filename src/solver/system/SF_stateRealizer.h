#pragma once

/// @file SF_stateRealizer.h
/// @brief STATE runtime: bind physical authorities and realize requested views.
#include "core/system/SF_stateRegistry.h"
#include "core/system/SF_stateViews.h"
#include "SF_runtimeRequirements.h"
#include "core/state/SF_stateBundle.h"
#include <memory>
#include <vector>

namespace SF::System {
/// @brief A stable non-owning handle to a base variable.
struct RealizedUnknown {
    const StateSymbol* descriptor=nullptr;
    std::vector<State::DistributedFieldView*> fields;
};
/// @brief A compiler-requested view, independent of equation identity or order.
struct RealizedStateView {
    CompiledStateView descriptor;
    std::vector<State::DistributedFieldView*> fields;
};
/// @brief Solver execution owns temporary storage; StateBundle remains the physical authority.
class StateRealization {
public:
    const RealizedUnknown& at(std::string_view id) const;
    const RealizedStateView& view(std::string_view symbol,StateViewKind kind,
        int stage=-1) const;
    bool requestsView(std::string_view symbol,StateViewKind kind,int stage=-1) const;
    /// @brief Bind a numerical provider's existing array, without allocating a second copy.
    void bindView(std::string_view symbol,StateViewKind kind,
        State::DistributedFieldView& field,int stage=-1);
    std::size_t size() const { return unknowns_.size(); }
    /// @brief Scope-qualified reads; binding is stable, availability is per stage.
    std::function<double(int,int)> stageReader(std::string_view symbol,int stages);
    void beginStage(int index,int count,double time);
    void stageReady();
    void rhsReady();
    void finishStage();
    void requireStageRead(int index,double time) const;
    void requireStageAdvance(int index,double time) const;
    /// @brief provider 已绑定输入版本后，在启动阶段绑定显式声明的派生视图。
    void bindDerivedViews();
    /// @brief 仅供已有 HOW Update callback 调用；读取绝不隐式更新。
    void updateView(std::string_view symbol,StateViewKind kind,int stage=-1);
private:
    enum class StagePhase { Closed, Preparing, Reading, Advancing };
    struct StageEpoch { int index=-1;double time=0;StagePhase phase=StagePhase::Closed; };
    std::shared_ptr<StageEpoch> epoch_=std::make_shared<StageEpoch>();
    State::StateBundle* state_=nullptr;
    const StateRegistry* symbols_=nullptr;
    struct DerivedView {
        CompiledStateView descriptor;
        State::DistributedFieldView field;
        std::function<void()> update;
        std::vector<unsigned long long> versions;
        unsigned long long outputVersion=0;
        bool ready=false;
    };
    std::vector<std::shared_ptr<DerivedView>> derived_;
    void bindDerived(const CompiledStateView& descriptor);
    friend StateRealization realizeState(const StateRegistry&,const RuntimeRequirements&,
        State::StateBundle&,const std::vector<CompiledStateView>&);
    struct ViewStorage {
        std::vector<double> values;
        State::DistributedFieldView field;
    };
    std::vector<RealizedUnknown> unknowns_;
    std::vector<RealizedStateView> views_;
    std::vector<std::unique_ptr<ViewStorage>> owned_;
};

/// @brief Startup-only allocation: absent Working/Correction demands allocate nothing.
StateRealization realizeState(const StateRegistry& symbols,
    const RuntimeRequirements& requirements,State::StateBundle& state,
    const std::vector<CompiledStateView>& views={});
} // namespace SF::System
