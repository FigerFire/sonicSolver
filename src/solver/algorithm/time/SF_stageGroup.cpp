#include "SF_stageGroup.h"
#include "solver/system/SF_methodObjects.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
namespace SF::Time {
namespace {
struct Epoch {
    State::StateBundle& state;
    System::StateRealization& views;
    const System::CompiledTimeRecipe& recipe;
    std::vector<TemporalCallbacks> callbacks;
    std::vector<int> status;
    int next=0, snapshots=0, publications=0;
    bool active=false, publicationReady=false;
    double time() const {return state.time+recipe.stage(next).abscissa*state.dt;}
    void context(const Run::ExecutionContext& c) const {
        if (!active || c.stageIndex!=next || c.stageCount!=recipe.stageCount())
            throw std::runtime_error("Temporal participant stage index/recipe mismatch.");
    }
    void all(int value) const {
        if (!std::all_of(status.begin(),status.end(),[&](int x){return x==value;}))
            throw std::runtime_error("Temporal barrier: a participant phase is missing or incomplete.");
    }
};
}
void bindStageGroup(Run::OpRegistry& ops,const System::CompiledSolvePlan& plan,
        const System::CompiledNumericalSystem& numerical,State::StateBundle& state,
        System::StateRealization& views,const double& limit,std::vector<TemporalCallbacks> callbacks) {
    const auto reason=System::validateTemporalPlan(plan,numerical.time.recipe.stageCount());
    if(!reason.empty())throw std::runtime_error(reason);
    const auto& participants=plan.compiledProgram.temporalParticipants;
    if (callbacks.size()!=participants.size() || callbacks.empty())
        throw std::runtime_error("Temporal participant callbacks incomplete.");
    for (std::size_t i=0;i<callbacks.size();++i)
        if (callbacks[i].identity!=participants[i].identity || !callbacks[i].snapshot
            || !callbacks[i].prepare || !callbacks[i].rhs || !callbacks[i].advance
            || !callbacks[i].stepSize || !callbacks[i].validatePublish || !callbacks[i].publish)
            throw std::runtime_error("Missing/mismatched temporal participant binding.");
    auto epoch=std::make_shared<Epoch>(Epoch{state,views,numerical.time.recipe,std::move(callbacks),
        std::vector<int>(participants.size())});
    ops.bind(System::TemporalOps::Dt,System::TemporalOps::Provider,[epoch,&limit,&numerical] {
        if (epoch->active) throw std::runtime_error("Cannot restart an incomplete temporal step.");
        double dt=std::min(limit,numerical.dt.maxDeltaT);
        for (const auto& p:epoch->callbacks) {
            const double proposal=p.stepSize();
            if (!std::isfinite(proposal) || proposal<=0) throw std::runtime_error("Invalid participant dt proposal.");
            dt=std::min(dt,proposal);
        }
        if (!std::isfinite(dt) || dt<=0) throw std::runtime_error("Invalid common temporal dt.");
        epoch->state.dt=dt;epoch->active=true;epoch->next=0;
        epoch->snapshots=epoch->publications=0;epoch->publicationReady=false;
        std::fill(epoch->status.begin(),epoch->status.end(),0);
    });
    for (std::size_t i=0;i<participants.size();++i) {
        const auto& p=participants[i];
        ops.bind(p.snapshot,p.provider,[epoch,i] {
            if (!epoch->active || epoch->status[i]!=0 || epoch->next!=0)
                throw std::runtime_error("Duplicate/out-of-order temporal snapshot.");
            epoch->callbacks[i].snapshot();epoch->status[i]=1;++epoch->snapshots;
        });
        ops.bind(p.prepareStage,p.provider,[epoch,i](const Run::ExecutionContext& c) {
            epoch->context(c);
            if (epoch->status[i]!=0) throw std::runtime_error("Duplicate stage preparation.");
            epoch->callbacks[i].prepare(c);epoch->status[i]=1;
        });
        ops.bind(p.rhs,p.provider,[epoch,i](const Run::ExecutionContext& c) {
            epoch->context(c);epoch->views.requireStageRead(c.stageIndex,epoch->time());
            if (epoch->status[i]!=1) throw std::runtime_error("Duplicate/missing participant RHS.");
            epoch->callbacks[i].rhs(c);epoch->status[i]=2;
        });
        ops.bind(p.advance,p.provider,[epoch,i](const Run::ExecutionContext& c) {
            epoch->context(c);epoch->views.requireStageAdvance(c.stageIndex,epoch->time());
            if (epoch->status[i]!=2) throw std::runtime_error("Advance before participant RHS or duplicate advance.");
            epoch->callbacks[i].advance(c);epoch->status[i]=3;
        });
        ops.bind(p.publish,p.provider,[epoch,i] {
            if (!epoch->active || !epoch->publicationReady || epoch->status[i]!=3)
                throw std::runtime_error("Early/duplicate physical publication.");
            epoch->callbacks[i].publish();epoch->status[i]=4;
            if (++epoch->publications==(int)epoch->callbacks.size()) epoch->active=false;
        });
    }
    ops.bind(System::TemporalOps::Open,System::TemporalOps::Provider,[epoch](const Run::ExecutionContext& c) {
        epoch->context(c);
        if (epoch->snapshots!=(int)epoch->callbacks.size()) throw std::runtime_error("Stage before all snapshots.");
        epoch->all(epoch->next==0?1:3);
        epoch->views.beginStage(c.stageIndex,c.stageCount,epoch->time());
        std::fill(epoch->status.begin(),epoch->status.end(),0);
    });
    ops.bind(System::TemporalOps::Ready,System::TemporalOps::Provider,[epoch] {epoch->all(1);epoch->views.stageReady();});
    ops.bind(System::TemporalOps::RhsReady,System::TemporalOps::Provider,[epoch] {epoch->all(2);epoch->views.rhsReady();});
    ops.bind(System::TemporalOps::Close,System::TemporalOps::Provider,[epoch] {epoch->all(3);epoch->views.finishStage();++epoch->next;});
    ops.bind(System::TemporalOps::PublishReady,System::TemporalOps::Provider,[epoch] {
        if (!epoch->active || epoch->publicationReady || epoch->next!=epoch->recipe.stageCount())
            throw std::runtime_error("Physical publication before complete recipe or repeated commit.");
        epoch->all(3);
        for (const auto& callback:epoch->callbacks) callback.validatePublish();
        epoch->publicationReady=true;
    });
}
}
