#include "SF_residual.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace SF::Run {
ResidualEvaluation::ResidualEvaluation(const System::SolvePlanNode& loop,const OpRegistry& operations,
        const ExecutionContext& execution):loop_(loop),operations_(operations),execution_(execution),metrics_(loop.residuals.size()) {
    for(std::size_t i=0;i<loop_.residuals.size();++i)if(loop_.residuals[i].kind==System::ResidualKind::StateChange
            && operations_.residualAvailable(loop_.residuals[i].id,loop_.residuals[i].provider,execution_)) {
        auto& previous=metrics_[i].previous;
        previous=operations_.sampleResidual(loop_.residuals[i].id,loop_.residuals[i].provider,execution_);
        if(loop_.residuals[i].component>=0) {
            const auto component=static_cast<std::size_t>(loop_.residuals[i].component);
            if(component>=previous.size())throw std::runtime_error("Missing initial residual component");
            previous={previous[component]};
        }
        if(previous.empty() || std::any_of(previous.begin(),previous.end(),[](double v){return !std::isfinite(v);}))
            throw std::runtime_error("Invalid initial residual state: "+loop_.residuals[i].id);
    }
}
bool ResidualEvaluation::sampleValid(std::size_t i,const Revisions& revisions) const {
    const auto at=revisions.find(loop_.residuals[i].producerOperation);
    return at!=revisions.end() && at->second>0
        && operations_.residualAvailable(loop_.residuals[i].id,loop_.residuals[i].provider,execution_);
}
bool ResidualEvaluation::evaluate(std::size_t i,int iteration) {
    const auto& definition=loop_.residuals[i];auto& metric=metrics_[i];
    auto sample=operations_.sampleResidual(definition.id,definition.provider,execution_);
    if(definition.component>=0) {
        if(static_cast<std::size_t>(definition.component)>=sample.size())throw std::runtime_error("Residual sample component is unavailable: "+definition.id);
        sample={sample[definition.component]};
    }
    if(definition.kind==System::ResidualKind::StateChange && metric.previous.empty()) {
        if(sample.empty() || std::any_of(sample.begin(),sample.end(),[](double v){return !std::isfinite(v);}))
            throw std::runtime_error("Invalid first residual state: "+definition.id);
        metric.previous=std::move(sample);
        std::cout<<"[SF RES] loop="<<loop_.id<<" metric="<<definition.id<<" state=pending-first-increment\n";
        return false;
    }
    if(sample.empty() || (definition.kind==System::ResidualKind::StateChange && sample.size()!=metric.previous.size()))
        throw std::runtime_error("Residual sampler changed layout or returned no values: "+definition.id);
    double norm=0;
    for(std::size_t j=0;j<sample.size();++j) {
        const double value=definition.kind==System::ResidualKind::StateChange?sample[j]-metric.previous[j]:sample[j];
        if(!std::isfinite(value) || !std::isfinite(sample[j]))throw std::runtime_error("Non-finite residual: "+definition.id);
        if(definition.norm==System::ResidualNorm::Maximum)norm=std::max(norm,std::abs(value));
        else norm=std::hypot(norm,value);
    }
    if(!metric.initialized){metric.initial=norm;metric.initialized=true;}
    if(definition.kind==System::ResidualKind::StateChange)metric.previous=std::move(sample);
    const double scale=definition.normalization==System::ResidualNormalization::Initial?metric.initial:definition.scale;
    const bool relative=definition.normalization!=System::ResidualNormalization::None
        && (scale==0?norm==0:norm/scale<=definition.relativeTolerance);
    const bool converged=norm<=definition.absoluteTolerance || relative;
    std::cout<<std::setprecision(17)<<"[SF RES] loop="<<loop_.id<<" iteration="<<iteration
        <<" metric="<<definition.id<<" occurrence="<<definition.occurrence<<" target="<<definition.target
        <<" view="<<System::targetText({definition.target,definition.view})
        <<" kind="<<System::toString(definition.kind)<<" normType="<<System::toString(definition.norm)<<" norm="<<norm
        <<" normalization="<<System::toString(definition.normalization)<<" reference="<<scale
        <<" absoluteTolerance="<<definition.absoluteTolerance<<" relativeTolerance="<<definition.relativeTolerance
        <<" accepted="<<converged<<'\n';
    return converged;
}
bool ResidualEvaluation::enter(const Revisions& revisions) {
    if(!loop_.checkAtEntry)return true;
    bool converged=loop_.residualAll;
    for(std::size_t i=0;i<loop_.residuals.size();++i) {
        // An increment has no initial pair. Other samplers are usable only
        // after their frozen producing operation has completed in this plan.
        if(!sampleValid(i,revisions) || loop_.residuals[i].kind==System::ResidualKind::StateChange) {
            std::cout<<"[SF RES] loop="<<loop_.id<<" metric="<<loop_.residuals[i].id<<" state=pending\n";
            return true;
        }
        const bool accepted=evaluate(i,0);
        converged=loop_.residualAll?converged && accepted:converged || accepted;
    }
    return !converged;
}
bool ResidualEvaluation::continueIteration(const Revisions& before,const Revisions& after,int iteration) {
    bool converged=loop_.residualAll;
    for(std::size_t i=0;i<loop_.residuals.size();++i) {
        const auto& definition=loop_.residuals[i];
        const auto old=before.find(definition.producerOperation),current=after.find(definition.producerOperation);
        if(current==after.end() || current->second==(old==before.end()?0:old->second))
            throw std::runtime_error("Residual '"+definition.id+"' was not updated in Loop '"+loop_.id+"'.");
        if(!operations_.residualAvailable(definition.id,definition.provider,execution_))
            throw std::runtime_error("Residual sample is stale/unpublished after its producing call: "+definition.id);
        const bool accepted=evaluate(i,iteration);
        converged=loop_.residualAll?converged && accepted:converged || accepted;
    }
    return !converged;
}
void ResidualEvaluation::finish(int iterations,bool converged) const {
    std::cout<<"[SF RES EXIT] loop="<<loop_.id<<" iterations="<<iterations
        <<" reason="<<(converged?"converged":loop_.failOnMaximum?"maximum-failed":"maximum-continue-explicit")<<'\n';
    if(!converged && loop_.failOnMaximum)
        throw std::runtime_error("Loop '"+loop_.id+"' reached maximum iterations without residual convergence.");
}
}
