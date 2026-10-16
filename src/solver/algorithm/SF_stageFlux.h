#pragma once

/// @file SF_stageFlux.h
/// @brief 借用 solver workspace 的 canonical stage 质量通量；无独立数值数组。
#include "core/field/SF_field.h"
#include "core/flux/SF_flux.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace SF::SolverAlgorithm {
/// @brief 物理步、RK stage 和时间共同标识一次数值通量生产。
struct StageFluxIdentity {
    int step=0,stage=0;
    double time=0;
    bool operator==(const StageFluxIdentity& b) const {
        return step==b.step && stage==b.stage && time==b.time;
    }
};
/// @brief 生命周期 token 仅保存元数据；canonical 数值仍唯一存于 PatchWorkspace。
struct StageFluxPublication {
    const Field* field=nullptr;
    const FluxField* flux=nullptr;
    const Physics::FluidStateModel::Model* model=nullptr;
    StageFluxIdentity identity;
    std::string recipe;
    int massComponent=-1,total=0,nvar=0;
    unsigned long generation=0;
    bool canonicalReady=false,finalized=false;
    void invalidate() { canonicalReady=finalized=false;++generation; }
    void finalizeCanonical() {if(canonicalReady)throw std::runtime_error("Duplicate canonical finalize.");canonicalReady=true;}
    void publish(const Field& f,const FluxField& storage,StageFluxIdentity stamp,std::string selectedRecipe) {
        if(!canonicalReady)throw std::runtime_error("Canonical mass flux publication before finalize/COPY barrier.");
        if (finalized) throw std::runtime_error("Duplicate canonical mass flux publication in one stage.");
        const auto activeModel=f.stateModel();
        if (!activeModel || activeModel->densityVariableCount()!=1 || activeModel->variableCount()!=f.NVar())
            throw std::runtime_error("Canonical mass flux requires one density in the active FluidStateModel.");
        int component=-1;
        for (int i=0;i<activeModel->variableCount();++i) if(activeModel->primaryVariables()[i].name=="rho") {
            if(component!=-1)throw std::runtime_error("Ambiguous FluidStateModel density component.");
            component=i;
        }
        if (component<0 || storage.variableCount()!=f.NVar()
            || storage.faceCount()!=3*static_cast<std::size_t>(f.TotalSize()) || selectedRecipe.empty())
            throw std::runtime_error("Canonical mass flux component/face layout/recipe mismatch.");
        if(stamp.step<0 || stamp.stage<0 || !std::isfinite(stamp.time))throw std::runtime_error("Canonical mass flux physical step/stage/time invalid.");
        field=&f;flux=&storage;model=activeModel.get();identity=stamp;recipe=std::move(selectedRecipe);
        massComponent=component;total=f.TotalSize();nvar=f.NVar();finalized=true;
    }
};
/// @brief 正计算方向面 cofactor 加权通量 F_hat；散度为 Jac*(F_right-F_left)。
/// 不额外乘面积或除物理 dx。只接受有效同 stage 的已完成 canonical COPY 数据。
class StageMassFluxView {
    std::weak_ptr<const StageFluxPublication> lifetime_;
    unsigned long generation_=0;
    StageFluxIdentity identity_;
    const Field* field_=nullptr;
    std::string recipe_;
public:
    StageMassFluxView(const std::shared_ptr<StageFluxPublication>& source,const Field& field,
                      StageFluxIdentity stamp,const std::string& recipe)
        :lifetime_(source),generation_(source?source->generation:0),identity_(stamp),field_(&field),recipe_(recipe) {
        validate();
    }
    std::shared_ptr<const StageFluxPublication> validate() const {
        auto s=lifetime_.lock();
        if(!s || !s->finalized || s->generation!=generation_ || !(s->identity==identity_))
            throw std::runtime_error("Canonical mass flux unavailable/stale: consumer requires finalized same physical step/stage epoch.");
        if(s->field!=field_ || s->model!=field_->stateModel().get() || s->total!=field_->TotalSize() || s->nvar!=field_->NVar()
            || s->flux->faceCount()!=3*static_cast<std::size_t>(s->total)
            || s->flux->variableCount()!=s->nvar || s->recipe!=recipe_)
            throw std::runtime_error("Canonical mass flux consumer patch/layout/provider recipe mismatch.");
        return s;
    }
    double at(int direction,int i,int j,int k) const {
        auto s=validate();
        if(direction<0 || direction>2 || i<0 || j<0 || k<0 || i>=field_->MX() || j>=field_->MY() || k>=field_->MZ())
            throw std::runtime_error("Canonical mass flux face orientation/index mismatch.");
        const int g=field_->NG();
        const int at[3]={i,j,k},shape[3]={field_->NX(),field_->NY(),field_->NZ()};
        for(int axis=0;axis<3;++axis) if(at[axis]<g-(axis==direction?1:0) || at[axis]>=g+shape[axis])
            throw std::runtime_error("Canonical mass flux read outside participating face layout.");
        const double value=(*s->flux)(static_cast<std::size_t>(direction)*s->total+field_->getIdx(i,j,k),s->massComponent);
        if(!std::isfinite(value))throw std::runtime_error("Non-finite canonical mass flux.");
        return value;
    }
};
}
