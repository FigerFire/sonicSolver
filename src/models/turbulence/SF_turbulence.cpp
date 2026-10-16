#include <set>
/// @file SF_turbulence.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_turbulence.h"

#include "SF_DNS.h"
#include "SF_LES.h"
#include "SF_RAS.h"
#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "core/mesh/SF_dimension.h"
#include "SF_kEpsilon.h"
#include "SF_kOmegaSST.h"
#include "SF_viscous.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>

namespace SF {
namespace Turbulence {

namespace {

std::unique_ptr<IModel> makeModel(const FDM::TurbulenceConfig& config,const FDM::IImmersedTurbulenceBoundary* wall) {
    if (!config.enabled) return nullptr;

    switch (config.family) {
        case FDM::TurbulenceFamily::RAS:
            if (config.model == FDM::TurbulenceModelKind::kOmegaSST) {
                return std::make_unique<RAS::KOmegaSSTModel>(wall);
            }
            return std::make_unique<RAS::KEpsilonModel>();

        case FDM::TurbulenceFamily::LES:
            return std::make_unique<LES::SmagorinskyModel>();

        case FDM::TurbulenceFamily::DNS:
            return std::make_unique<DNS::DirectNumericalSimulationModel>();

        case FDM::TurbulenceFamily::None:
            break;
    }

    return nullptr;
}

int clampToInterior(int index, int lower, int upper) {
    return std::max(lower, std::min(index, upper));
}


} // namespace

void ScalarFields::resizeLike(const Field& field) {
    mx_ = field.MX();
    my_ = field.MY();
    mz_ = field.MZ();
    ng_ = field.NG();
    totalSize_ = mx_ * my_ * mz_;

    k_.assign((size_t)totalSize_, 0.0);
    epsilon_.assign((size_t)totalSize_, 0.0);
    omega_.assign((size_t)totalSize_, 0.0);
    eddyMu_.assign((size_t)totalSize_, 0.0);
}

int ScalarFields::idx(int i, int j, int k) const {
    return (k * my_ + j) * mx_ + i;
}

double& ScalarFields::K(int i, int j, int k) { return k_[(size_t)idx(i, j, k)]; }
double ScalarFields::K(int i, int j, int k) const { return k_[(size_t)idx(i, j, k)]; }
double& ScalarFields::Epsilon(int i, int j, int k) { return epsilon_[(size_t)idx(i, j, k)]; }
double ScalarFields::Epsilon(int i, int j, int k) const { return epsilon_[(size_t)idx(i, j, k)]; }
double& ScalarFields::Omega(int i, int j, int k) { return omega_[(size_t)idx(i, j, k)]; }
double ScalarFields::Omega(int i, int j, int k) const { return omega_[(size_t)idx(i, j, k)]; }
double& ScalarFields::EddyMu(int i, int j, int k) { return eddyMu_[(size_t)idx(i, j, k)]; }
double ScalarFields::EddyMu(int i, int j, int k) const { return eddyMu_[(size_t)idx(i, j, k)]; }

double ScalarFields::clamped(ScalarSlot slot, int i, int j, int k) const {
    int ii = clampToInterior(i, ng_, mx_ - ng_ - 1);
    int jj = clampToInterior(j, ng_, my_ - ng_ - 1);
    int kk = clampToInterior(k, ng_, mz_ - ng_ - 1);
    return this->slot(slot, ii, jj, kk);
}

double& ScalarFields::slot(ScalarSlot slot, int i, int j, int k) {
    const size_t flatIndex = (size_t)idx(i, j, k);
    switch (slot) {
        case ScalarSlot::K: return k_[flatIndex];
        case ScalarSlot::Epsilon: return epsilon_[flatIndex];
        case ScalarSlot::Omega: return omega_[flatIndex];
        case ScalarSlot::EddyMu: return eddyMu_[flatIndex];
    }
    return k_[flatIndex];
}

double ScalarFields::slot(ScalarSlot slot, int i, int j, int k) const {
    const size_t flatIndex = (size_t)idx(i, j, k);
    switch (slot) {
        case ScalarSlot::K: return k_[flatIndex];
        case ScalarSlot::Epsilon: return epsilon_[flatIndex];
        case ScalarSlot::Omega: return omega_[flatIndex];
        case ScalarSlot::EddyMu: return eddyMu_[flatIndex];
    }
    return k_[flatIndex];
}

std::vector<double>& ScalarFields::values(ScalarSlot slot) {
    switch (slot) {
        case ScalarSlot::K: return k_;
        case ScalarSlot::Epsilon: return epsilon_;
        case ScalarSlot::Omega: return omega_;
        case ScalarSlot::EddyMu: return eddyMu_;
    }
    throw std::runtime_error("Unknown turbulence scalar slot.");
}

const std::vector<double>& ScalarFields::values(ScalarSlot slot) const {
    switch (slot) {
        case ScalarSlot::K: return k_;
        case ScalarSlot::Epsilon: return epsilon_;
        case ScalarSlot::Omega: return omega_;
        case ScalarSlot::EddyMu: return eddyMu_;
    }
    throw std::runtime_error("Unknown turbulence scalar slot.");
}

Manager::Manager(FDM::TurbulenceConfig config,const FDM::IImmersedTurbulenceBoundary* wall)
    : wall_(wall), config_(std::move(config)) {
    if (wall_ && (config_.model!=FDM::TurbulenceModelKind::kOmegaSST || config_.laminarDynamicViscosity<=0))
        throw std::runtime_error("Bound immersed turbulence wall implements only SST with positive molecular viscosity.");
}

bool Manager::initialize(const Field& field) {
    model_ = makeModel(config_,wall_);
    if (!model_) return false;

    state_.resizeLike(field);
    model_->initialize(field, state_, config_);
    return true;
}

bool Manager::active() const {
    return model_ != nullptr;
}

std::string Manager::description() const {
    if (!model_) return "None";
    return std::string(model_->familyName()) + "/" + model_->modelName();
}

void Manager::applyBoundary(const Field& field) {
    if (model_) model_->applyPhysicalBoundary(field, state_, config_);
}

void Manager::refreshClosure(const Field& field) {
    if (model_) model_->refreshClosure(field, state_, config_);
}
FDM::OnceStepTransportData Manager::onceStepTransport() {
    FDM::OnceStepTransportData data;
    if (!model_ || config_.family!=FDM::TurbulenceFamily::RAS) return data;
    data.variables={&state_.values(ScalarSlot::K),&state_.values(config_.model==FDM::TurbulenceModelKind::kEpsilon ? ScalarSlot::Epsilon : ScalarSlot::Omega)};
    data.eddyViscosity=&state_.values(ScalarSlot::EddyMu);
    data.boundary={&config_.scalars.kBoundary,config_.model==FDM::TurbulenceModelKind::kEpsilon ? &config_.scalars.epsilonBoundary : &config_.scalars.omegaBoundary};
    data.floors={config_.coefficients.kFloor,config_.model==FDM::TurbulenceModelKind::kEpsilon ? config_.coefficients.epsilonFloor : config_.coefficients.omegaFloor};
    data.molecularDiffusion=config_.model==FDM::TurbulenceModelKind::kOmegaSST;
    data.molecularViscosity=config_.laminarDynamicViscosity;
    data.cellClosure=[this](const Field& flow,int i,int j,int k) { return model_->onceStepClosure(flow,state_,config_,i,j,k); };
    return data;
}

double Manager::dynamicViscosity(const Field& field,
                                 int i, int j, int k,
                                 double laminarMu) const {
    if (!model_) return laminarMu;
    int ii = clampToInterior(i, field.NG(), field.NG() + field.NX() - 1);
    int jj = clampToInterior(j, field.NG(), field.NG() + field.NY() - 1);
    int kk = clampToInterior(k, field.NG(), field.NG() + field.NZ() - 1);
    if (wall_ || stageClosurePublished_) return laminarMu+state_.EddyMu(ii,jj,kk);
    return laminarMu + std::max(model_->eddyDynamicViscosity(field, state_, config_, ii, jj, kk), 0.0);
}

std::vector<std::string> Manager::stageScalarTargets() const {
    if(model_ && config_.model==FDM::TurbulenceModelKind::kOmegaSST && !wall_)return {"k","omega"};
    return {};
}
FDM::TransportStageCoefficients Manager::stageCoefficients(const Field& field,
        const FDM::TransportStageSamples& samples,int i,int j,int k) const {
    if(stageScalarTargets().empty())throw std::runtime_error("Unsupported stage turbulence closure/model/immersed wall.");
    return model_->stageCoefficients(field,samples,config_,i,j,k);
}
void Manager::publishStageViscosity(const Field& field,const std::vector<double>& values) {
    if(stageScalarTargets().empty() || values.size()!=static_cast<std::size_t>(field.TotalSize()))
        throw std::runtime_error("Stage viscosity publication model/layout mismatch.");
    for(double v:values)if(!std::isfinite(v) || v<0)throw std::runtime_error("Stage eddy viscosity invalid.");
    std::copy(values.begin(),values.end(),state_.values(ScalarSlot::EddyMu).begin());stageClosurePublished_=true;
}
FDM::ScalarEquationData Manager::stageScalarData(const Field& field,const std::string& name) const {
    if(stageScalarTargets().empty() || (name!="k" && name!="omega"))throw std::runtime_error("Unsupported stage scalar boundary target.");
    FDM::ScalarEquationData data;data.target=name;
    const auto& settings=name=="k"?config_.scalars.kBoundary:config_.scalars.omegaBoundary;
    if(settings.empty())throw std::runtime_error("Native stage turbulence requires explicit scalar boundary settings.");
    const int g=field.NG(),shape[3]={field.NX(),field.NY(),field.NZ()};
    for(const auto& setting:settings) {
        if(setting.type==BCType::EMPTY)continue;
        if(setting.type!=BCType::FIXED_VALUE && setting.type!=BCType::ZERO_GRADIENT)
            throw std::runtime_error("Stage turbulence boundary supports fixedValue/zeroGradient only.");
        auto found=field.getAllSets().find(setting.name);
        if(found==field.getAllSets().end())throw std::runtime_error("Stage turbulence boundary set absent: "+setting.name);
        std::set<int> covered;
        for(int axis=0;axis<3;++axis)if(shape[axis]>1 && Math::isDirectionActiveIndex(axis))for(int side=0;side<2;++side) {
            std::set<int> face;
            for(const auto& point:found->second) {
                int i,j,k;field.getIJK(point,i,j,k);const int at[3]={i,j,k};
                bool physical=true;for(int a=0;a<3;++a)physical=physical && at[a]>=g && at[a]<g+shape[a];
                if(!physical)continue; // Mesh tags include inherited halo replicas; only physical face nodes define the law.
                if(at[axis]==g+(side?shape[axis]-1:0))face.insert(point);
            }
            const std::size_t expected=static_cast<std::size_t>(shape[(axis+1)%3])*shape[(axis+2)%3];
            if(face.size()!=expected)continue;
            auto& bc=data.boundary[2*axis+side];
            if(bc.kind!=FDM::ScalarBoundaryKind::Unspecified)throw std::runtime_error("Overlapping native turbulence boundary face.");
            bc.kind=setting.type==BCType::FIXED_VALUE?FDM::ScalarBoundaryKind::FixedValue:FDM::ScalarBoundaryKind::ZeroGradient;
            if(setting.type==BCType::FIXED_VALUE)bc.value=[v=setting.value](double,double,double,double){return v;};
            covered.insert(face.begin(),face.end());
        }
        std::set<int> physicalSet;
        for(int point:found->second){int i,j,k;field.getIJK(point,i,j,k);if(i>=g && i<g+shape[0] && j>=g && j<g+shape[1] && k>=g && k<g+shape[2])physicalSet.insert(point);}
        if(covered.empty() || covered!=physicalSet)
            throw std::runtime_error("Native stage turbulence boundary requires complete axis-aligned face(s): "+setting.name);
    }
    for(int axis=0;axis<3;++axis)if(shape[axis]>1 && Math::isDirectionActiveIndex(axis))for(int side=0;side<2;++side)
        if(data.boundary[2*axis+side].kind==FDM::ScalarBoundaryKind::Unspecified)throw std::runtime_error("Missing native stage turbulence boundary face.");
    return data;
}

std::vector<std::string> Manager::distributedReadFields() const {
    if (!model_ || config_.family != FDM::TurbulenceFamily::RAS) return {};
    if (config_.model == FDM::TurbulenceModelKind::kEpsilon) {
        return {"k", "epsilon", "mu_t"};
    }
    if (config_.model == FDM::TurbulenceModelKind::kOmegaSST) {
        return {"k", "omega", "mu_t"};
    }
    return {};
}

std::vector<std::string> Manager::distributedWriteFields() const {
    if (!model_) return {};
    auto result = distributedReadFields();
    if (result.empty()) result.push_back("mu_t");
    return result;
}

int Manager::distributedHaloDepth() const {
    return distributedReadFields().empty() ? 0 : 1;
}

void Manager::registerDistributed(State::DistributedFieldRegistry& registry,
                                 int blockId, Field& geometry) {
    for (const auto& name:distributedWriteFields()) {
        const auto slot=name=="k" ? ScalarSlot::K
            : name=="epsilon" ? ScalarSlot::Epsilon
            : name=="omega" ? ScalarSlot::Omega : ScalarSlot::EddyMu;
        registry.add(State::workspaceView(name,blockId,geometry,state_.values(slot),1));
    }
}

void applyScalarInitialConditions(const Field& flow,
                                  const std::vector<BCSetting<double>>& settings,
                                  ScalarFields& state,
                                  ScalarSlot slot) {
    for (const auto& ic : settings) {
        const auto& indices = flow.getSet(ic.name);
        for (int flatIndex : indices) {
            int i, j, k;
            flow.getIJK(flatIndex, i, j, k);
            state.slot(slot, i, j, k) = ic.value;
        }
    }
}


double clampPositive(double value, double floorValue) {
    return std::max(value, floorValue);
}

double strainRateMagnitude(const Field& flow, int i, int j, int k) {
    using Viscous::CentralOrder;

    // Use physical gradients (chain-rule with metric terms) instead of
    // raw computational-space derivatives, so the strain-rate tensor is
    // correct on stretched / curvilinear grids.
    const SymmTensor3 strain = Math::symm(Viscous::velocityGradientAt(
        flow, i, j, k, CentralOrder::SECOND));
    return std::sqrt(std::max(
        2.0 * Math::doubleDot(strain, strain), 0.0));
}

double characteristicFilterWidth(const Field& flow,
                                 int i, int j, int k,
                                 double filterScale) {
    double jac = std::abs(flow.Jac(i, j, k));
    double base = std::cbrt(1.0 / std::max(jac, 1.0e-12));
    return std::max(filterScale, 1.0e-12) * base;
}


} // namespace Turbulence
} // namespace SF
