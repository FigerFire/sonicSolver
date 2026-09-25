/// @file SF_field.cpp
/// @brief canonical 场数据、索引与状态存储实现。

#include "SF_field.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace SF {

    Field::Field(FieldStorageLayout layout)
        : mx(0), my(0), mz(0), ng(0), total_size(0),
          storageLayout_(layout) {}

    void Field::setStorageLayout(FieldStorageLayout layout) {
        if (total_size != 0) {
            throw std::runtime_error(
                "Field storage layout must be selected before setup; "
                "implicit state reordering is forbidden.");
        }
        storageLayout_ = layout;
    }

    double* Field::componentData(int variable) {
        return const_cast<double*>(
            static_cast<const Field&>(*this).componentData(variable));
    }

    const double* Field::componentData(int variable) const {
        if (storageLayout_ != FieldStorageLayout::SoA) {
            throw std::runtime_error(
                "Field::componentData requires SoA storage.");
        }
        if (variable < 0 || variable >= nVar_) {
            throw std::out_of_range(
                "Field::componentData variable is outside NVar.");
        }
        return data.data() + (size_t)variable * (size_t)total_size;
    }

    void Field::setup(int nx, int ny, int nz, int ghost, int nVar) {
        if (nVar == 0) nVar = stateModel_ ? stateModel_->variableCount() : 5;
        if (nVar < 5) {
            throw std::runtime_error("Field::setup requires at least five conserved variables.");
        }
        nVar_ = nVar;
        this->ng = ghost;
        this->mx = nx + 2 * ghost;
        this->my = ny + 2 * ghost;
        this->mz = nz + 2 * ghost;
        this->total_size = mx * my * mz;
        size_t total_vars = (size_t)total_size * (size_t)nVar_;

        data.assign(total_vars, 0.0);
        geometry_.x.assign(total_size, 0.0);
        geometry_.y.assign(total_size, 0.0);
        geometry_.z.assign(total_size, 0.0);
        boundary_.cellType.assign(total_size, FLUID_CELL);
        boundary_.solverMask.assign(total_size, 0);
        boundary_.communicationHaloMask.assign(total_size, 0);
        boundary_.globalDofId.assign(total_size, -1);
        boundary_.globalDofOwnerRank.assign(total_size, -1);
        boundary_.globalDofOwnerMask.assign(total_size, 1);
        ibmFluidMask_.assign(total_size, 1);
        geometry_.wallDistance.assign(total_size, 0.0);
        
        geometry_.jac.resize(total_size, 1.0);
        
        geometry_.xiX.resize(total_size, 0.0);
        geometry_.xiY.resize(total_size, 0.0);
        geometry_.xiZ.resize(total_size, 0.0);
        
        geometry_.etaX.resize(total_size, 0.0);
        geometry_.etaY.resize(total_size, 0.0);
        geometry_.etaZ.resize(total_size, 0.0);
        
        geometry_.zetaX.resize(total_size, 0.0);
        geometry_.zetaY.resize(total_size, 0.0);
        geometry_.zetaZ.resize(total_size, 0.0);
        geometry_.canonicalFaceMetrics.assign(
            (size_t)3 * total_size * 4, 0.0);
        geometry_.canonicalFaceMetricMask.assign(
            (size_t)3 * total_size, 0);
        geometry_.canonicalFaceOwnerMask.assign(
            (size_t)3 * total_size, 1);

        thermodynamicCache_.reset();
        
    }

    void Field::setStateModel(
            std::shared_ptr<const Physics::FluidStateModel::Model> equations) {
        if (!equations) throw std::runtime_error("Field::setStateModel requires a model.");
        if (total_size > 0 && nVar_ != equations->variableCount()) {
            throw std::runtime_error(
                "Field::setStateModel variable count differs from allocated state; "
                "call resizeConservedVariables explicitly before binding.");
        }
        stateModel_ = std::move(equations);
        thermodynamicCache_.bind(stateModel_);
        thermodynamicCache_.setCapacity(
            (int)std::min(total_size, 4096));
        thermodynamicCache_.invalidate();
    }

    void Field::resizeConservedVariables(int nVar) {
        if (nVar < 5) throw std::runtime_error("Field requires at least five conserved variables.");
        if (stateModel_ && nVar != stateModel_->variableCount()) {
            throw std::runtime_error("Field resize conflicts with bound FluidStateModel.");
        }
        nVar_ = nVar;
        data.assign((size_t)total_size * (size_t)nVar_, 0.0);
        thermodynamicCache_.reset();
        thermodynamicCache_.invalidate();
    }

    Physics::FluidStateModel::ThermodynamicState Field::thermodynamicState(
            int i, int j, int k) const {
        if (!stateModel_) {
            throw std::runtime_error("Field has no FluidStateModel thermodynamic closure.");
        }
        const int cell = getIdx(i, j, k);
        // cache hit 路径不得再收集守恒状态：只有 miss 才读取并闭合 EOS。
        if (const auto* cached = thermodynamicCache_.find(cell)) {
            return *cached;
        }
        // thread_local scratch 让普通 timestep 不再为每次查询分配 vector，
        // 同时保持按线程隔离（OpenMP 区间内并发查询的安全前提不变）。
        static thread_local std::vector<double> scratch;
        if (scratch.size() != (size_t)nVar_) scratch.resize((size_t)nVar_);
        for (int v = 0; v < nVar_; ++v) {
            scratch[(size_t)v] = (*this)(i, j, k, v);
        }
        return thermodynamicCache_.computeAndStore(cell, scratch.data(), nVar_);
    }

    void Field::invalidateThermodynamicCache() {
        thermodynamicCache_.invalidate();
    }

    void Field::setCanonicalFaceMetrics(
            int dir, int i, int j, int k,
            const std::array<double, 4>& metrics) {
        const size_t face = idxFace(dir, i, j, k);
        geometry_.canonicalFaceMetricMask[face] = 1;
        for (int c = 0; c < 4; ++c) {
            geometry_.canonicalFaceMetrics[face * 4 + (size_t)c] =
                metrics[(size_t)c];
        }
    }

    bool Field::hasCanonicalFaceMetrics(
            int dir, int i, int j, int k) const {
        return geometry_.canonicalFaceMetricMask[
            idxFace(dir, i, j, k)] != 0;
    }

    void Field::setCanonicalFaceOwner(
            int dir, int i, int j, int k, bool isOwner) {
        geometry_.canonicalFaceOwnerMask[idxFace(dir, i, j, k)] =
            isOwner ? 1 : 0;
    }

    bool Field::isCanonicalFaceOwner(
            int dir, int i, int j, int k) const {
        return geometry_.canonicalFaceOwnerMask[
            idxFace(dir, i, j, k)] != 0;
    }

    void Field::canonicalFaceMetrics(
            int dir, int i, int j, int k, double metrics[4]) const {
        const size_t face = idxFace(dir, i, j, k);
        if (!geometry_.canonicalFaceMetricMask[face]) {
            throw std::runtime_error(
                "canonicalFaceMetrics called for a non-canonical face.");
        }
        for (int c = 0; c < 4; ++c) {
            metrics[c] =
                geometry_.canonicalFaceMetrics[face * 4 + (size_t)c];
        }
    }

    void Field::clearCellFlags(int marker) {
        std::fill(
            boundary_.cellType.begin(), boundary_.cellType.end(), marker);
    }

    void Field::clearCommunicationHaloMask() {
        std::fill(boundary_.communicationHaloMask.begin(),
                  boundary_.communicationHaloMask.end(), 0);
    }

    void Field::getIJK(int idx, int& i, int& j, int& k) const {
        i = idx % mx;
        j = (idx / mx) % my;
        k = idx / (mx * my);
    }

    void Field::setBoundarySets(const std::map<std::string, std::vector<int>>& sets) {
        boundary_.sets = sets;
    }

    void Field::clearSolverBoundaryMask() {
        std::fill(
            boundary_.solverMask.begin(), boundary_.solverMask.end(), 0);
    }

    void Field::markSolverBoundarySet(const std::string& name) {
        auto it = boundary_.sets.find(name);
        if (it == boundary_.sets.end()) return;
        for (int idx : it->second) {
            if (idx >= 0 && idx < total_size) {
                boundary_.solverMask[(size_t)idx] = 1;
            }
        }
    }

    bool Field::isInSet(const std::string& name, int idx) const {
        auto it = boundary_.sets.find(name);
        if (it == boundary_.sets.end()) return false;
        const std::vector<int>& values = it->second;
        return std::find(values.begin(), values.end(), idx) != values.end();
    }

    void Field::markSolverBoundaryAxis(int axis) {
        const int i0 = ng;
        const int i1 = ng + NX() - 1;
        const int j0 = ng;
        const int j1 = ng + NY() - 1;
        const int k0 = ng;
        const int k1 = ng + NZ() - 1;

        auto mark = [&](int i, int j, int k) {
            int idx = getIdx(i, j, k);
            if (idx >= 0 && idx < total_size) {
                boundary_.solverMask[(size_t)idx] = 1;
            }
        };

        if (axis == 0) {
            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j) {
                    mark(i0, j, k);
                    if (i1 != i0) mark(i1, j, k);
                }
        } else if (axis == 1) {
            for (int k = k0; k <= k1; ++k)
                for (int i = i0; i <= i1; ++i) {
                    mark(i, j0, k);
                    if (j1 != j0) mark(i, j1, k);
                }
        } else if (axis == 2) {
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    mark(i, j, k0);
                    if (k1 != k0) mark(i, j, k1);
                }
        }
    }

    bool Field::isSolverBoundaryPoint(int i, int j, int k) const {
        int idx = getIdx(i, j, k);
        if (idx < 0 || idx >= total_size) return false;
        return boundary_.solverMask[(size_t)idx] != 0;
    }

    const std::vector<int>& Field::getSet(const std::string& name) const {
        // 使用 find 而不是 at，手动处理报错
        auto it = boundary_.sets.find(name);
        
        if (it == boundary_.sets.end()) {
            std::cerr << "\n" << std::string(50, '*') << std::endl;
            std::cerr << "[SF FATAL ERROR] Map Key Not Found!" << std::endl;
            std::cerr << "[SF DEBUG] The code is looking for: '" << name << "'" << std::endl;
            std::cerr << "[SF DEBUG] But the Mesh file only contains: ";
            for (auto const& [key, val] : boundary_.sets) {
                std::cerr << "'" << key << "' ";
            }
            std::cerr << "\n" << std::string(50, '*') << std::endl;
            
            // 抛出异常或退出
            throw std::out_of_range("Key missing: " + name);
        }
        
        return it->second;
    }

}
