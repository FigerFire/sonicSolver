/// @file SF_kktDofs.cpp
/// @brief 构造互不碰撞的 pressure/velocity/constraint/solid DOF 空间。

#include "SF_kktDofs.h"

#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace SF::PressureBased {

KKTSubspace activeKKTSubspace(const FDM::ImmersedSolidEquation& solid) {
    KKTSubspace result;
    const bool solveGeneralizedVelocity=solid.solveGeneralizedVelocity;
    if(solveGeneralizedVelocity && solid.activeComponents.empty())
        throw std::runtime_error("KKT solid equation has no declared active components.");
    const int dimensions=Math::activeDimensionCount();
    if(dimensions==3) {
        result.velocityComponents={0,1,2};
        if(solveGeneralizedVelocity)result.solidComponents=solid.activeComponents;
        return result;
    }
    if(dimensions!=2) {
        throw std::runtime_error(
            "surface KKT supports exactly two or three active dimensions.");
    }
    const auto normal=Math::inactiveDirectionNormalVector();
    int normalComponent=-1;
    bool oblique=false;
    for(int component=0;component<3;++component) {
        if(std::abs(normal[(size_t)component])>1.0-1.0e-12) {
            if(normalComponent>=0)oblique=true;
            else normalComponent=component;
        } else if(std::abs(normal[(size_t)component])>1.0e-12) {
            oblique=true;
        }
    }
    if(normalComponent<0 || oblique) {
        throw std::runtime_error(
            "surface KKT planar reduction currently requires the empty "
            "physical normal to align with a Cartesian component.");
    }
    for(int component=0;component<3;++component)
        if(component!=normalComponent)result.velocityComponents.push_back(component);
    if(solveGeneralizedVelocity) {
        for(int c:solid.activeComponents)
            if((c<3 && c!=normalComponent) || c==3+normalComponent)
                result.solidComponents.push_back(c);
            else if(solid.initialGuess[(size_t)c]!=0.0)
                throw std::runtime_error("Rigid initial velocity violates the empty-plane constraint.");
    }
    return result;
}

KKTLayout::KKTLayout(
        int cellCount,
        std::vector<std::int64_t> constraintEntities,
        std::vector<int> velocityComponents,
        std::vector<int> solidComponents)
    : cells_(cellCount), constraints_(std::move(constraintEntities)),
      velocityComponents_(std::move(velocityComponents)),
      solidComponents_(std::move(solidComponents)) {
    if (cells_ <= 0 || constraints_.empty()
        || velocityComponents_.empty()) {
        throw std::runtime_error("KKTLayout received invalid block sizes.");
    }
    std::set<std::int64_t> unique;
    for (std::int64_t entity : constraints_) {
        if (entity < 0 || !unique.insert(entity).second) {
            throw std::runtime_error(
                "KKTLayout requires unique non-negative constraint ids.");
        }
    }
    std::set<int> vectorComponents;
    for (int component : velocityComponents_) {
        requireComponent(component);
        if (!vectorComponents.insert(component).second) {
            throw std::runtime_error(
                "KKTLayout requires unique velocity components.");
        }
    }
    std::set<int> uniqueSolidComponents;
    for (int component : solidComponents_) {
        if (component < 0 || component >= 6
            || !uniqueSolidComponents.insert(component).second) {
            throw std::runtime_error(
                "KKTLayout requires unique solid components in [0,5].");
        }
    }
    const std::size_t vectorCount=velocityComponents_.size();
    orderedDofs_.reserve(
        static_cast<std::size_t>(cells_)
        +vectorCount*static_cast<std::size_t>(cells_+constraints_.size())
        +solidComponents_.size());
    for (int cell=0; cell<cells_; ++cell) {
        orderedDofs_.push_back(pressure(cell));
    }
    for (int cell=0; cell<cells_; ++cell) {
        for (int component : velocityComponents_) {
            orderedDofs_.push_back(velocity(cell,component));
        }
    }
    for (int marker=0; marker<static_cast<int>(constraints_.size()); ++marker) {
        for (int component : velocityComponents_) {
            orderedDofs_.push_back(constraint(marker,component));
        }
    }
    for (int component : solidComponents_) {
        orderedDofs_.push_back(solid(component));
    }
}

LinearAlgebra::GlobalDofId KKTLayout::pressure(int cell) const {
    requireCell(cell);
    return LinearAlgebra::GlobalDofId::make(
        LinearAlgebra::GlobalDofSpace::Pressure,cell);
}

LinearAlgebra::GlobalDofId KKTLayout::velocity(
        int cell,int component) const {
    requireCell(cell); requireComponent(component);
    return LinearAlgebra::GlobalDofId::make(
        LinearAlgebra::GlobalDofSpace::Velocity,cell,component);
}

LinearAlgebra::GlobalDofId KKTLayout::constraint(
        int marker,int component) const {
    requireMarker(marker); requireComponent(component);
    return LinearAlgebra::GlobalDofId::make(
        LinearAlgebra::GlobalDofSpace::Constraint,
        constraints_[static_cast<std::size_t>(marker)],component);
}

LinearAlgebra::GlobalDofId KKTLayout::solid(int component) const {
    componentSlot(solidComponents_,component,"solid");
    return LinearAlgebra::GlobalDofId::make(
        LinearAlgebra::GlobalDofSpace::Solid,0,component);
}

std::size_t KKTLayout::pressureSlot(int cell) const {
    requireCell(cell); return static_cast<std::size_t>(cell);
}

std::size_t KKTLayout::velocitySlot(int cell,int component) const {
    requireCell(cell); requireComponent(component);
    return static_cast<std::size_t>(cells_)
        +static_cast<std::size_t>(cell)*velocityComponents_.size()
        +componentSlot(velocityComponents_,component,"velocity");
}

std::size_t KKTLayout::constraintSlot(int marker,int component) const {
    requireMarker(marker); requireComponent(component);
    return static_cast<std::size_t>(cells_)
        +static_cast<std::size_t>(cells_)*velocityComponents_.size()
        +static_cast<std::size_t>(marker)*velocityComponents_.size()
        +componentSlot(velocityComponents_,component,"constraint");
}

std::size_t KKTLayout::solidSlot(int component) const {
    return static_cast<std::size_t>(cells_)
        +static_cast<std::size_t>(cells_+constraints_.size())
            *velocityComponents_.size()
        +componentSlot(solidComponents_,component,"solid");
}

void KKTLayout::requireCell(int cell) const {
    if (cell < 0 || cell >= cells_) {
        throw std::runtime_error("KKT cell is outside its DOF block.");
    }
}

void KKTLayout::requireMarker(int marker) const {
    if (marker < 0 || marker >= static_cast<int>(constraints_.size())) {
        throw std::runtime_error("KKT marker is outside its DOF block.");
    }
}

void KKTLayout::requireComponent(int component) {
    if (component < 0 || component >= 3) {
        throw std::runtime_error("KKT vector component must be 0, 1, or 2.");
    }
}

std::size_t KKTLayout::componentSlot(
        const std::vector<int>& components,
        int component,
        const char* block) {
    const auto found=std::find(components.begin(),components.end(),component);
    if(found==components.end()) {
        throw std::runtime_error(
            std::string("KKT ")+block+" component is outside its block.");
    }
    return static_cast<std::size_t>(found-components.begin());
}

} // namespace SF::PressureBased
