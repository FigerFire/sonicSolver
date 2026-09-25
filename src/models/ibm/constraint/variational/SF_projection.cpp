/// @file SF_projection.cpp
/// @brief 离散最小作用量/质量范数投影的局部 KKT 解。

#include "constraint/variational/SF_projection.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::IBM::Variational {
namespace {

bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

double component(const Vector3& value, int direction) {
    if (direction==0) return value.x;
    if (direction==1) return value.y;
    return value.z;
}

void setComponent(Vector3& value, int direction, double componentValue) {
    if (direction==0) value.x=componentValue;
    else if (direction==1) value.y=componentValue;
    else value.z=componentValue;
}

double dotProduct(const std::vector<double>& a,
                  const std::vector<double>& b,
                  const std::vector<unsigned char>& ownedMarkers,
                  const std::function<void(std::vector<double>&)>& globalSum) {
    if (a.size()!=b.size()
        || (!ownedMarkers.empty() && ownedMarkers.size()!=a.size())) {
        throw std::runtime_error(
            "Surface Schur CG received inconsistent vector ownership.");
    }
    double result=0.0;
    for (std::size_t i=0;i<a.size();++i) {
        if (!ownedMarkers.empty() && ownedMarkers[i]==0) continue;
        result+=a[i]*b[i];
    }
    if (globalSum) {
        std::vector<double> reduced{result};
        globalSum(reduced);
        result=reduced.front();
    }
    return result;
}

void applySurfaceSchur(
        const SurfaceSchurProjectionProblem& problem,
        const std::vector<double>& input,
        std::vector<double>& output) {
    const auto& points=problem.surface->points;
    const auto& inverseMass=*problem.inverseEulerianMass;
    std::vector<double> eulerian(inverseMass.size(),0.0);
    for (std::size_t marker=0;marker<points.size();++marker) {
        const double scaled=std::sqrt(points[marker].measure)*input[marker];
        for (const auto& edge:points[marker].interpolation) {
            if (edge.cell<0
                || static_cast<std::size_t>(edge.cell)>=inverseMass.size()) {
                throw std::runtime_error(
                    "Surface Schur operator received an invalid cell id.");
            }
            eulerian[static_cast<std::size_t>(edge.cell)]
                +=edge.value*scaled;
        }
    }
    output.assign(points.size(),0.0);
    for (std::size_t marker=0;marker<points.size();++marker) {
        double value=0.0;
        for (const auto& edge:points[marker].interpolation) {
            const std::size_t cell=static_cast<std::size_t>(edge.cell);
            value+=edge.value*inverseMass[cell]*eulerian[cell];
        }
        output[marker]=problem.dt*std::sqrt(points[marker].measure)*value;
    }
    if (problem.globalSum) problem.globalSum(output);
}

int solveCG(const SurfaceSchurProjectionProblem& problem,
            const std::vector<double>& rhs,
            std::vector<double>& solution,
            double& relativeResidual) {
    const double initialSquared=dotProduct(
        rhs,rhs,problem.ownedMarkers,problem.globalSum);
    solution.assign(rhs.size(),0.0);
    if (initialSquared==0.0) {
        relativeResidual=0.0;
        return 0;
    }
    if (!std::isfinite(initialSquared) || initialSquared<0.0) {
        throw std::runtime_error("Surface Schur CG received non-finite rhs.");
    }
    std::vector<double> residual=rhs;
    std::vector<double> direction=residual;
    std::vector<double> image(rhs.size(),0.0);
    double residualSquared=initialSquared;
    for (int iteration=1;iteration<=problem.maxIterations;++iteration) {
        applySurfaceSchur(problem,direction,image);
        const double denominator=dotProduct(
            direction,image,problem.ownedMarkers,problem.globalSum);
        if (!std::isfinite(denominator) || denominator<=0.0) {
            throw std::runtime_error(
                "Surface Schur CG lost positive definiteness at iteration "
                +std::to_string(iteration)+".");
        }
        const double alpha=residualSquared/denominator;
        for (std::size_t i=0;i<rhs.size();++i) {
            solution[i]+=alpha*direction[i];
            residual[i]-=alpha*image[i];
        }
        const double nextSquared=dotProduct(
            residual,residual,problem.ownedMarkers,problem.globalSum);
        relativeResidual=std::sqrt(nextSquared/initialSquared);
        if (!std::isfinite(relativeResidual)) {
            throw std::runtime_error(
                "Surface Schur CG produced a non-finite residual.");
        }
        if (relativeResidual<=problem.relativeTolerance) return iteration;
        const double beta=nextSquared/residualSquared;
        for (std::size_t i=0;i<rhs.size();++i) {
            direction[i]=residual[i]+beta*direction[i];
        }
        residualSquared=nextSquared;
    }
    throw std::runtime_error(
        "Surface Schur CG did not reach relativeTolerance="
        +std::to_string(problem.relativeTolerance)+" in "
        +std::to_string(problem.maxIterations)+" iterations.");
}

} // namespace

PointProjectionSolution solvePointProjection(
        const PointProjectionProblem& problem) {
    if (!finite(problem.predictedVelocity)
        || !finite(problem.targetVelocity)
        || !std::isfinite(problem.density) || problem.density <= 0.0
        || !std::isfinite(problem.measure) || problem.measure <= 0.0
        || !std::isfinite(problem.dt) || problem.dt <= 0.0) {
        throw std::runtime_error(
            "Variational IBM point projection requires finite velocities, "
            "positive density/measure and positive dt.");
    }

    PointProjectionSolution result;
    result.correctedVelocity = problem.targetVelocity;
    const Vector3 velocityIncrement =
        result.correctedVelocity - problem.predictedVelocity;
    result.multiplier = velocityIncrement*(problem.density/problem.dt);
    result.kineticIncrementFunctional = 0.5*problem.density
        *dot(velocityIncrement,velocityIncrement)
        *problem.measure/problem.dt;
    result.constraintResidual = norm(
        result.correctedVelocity-problem.targetVelocity);
    result.stationarityResidual = norm(
        velocityIncrement*(problem.density/problem.dt)-result.multiplier);
    if (!finite(result.multiplier)
        || !std::isfinite(result.kineticIncrementFunctional)
        || !std::isfinite(result.constraintResidual)
        || !std::isfinite(result.stationarityResidual)) {
        throw std::runtime_error(
            "Variational IBM point KKT produced a non-finite solution.");
    }
    return result;
}

double mechanicalEnergyIncrement(
        const Vector3& momentumIncrement,
        const Vector3& oldVelocity,
        const Vector3& newVelocity) {
    if (!finite(momentumIncrement) || !finite(oldVelocity)
        || !finite(newVelocity)) {
        throw std::runtime_error(
            "Variational IBM mechanical work received non-finite state.");
    }
    const double result = dot(
        momentumIncrement,(oldVelocity+newVelocity)*0.5);
    if (!std::isfinite(result)) {
        throw std::runtime_error(
            "Variational IBM mechanical work is non-finite.");
    }
    return result;
}

SurfaceSchurProjectionSolution solveSurfaceSchurProjection(
        const SurfaceSchurProjectionProblem& problem) {
    if (!problem.surface || !problem.inverseEulerianMass
        || problem.surface->points.empty()
        || problem.constraintError.size()!=problem.surface->points.size()
        || !std::isfinite(problem.dt) || problem.dt<=0.0
        || problem.maxIterations<=0
        || !std::isfinite(problem.relativeTolerance)
        || problem.relativeTolerance<=0.0) {
        throw std::runtime_error(
            "Surface Schur projection received an invalid system or "
            "solver control.");
    }
    if (!problem.ownedMarkers.empty()
        && problem.ownedMarkers.size()!=problem.surface->points.size()) {
        throw std::runtime_error(
            "Surface Schur projection marker ownership size differs from "
            "the constraint layout.");
    }
    for (double value:*problem.inverseEulerianMass) {
        if (!std::isfinite(value) || value<0.0) {
            throw std::runtime_error(
                "Surface Schur projection requires non-negative finite "
                "inverse Eulerian masses.");
        }
    }
    SurfaceSchurProjectionSolution result;
    result.multiplier.assign(problem.surface->points.size(),Vector3());
    for (int direction=0;direction<3;++direction) {
        std::vector<double> rhs(problem.surface->points.size(),0.0);
        for (std::size_t marker=0;marker<rhs.size();++marker) {
            const auto& point=problem.surface->points[marker];
            if (!std::isfinite(point.measure) || point.measure<=0.0
                || !finite(problem.constraintError[marker])) {
                throw std::runtime_error(
                    "Surface Schur projection received invalid marker data.");
            }
            rhs[marker]=std::sqrt(point.measure)
                *component(problem.constraintError[marker],direction);
        }
        std::vector<double> scaledMultiplier;
        double relativeResidual=0.0;
        const int iterations=solveCG(
            problem,rhs,scaledMultiplier,relativeResidual);
        result.iterations=std::max(result.iterations,iterations);
        result.maximumRelativeResidual=std::max(
            result.maximumRelativeResidual,relativeResidual);
        for (std::size_t marker=0;marker<rhs.size();++marker) {
            setComponent(result.multiplier[marker],direction,
                scaledMultiplier[marker]
                /std::sqrt(problem.surface->points[marker].measure));
        }
    }
    return result;
}

} // namespace SF::IBM::Variational
