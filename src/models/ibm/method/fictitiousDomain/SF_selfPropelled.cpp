/// @file SF_selfPropelled.cpp
/// @brief 由虚拟域零合力/零合力矩条件装配并求解三自由度小系统。

#include "method/fictitiousDomain/SF_selfPropelled.h"

#include "discrete/SF_linearSystem.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::IBM::Forcing::SelfPropelled {

Vector3 targetVelocity(
        FDM::IBMRigidMotionMode mode,
        const Vector3& generalizedVelocity,
        const Vector3& relativePosition,
        const Vector3& deformationVelocity) {
    return (mode == FDM::IBMRigidMotionMode::Rotate
        ? cross(generalizedVelocity,relativePosition)
        : generalizedVelocity)+deformationVelocity;
}

GeneralizedSystem assembleGeneralizedSystem(
        const std::vector<Sample>& samples,
        FDM::IBMRigidMotionMode mode) {
    GeneralizedSystem result;
    for (const auto& sample : samples) {
        if (!std::isfinite(sample.phaseMass) || sample.phaseMass <= 0.0) {
            throw std::runtime_error(
                "DFM self-propulsion found non-positive virtual mass.");
        }
        const Vector3 relativeFluid =
            sample.predictedVelocity-sample.deformationVelocity;
        if (mode == FDM::IBMRigidMotionMode::Motivation) {
            for (int d=0; d<3; ++d) result.matrix[(size_t)(3*d+d)]
                += sample.phaseMass;
            result.rhs = result.rhs+relativeFluid*sample.phaseMass;
            continue;
        }
        const Vector3& r = sample.relativePosition;
        const double rr = dot(r,r);
        const double coordinate[3] = {r.x,r.y,r.z};
        for (int row=0; row<3; ++row) {
            for (int column=0; column<3; ++column) {
                result.matrix[(size_t)(3*row+column)] += sample.phaseMass
                    *((row == column ? rr : 0.0)
                      -coordinate[row]*coordinate[column]);
            }
        }
        result.rhs = result.rhs+cross(r,relativeFluid)*sample.phaseMass;
    }
    return result;
}

Vector3 solveGeneralizedSystem(
        const GeneralizedSystem& system,
        FDM::IBMRigidMotionMode mode,
        const Vector3& externalForce,
        const Vector3& externalTorque,
        double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "DFM self-propulsion requires a positive finite dt.");
    }
    std::vector<std::vector<double>> matrix(3,std::vector<double>(3,0.0));
    for (int row=0;row<3;++row)
        for (int column=0;column<3;++column)
            matrix[(size_t)row][(size_t)column]=
                system.matrix[(size_t)(3*row+column)];
    const Vector3 rhs=system.rhs+(mode == FDM::IBMRigidMotionMode::Rotate
        ? externalTorque : externalForce)*dt;
    std::vector<double> solution;
    if (!Math::LinearSystem::solve(
            matrix,{rhs.x,rhs.y,rhs.z},solution)) {
        throw std::runtime_error(
            "DFM self-propulsion generalized mass matrix is singular; "
            "check body geometry and selected rigid motion DOFs.");
    }
    const Vector3 result{solution[0],solution[1],solution[2]};
    if (!std::isfinite(result.x) || !std::isfinite(result.y)
        || !std::isfinite(result.z)) {
        throw std::runtime_error(
            "DFM self-propulsion produced non-finite generalized velocity.");
    }
    return result;
}

Vector3 solveGeneralizedVelocity(
        const std::vector<Sample>& samples,
        FDM::IBMRigidMotionMode mode,
        const Vector3& externalForce,
        const Vector3& externalTorque,
        double dt) {
    return solveGeneralizedSystem(
        assembleGeneralizedSystem(samples,mode),mode,
        externalForce,externalTorque,dt);
}

} // namespace SF::IBM::Forcing::SelfPropelled
