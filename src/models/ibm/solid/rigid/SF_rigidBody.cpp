/// @file SF_rigidBody.cpp
/// @brief 刚体构型、速度、惯量方程及 Bhalla 广义速度视图的实现。

#include "solid/rigid/SF_rigidBody.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing {
namespace {

std::array<double,9> rotation(const Vector3& omega, double dt) {
    const double magnitude = norm(omega);
    if (magnitude == 0.0) return {{1,0,0,0,1,0,0,0,1}};
    const Vector3 axis = omega*(1.0/magnitude);
    const double c = std::cos(magnitude*dt);
    const double s = std::sin(magnitude*dt);
    const double one = 1.0-c;
    return {{
        c+axis.x*axis.x*one,
        axis.x*axis.y*one-axis.z*s,
        axis.x*axis.z*one+axis.y*s,
        axis.y*axis.x*one+axis.z*s,
        c+axis.y*axis.y*one,
        axis.y*axis.z*one-axis.x*s,
        axis.z*axis.x*one-axis.y*s,
        axis.z*axis.y*one+axis.x*s,
        c+axis.z*axis.z*one}};
}

Vector3 applyRotation(const std::array<double,9>& matrix,
                      const Vector3& value) {
    return {
        matrix[0]*value.x+matrix[1]*value.y+matrix[2]*value.z,
        matrix[3]*value.x+matrix[4]*value.y+matrix[5]*value.z,
        matrix[6]*value.x+matrix[7]*value.y+matrix[8]*value.z};
}

Vector3 applyTranspose(const std::array<double,9>& matrix,
                       const Vector3& value) {
    return {
        matrix[0]*value.x+matrix[3]*value.y+matrix[6]*value.z,
        matrix[1]*value.x+matrix[4]*value.y+matrix[7]*value.z,
        matrix[2]*value.x+matrix[5]*value.y+matrix[8]*value.z};
}

std::array<double,9> product(const std::array<double,9>& first,
                             const std::array<double,9>& second) {
    std::array<double,9> result{};
    for (int row=0; row<3; ++row) {
        for (int column=0; column<3; ++column) {
            for (int inner=0; inner<3; ++inner) {
                result[(size_t)(3*row+column)] +=
                    first[(size_t)(3*row+inner)]
                    *second[(size_t)(3*inner+column)];
            }
        }
    }
    return result;
}

Vector3 vector(const GeoProcessing::Point& point) {
    return {point.x,point.y,point.z};
}

GeoProcessing::Point point(const Vector3& value) {
    return {value.x,value.y,value.z};
}

} // namespace

void RigidBodyState::configure(const FDM::IBMForcingConfig& config) {
    config_ = config;
    center_ = config.centerOfMass;
    linearVelocity_ = config.linearVelocity;
    angularVelocity_ = config.angularVelocity;
    orientation_ = {{1,0,0,0,1,0,0,0,1}};
    if (config.rigidMotionMode == FDM::IBMRigidMotionMode::Motivation
        && norm(config.angularVelocity) != 0.0) {
        throw std::runtime_error(
            "motionMode motivation requires angularVelocity (0 0 0).");
    }
    if (config.rigidMotionMode == FDM::IBMRigidMotionMode::Rotate
        && norm(config.linearVelocity) != 0.0) {
        throw std::runtime_error(
            "motionMode rotate requires linearVelocity (0 0 0).");
    }
}

Vector3 RigidBodyState::transformed(
        const Vector3& reference, double targetTime) const {
    const Vector3 relative = reference-config_.centerOfMass;
    if (config_.motion == FDM::IBMSolidMotion::PrescribedRigid) {
        if (config_.rigidMotionMode
            == FDM::IBMRigidMotionMode::Motivation) {
            return reference+config_.linearVelocity*targetTime;
        }
        return config_.centerOfMass
            + applyRotation(
                rotation(config_.angularVelocity,targetTime),relative);
    }
    return center_+applyRotation(orientation_,relative);
}

std::vector<GeoProcessing::Triangle> RigidBodyState::triangles(
        const GeoProcessing::STLGeometry& geometry,
        double targetTime) const {
    std::vector<GeoProcessing::Triangle> result = geometry.triangles();
    for (auto& triangle : result) {
        for (auto& vertex : triangle.v) {
            vertex = point(transformed(vector(vertex),targetTime));
        }
        const Vector3 first = vector(triangle.v[0]);
        const Vector3 second = vector(triangle.v[1]);
        const Vector3 third = vector(triangle.v[2]);
        triangle.normal = point(normalize(cross(second-first,third-first)));
    }
    return result;
}

GeoProcessing::Point RigidBodyState::referencePoint(
        const Vector3& world, double targetTime) const {
    if (!std::isfinite(targetTime)) {
        throw std::runtime_error(
            "RigidBodyState requires a finite targetTime.");
    }
    Vector3 referenceRelative;
    if (config_.motion == FDM::IBMSolidMotion::PrescribedRigid) {
        if (config_.rigidMotionMode
            == FDM::IBMRigidMotionMode::Motivation) {
            referenceRelative = world
                - (config_.centerOfMass
                   + config_.linearVelocity*targetTime);
        } else {
            referenceRelative = applyTranspose(
                rotation(config_.angularVelocity,targetTime),
                world-config_.centerOfMass);
        }
    } else {
        referenceRelative = applyTranspose(
            orientation_,world-center_);
    }
    return point(config_.centerOfMass+referenceRelative);
}

Vector3 RigidBodyState::velocityAt(
        const Vector3& pointValue, double targetTime) const {
    const Vector3 currentCenter = center(targetTime);
    if (config_.rigidMotionMode == FDM::IBMRigidMotionMode::Motivation) {
        return linearVelocity_;
    }
    return cross(angularVelocity_,pointValue-currentCenter);
}

Vector3 RigidBodyState::deformationVelocityAt(
        const Vector3&, double targetTime) const {
    if (!std::isfinite(targetTime)) {
        throw std::runtime_error(
            "RigidBodyState requires a finite targetTime.");
    }
    return {};
}

void RigidBodyState::prepareEquationView(
        FDM::ImmersedSurfaceSystem& surface,
        double targetTime,
        double dt) const {
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "RigidBodyState equation view requires finite targetTime and "
            "positive dt.");
    }
    const Vector3 unit[3] = {
        {1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    const bool selfPropelled =
        config_.solidModel == FDM::IBMSolidModel::SelfPropelledRigid;
    const bool unknownRigidVelocity =
        config_.motion == FDM::IBMSolidMotion::CoupledRigid;
    for (auto& marker : surface.points) {
        marker.prescribedVelocity = selfPropelled
            ? deformationVelocityAt(marker.position,targetTime)
            : unknownRigidVelocity ? Vector3()
                                   : velocityAt(marker.position,targetTime);
        for (int q=0; q<3; ++q) {
            marker.solidVelocityBasis[(size_t)q] = unit[q];
            marker.solidVelocityBasis[(size_t)(q+3)] =
                cross(unit[q],marker.relativePosition);
        }
    }

    auto& equation = surface.solidEquation;
    equation = {};
    equation.solveGeneralizedVelocity =
        config_.motion == FDM::IBMSolidMotion::CoupledRigid;
    if(equation.solveGeneralizedVelocity) {
        equation.activeComponents = config_.rigidMotionMode
            == FDM::IBMRigidMotionMode::Motivation
            ? std::vector<int>{0,1,2} : std::vector<int>{3,4,5};
    }
    equation.generalizedDofs = (int)equation.activeComponents.size();
    equation.initialGuess = {
        linearVelocity_.x,linearVelocity_.y,linearVelocity_.z,
        angularVelocity_.x,angularVelocity_.y,angularVelocity_.z};
    if (!equation.solveGeneralizedVelocity) return;

    if (!std::isfinite(config_.rigidMass) || config_.rigidMass <= 0.0) {
        throw std::runtime_error(
            "Coupled six-DOF rigid body requires positive mass.");
    }
    for (int component=0; component<3; ++component) {
        equation.lhs[(size_t)(6*component+component)] =
            config_.rigidMass/dt;
    }
    equation.rhs[0]=config_.rigidMass*linearVelocity_.x/dt
        +config_.externalForce.x;
    equation.rhs[1]=config_.rigidMass*linearVelocity_.y/dt
        +config_.externalForce.y;
    equation.rhs[2]=config_.rigidMass*linearVelocity_.z/dt
        +config_.externalForce.z;

    const auto inertia = worldInertia(targetTime);
    Vector3 angularMomentum;
    for (int row=0; row<3; ++row) {
        double value = 0.0;
        for (int column=0; column<3; ++column) {
            const double entry = inertia[(size_t)(3*row+column)];
            if (!std::isfinite(entry)) {
                throw std::runtime_error(
                    "Coupled rotating rigid body has invalid inertia.");
            }
            equation.lhs[(size_t)(6*(row+3)+column+3)] = entry/dt;
            value += entry*(column == 0 ? angularVelocity_.x
                            : column == 1 ? angularVelocity_.y
                                          : angularVelocity_.z);
        }
        if (row == 0) angularMomentum.x = value;
        else if (row == 1) angularMomentum.y = value;
        else angularMomentum.z = value;
    }
    const Vector3 angularRhs = angularMomentum*(1.0/dt)
        +config_.externalTorque-cross(angularVelocity_,angularMomentum);
    equation.rhs[3]=angularRhs.x;
    equation.rhs[4]=angularRhs.y;
    equation.rhs[5]=angularRhs.z;
}

Vector3 RigidBodyState::center(double targetTime) const {
    if (config_.motion == FDM::IBMSolidMotion::PrescribedRigid
        && config_.rigidMotionMode
            == FDM::IBMRigidMotionMode::Motivation) {
        return config_.centerOfMass+config_.linearVelocity*targetTime;
    }
    return center_;
}

std::array<double,9> RigidBodyState::worldInertia(
        double targetTime) const {
    std::array<double,9> current=orientation_;
    if (config_.motion == FDM::IBMSolidMotion::PrescribedRigid
        && config_.rigidMotionMode == FDM::IBMRigidMotionMode::Rotate) {
        current=rotation(config_.angularVelocity,targetTime);
    }
    const double principal[3]={config_.principalInertia.x,
                               config_.principalInertia.y,
                               config_.principalInertia.z};
    std::array<double,9> result{};
    for(int row=0;row<3;++row)
        for(int column=0;column<3;++column)
            for(int axis=0;axis<3;++axis)
                result[(size_t)(3*row+column)]+=
                    current[(size_t)(3*row+axis)]*principal[axis]
                    *current[(size_t)(3*column+axis)];
    return result;
}

void RigidBodyState::advanceCoupled(
        const Vector3& linearVelocity,
        const Vector3& angularVelocity,
        double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "RigidBodyState requires finite positive dt.");
    }
    if ((config_.rigidMotionMode == FDM::IBMRigidMotionMode::Motivation
         && norm(angularVelocity) != 0.0)
        || (config_.rigidMotionMode == FDM::IBMRigidMotionMode::Rotate
            && norm(linearVelocity) != 0.0)) {
        throw std::runtime_error("Coupled rigid velocity violates motionMode.");
    }
    linearVelocity_ = linearVelocity;
    angularVelocity_ = angularVelocity;
    center_ = center_+linearVelocity_*dt;
    orientation_ = product(rotation(angularVelocity_,dt),orientation_);
}

} // namespace SF::IBM::Forcing
