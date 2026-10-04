#pragma once
/// @brief Immutable selected IBM mathematical/storage implementation; no runtime routing.
#include "core/interfaces/SF_immersedSystem.h"
#include "SF_methodObjects.h"
namespace SF::System {
struct CompiledImmersedContract {
    FDM::IBMForcingAlgorithm algorithm;
    FDM::IBMEnforcement enforcement;
    FDM::IBMConstraintSupport support;
    FDM::IBMSolidModel solid;
    FDM::IBMRepresentation representation;
    FDM::IBMRigidMotionMode rigidMotionMode;
    TargetKind kind;
    std::string equation;
    std::string target;
    std::vector<CompiledMathRef> members;
};
/// @brief 注册已有 impulse/projection/KKT 实现，不创建第二个 timestep。
void addImmersedMethods(ProviderRegistry& registry);
/// @brief 编译期拒绝缺少相端口/边界与组合 schedule 的数学系统。
void validateImmersedComposition(const RawEquationSystem& raw);
/// @brief 校验 bound IBM port 与已冻结数学/STATE/method 来自同一实现。
void validateImmersedBindings(const CompiledSolvePlan& plan,const FDM::IImmersedSystem* system,
                              bool boundaryPresent);
}
