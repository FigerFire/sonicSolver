#pragma once

/// @file SF_transformation.h
/// @brief RawEquationSystem 到 ExecutableEquationSystem 的统一 transformation contract。
///        Transformation IR（descriptor/match/record）也在此声明。

#include "core/system/SF_transformationTypes.h"
#include "core/system/SF_solveProgram.h"

#include <memory>
#include <string>
#include <vector>

namespace SF::System {

class ExecutableEquationSystemBuilder {
public:
    explicit ExecutableEquationSystemBuilder(const RawEquationSystem& raw);

    void addUnknown(UnknownDescriptor unknown);
    void addEquation(
        EquationDescriptor descriptor, Equation::Definition definition);
    void addOperator(GeneratedOperatorDescriptor operation);
    /// @brief 声明一个 executable operation（runtime provider 必须实现它）。
    void addExecutableOperation(ExecutableOperation operation);
    void addCompiledEquation(CompiledEquation equation);

    const ExecutableEquationSystem& system() const { return system_; }
    /// @brief 已声明的 executable operations（transformer 只追加）。
    const std::vector<ExecutableOperation>& operations() const {
        return system_.operations;
    }
    ExecutableEquationSystem finish() { return std::move(system_); }

private:
    ExecutableEquationSystem system_;
};

class IEquationSystemTransformer {
public:
    virtual ~IEquationSystemTransformer() = default;
    virtual const TransformationDescriptor& descriptor() const = 0;
    virtual TransformationMatch match(const RawEquationSystem& system) const = 0;
    virtual void transform(
        const RawEquationSystem& raw,
        ExecutableEquationSystemBuilder& executable,
        std::vector<ExecutionPolicy>& policies,
        TransformationRecord& record) const = 0;
};

class TransformerRegistry {
public:
    void registerTransformer(
        std::unique_ptr<IEquationSystemTransformer> transformer);
    const IEquationSystemTransformer* find(std::string_view id) const;

private:
    std::vector<std::unique_ptr<IEquationSystemTransformer>> transformers_;
};

class TransformationPipeline {
public:
    static ExecutableEquationSystem apply(
        const RawEquationSystem& raw,
        const std::vector<TransformationDescriptor>& requests,
        const TransformerRegistry& registry,
        std::vector<ExecutionPolicy>& policies,
        std::vector<TransformationRecord>& records);
};

/// @brief Builtin pressure-constraint transformer；applicability 只检查数学 contract。
std::unique_ptr<IEquationSystemTransformer> makePressureConstraintTransformer();

/// @brief Eulerian shared-pressure constraint 使用相同 registry/pipeline。
std::unique_ptr<IEquationSystemTransformer> makeSharedPressureTransformer();

/// @brief Variational/KKT IBM constraint 使用相同 registry/pipeline。
std::unique_ptr<IEquationSystemTransformer> makeImmersedConstraintTransformer();

} // namespace SF::System
