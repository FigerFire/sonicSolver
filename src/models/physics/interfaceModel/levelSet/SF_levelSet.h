#pragma once

/// @file SF_levelSet.h
/// @brief Level Set 界面模型调度器。

#include "SF_interfaceModel.h"
#include "SF_state.h"

namespace SF::Physics::InterfaceModels {

/// @brief 拥有 phi、几何和物性缓存的求解器无关 Level Set 模型。
class LevelSetModel final : public Model {
public:
    /// @brief 保存并校验 Level Set 配置。
    explicit LevelSetModel(const Multiphase::MultiPhaseConfig& config);

    Representation representation() const noexcept override {
        return Representation::LevelSet;
    }
    bool initialized() const noexcept override { return initialized_; }
    const Multiphase::MultiPhaseConfig& config() const noexcept override {
        return config_;
    }
    void initialize(const Field& field) override;
    void beginTimeStep(double dt) override;
    void assembleTransportRHS(const Field& field) override;
    void completeTimeStep(Field& field, double dt) override;
    void refresh(const Field& field) override;
    void addSourceTerms(Field& field, Residual& residual) const override;
    ScalarField& primaryScalar() override { return state_.scalar(); }
    const ScalarField& primaryScalar() const override {
        return state_.scalar();
    }
    std::vector<double>& primaryScalarRHS() override { return rhs_; }
    const std::vector<double>& primaryScalarRHS() const override {
        return rhs_;
    }
    Multiphase::LevelSetField* levelSetState() noexcept override {
        return &state_;
    }
    const Multiphase::LevelSetField* levelSetState() const noexcept override {
        return &state_;
    }
    const Multiphase::ConservationDiagnostics& conservation() const override {
        return diagnostics_;
    }
    FDM::InterfacePressureJump pressureJump(
        const Field& field, int i, int j, int k, int axis) const override;

    void applyBoundary(const Field& field) override;
    void correct(const Field& field, double dt) override;
    double dynamicViscosity(const Field& field,
                            int i, int j, int k,
                            double laminarMu) const override;

private:
    Multiphase::MultiPhaseConfig config_;
    Multiphase::LevelSetField state_;
    std::vector<double> rhs_;
    Multiphase::ConservationDiagnostics diagnostics_;
    bool initialized_ = false;

    void requireInitialized(const char* operation) const;
    void updateDiagnostics(const Field& field, bool resetReference);
};

} // namespace SF::Physics::InterfaceModels
