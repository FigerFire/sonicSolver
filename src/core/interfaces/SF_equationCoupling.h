#pragma once

/// @file SF_equationCoupling.h
/// @brief 附加方程、界面跳跃与主 Navier–Stokes 方程的耦合接口。

#include "core/field/SF_field.h"
#include "core/residual/SF_residual.h"
#include "SF_transportModel.h"
#include "core/state/SF_state.h"

#include <vector>

namespace SF::FDM {

/// @brief 沿结构正方向单个面的已知压力跳跃。
struct InterfacePressureJump {
    bool crossesInterface = false;
    double fractionFromLeft = 0.0;
    double targetRightMinusLeft = 0.0;
};

/// @brief 与具体压力算法无关的锐界面压力跳跃查询接口。
class IInterfaceJumpCondition {
public:
    virtual ~IInterfaceJumpCondition() = default;
    virtual InterfacePressureJump pressureJump(
        const Field& field, int i, int j, int k, int axis) const = 0;
};

/// @brief 附加方程系统与主 NS 方程的窄耦合生命周期。
class IEquationSystemCoupling {
public:
    virtual ~IEquationSystemCoupling() = default;
    virtual void beginStep(Field&, double) {}
    virtual void prepareRHS(Field&, double) {}
    virtual void assembleRHS(Field&, Residual&, double) {}
    virtual void preparePressureCorrection(Field&) {}
    virtual void commitStep(Field&, double) {}

    virtual void beginStep(const std::vector<Field*>& fields, double dt) {
        for (Field* field : fields) if (field) beginStep(*field, dt);
    }

    virtual void prepareRHS(const std::vector<Field*>& fields, double dt) {
        for (Field* field : fields) if (field) prepareRHS(*field, dt);
    }

    virtual void assembleRHS(const std::vector<Field*>& fields,
                             const std::vector<Residual*>& residuals,
                             double dt) {
        if (fields.size() != residuals.size()) {
            throw std::runtime_error("equation coupling patch/workspace size mismatch.");
        }
        for (size_t index = 0; index < fields.size(); ++index) {
            if (fields[index] && residuals[index]) {
                assembleRHS(*fields[index], *residuals[index], dt);
            }
        }
    }

    virtual void preparePressureCorrection(
            const std::vector<Field*>& fields) {
        for (Field* field : fields) {
            if (field) preparePressureCorrection(*field);
        }
    }

    virtual void commitStep(const std::vector<Field*>& fields, double dt) {
        for (Field* field : fields) if (field) commitStep(*field, dt);
    }

    virtual State::VariableRegistry* variables(Field&) { return nullptr; }

    virtual const ITransportModel* transportModel(const Field&) const {
        return nullptr;
    }

    virtual const IInterfaceJumpCondition* interfaceJumpCondition(
            const Field&) const {
        return nullptr;
    }
};

} // namespace SF::FDM
