#pragma once

/// @file SF_transferLedger.h
/// @brief 相间质量、动量、能量与外部热输入的守恒账本。

#include <array>
#include <cstddef>
#include <vector>

namespace SF::Physics::InterphaseTransfer {

/// @brief 与 Field 解耦的守恒源项缓存，存储顺序为 cell-major。
class SourceBundle {
public:
    SourceBundle() = default;
    SourceBundle(int cellCount, int variableCount);

    int cellCount() const { return cellCount_; }
    int variableCount() const { return variableCount_; }
    double& operator()(int cell, int variable);
    double operator()(int cell, int variable) const;

private:
    int cellCount_ = 0;
    int variableCount_ = 0;
    std::vector<double> values_;
};

/// @brief 每个离散点上的 phase-transfer 守恒账本。
///
/// 内部相间交换写入成对的正负源；wall/external energy 是混合物外部输入。
/// 账本只保存物理交换，不知道 FluidStateModel 的变量索引或 Field 数据布局。
class TransferLedger {
public:
    TransferLedger() = default;
    TransferLedger(int cellCount, int phaseCount);

    int cellCount() const { return cellCount_; }
    int phaseCount() const { return phaseCount_; }

    /// @brief 记录 donor -> receiver 的非负质量传递率。
    void addInternalMassTransfer(int cell, int donorPhase,
                                 int receiverPhase, double rate);

    /// @brief 记录成对守恒的相间动量传递。
    void addInternalMomentumTransfer(
        int cell, int donorPhase, int receiverPhase,
        const std::array<double, 3>& rate);

    /// @brief 记录成对守恒的相间能量传递。
    void addInternalEnergyTransfer(int cell, int donorPhase,
                                   int receiverPhase, double rate);

    /// @brief 记录相间力把宏观动能转换为相内能的贡献。
    void addInternalMechanicalEnergy(int cell,int phase,double rate);

    /// @brief 记录壁面输入，并同时写入 mixture external energy。
    void addWallEnergyInput(int cell, double rate);
    /// @brief 将壁热输入分配给指定相并写入总体外部能量。
    void addWallEnergyInput(int cell,int phase,double rate);

    /// @brief 记录非壁面外部输入，并同时写入 mixture external energy。
    void addExternalEnergyInput(int cell, double rate);

    double phaseMassSource(int cell, int phase) const;
    double phaseMomentumSource(int cell, int phase, int component) const;
    double phaseEnergySource(int cell, int phase) const;
    double mixtureEnergySource(int cell) const;
    double wallEnergyInput(int cell) const;
    double externalEnergyInput(int cell) const;
    double mechanicalEnergyConversion(int cell) const;

    /// @brief 检查内部质量/动量/能量严格成对守恒及外部能量账本闭合。
    void validate(double relativeTolerance = 1.0e-12) const;

private:
    int cellCount_ = 0;
    int phaseCount_ = 0;
    std::vector<double> mass_;
    std::vector<double> momentum_;
    std::vector<double> energy_;
    std::vector<double> mixtureEnergy_;
    std::vector<double> wallEnergy_;
    std::vector<double> externalEnergy_;
    std::vector<double> mechanicalEnergy_;
    std::vector<double> phaseExternalEnergy_;

    size_t phaseIndex(int cell, int phase) const;
    size_t momentumIndex(int cell, int phase, int component) const;
    void requireCell(int cell) const;
    void requirePhase(int phase) const;
};

} // namespace SF::Physics::InterphaseTransfer
