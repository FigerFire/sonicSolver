/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once
#include "core/interfaces/SF_immersedTurbulenceBoundary.h"

/// @file SF_turbulence.h
/// @brief 湍流模型调度层与公共状态定义。
///
/// 设计目标:
/// 1. Solver Algorithm 只依赖一个小接口 `ITransportModel`
/// 2. RAS / LES / DNS 具体公式留在各自子目录
/// 3. 湍流标量状态独立于 `SF::Field`，避免污染主守恒变量存储

#include "SF_field.h"
#include "SF_config.h"
#include "core/interfaces/SF_transportModel.h"
#include "core/state/SF_distributedField.h"

#include <memory>
#include <string>
#include <vector>

namespace SF {
namespace Turbulence {

/// @brief 湍流标量槽位枚举。
enum class ScalarSlot {
    K,
    Epsilon,
    Omega,
    EddyMu
};

/// @brief 湍流模块自管的标量场存储。
class ScalarFields {
public:
    /// @brief 按主流场尺寸分配全部湍流标量数组。
    /// @param field 主流场，提供维度与虚胞层数。
    void resizeLike(const Field& field);

    /// @brief 判断是否已经完成尺寸分配。
    /// @return 已分配返回true。
    bool empty() const { return totalSize_ == 0; }

    /// @brief 读取/写入 k。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位引用。
    double& K(int i, int j, int k);
    /// @brief 读取 k。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位值。
    double K(int i, int j, int k) const;

    /// @brief 读取/写入 epsilon。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位引用。
    double& Epsilon(int i, int j, int k);
    /// @brief 读取 epsilon。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位值。
    double Epsilon(int i, int j, int k) const;

    /// @brief 读取/写入 omega。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位引用。
    double& Omega(int i, int j, int k);
    /// @brief 读取 omega。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位值。
    double Omega(int i, int j, int k) const;

    /// @brief 读取/写入湍流附加动力粘度缓存。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位引用。
    double& EddyMu(int i, int j, int k);
    /// @brief 读取湍流附加动力粘度缓存。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位值。
    double EddyMu(int i, int j, int k) const;

    /// @brief 读取指定槽位，并在越界时夹回最近内部实胞。
    /// @param slot 标量槽位。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 夹取后的标量值。
    double clamped(ScalarSlot slot, int i, int j, int k) const;

    /// @brief 读取指定槽位的可写引用。
    /// @param slot 标量槽位。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位引用。
    double& slot(ScalarSlot slot, int i, int j, int k);

    /// @brief 读取指定槽位。
    /// @param slot 标量槽位。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 对应槽位值。
    double slot(ScalarSlot slot, int i, int j, int k) const;

    /// @brief 返回指定槽位的连续存储，供统一分布式状态注册。
    std::vector<double>& values(ScalarSlot slot);
    const std::vector<double>& values(ScalarSlot slot) const;

    /// @brief 总单元数。
    /// @return `MX*MY*MZ`。
    int totalSize() const { return totalSize_; }

    /// @brief 虚胞层数。
    /// @return NG。
    int NG() const { return ng_; }

    /// @brief 含虚胞总维度。
    /// @return MX / MY / MZ。
    int MX() const { return mx_; }
    int MY() const { return my_; }
    int MZ() const { return mz_; }

private:
    int mx_ = 0;
    int my_ = 0;
    int mz_ = 0;
    int ng_ = 0;
    int totalSize_ = 0;

    std::vector<double> k_;
    std::vector<double> epsilon_;
    std::vector<double> omega_;
    std::vector<double> eddyMu_;

    /// @brief 计算平铺索引。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 0-based 平铺索引。
    int idx(int i, int j, int k) const;
};

/// @brief 湍流具体模型接口。
class IModel {
public:
    virtual ~IModel() = default;

    /// @brief 返回模型所属大类名。
    /// @return `"RAS"` / `"LES"` / `"DNS"`。
    virtual const char* familyName() const = 0;

    /// @brief 返回具体模型名。
    /// @return 例如 `"kEpsilon"`。
    virtual const char* modelName() const = 0;

    /// @brief 初始化模型自管状态。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    virtual void initialize(const Field& flow,
                            ScalarFields& state,
                            const FDM::TurbulenceConfig& config) = 0;

    /// @brief 模型floor与immersed wall数学；不执行通用节点BC。
    virtual void applyPhysicalBoundary(const Field&,ScalarFields&,const FDM::TurbulenceConfig&) = 0;
    /// @brief 发布代数黏度，不接受dt、不推进输运变量。
    virtual void refreshClosure(const Field&,ScalarFields&,const FDM::TurbulenceConfig&) = 0;
    /// @brief 只读当前live数组计算局部闭合，不写回或调用输运离散。
    virtual FDM::TransportCellClosure onceStepClosure(const Field&,const ScalarFields&,
        const FDM::TurbulenceConfig&,int,int,int) const {
        throw std::runtime_error("Selected model has no once-step transport closure.");
    }

    virtual FDM::TransportStageCoefficients stageCoefficients(const Field&,const FDM::TransportStageSamples&,
        const FDM::TurbulenceConfig&,int,int,int) const {
        throw std::runtime_error("Selected turbulence model has no stage closure.");
    }
    /// @brief 读取当前单元的湍流附加动力粘度。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 湍流附加动力粘度 mu_t。
    virtual double eddyDynamicViscosity(const Field& flow,
                                        const ScalarFields& state,
                                        const FDM::TurbulenceConfig& config,
                                        int i, int j, int k) const = 0;
};

/// @brief 统一的湍流管理器，对外实现 `ITransportModel`。
class Manager final : public FDM::ITransportModel {
public:
    /// @brief 从强类型配置构造管理器。
    /// @param config 湍流配置值对象。
    explicit Manager(FDM::TurbulenceConfig config,const FDM::IImmersedTurbulenceBoundary* wall=nullptr);

    /// @brief 根据配置创建模型及原数组初值；随后由 solver 初始化 BC/闭合。
    /// @param field 主流场。
    /// @return 成功创建有效模型时返回true。
    bool initialize(const Field& field);

    /// @brief 当前是否真的启用了某个湍流模型。
    /// @return 仅当 `enabled=true` 且模型创建成功时返回true。
    bool active() const;

    /// @brief 返回 `"RAS/kEpsilon"` 这类可读描述。
    /// @return 便于日志输出的模型说明。
    std::string description() const;

    /// @brief 应用模型 floors 和壁面关系；通用 BC 由 solver 执行。
    /// @param field 主流场。
    void applyBoundary(const Field& field) override;

    /// @brief 仅刷新原mu_t闭合数组。
    void refreshClosure(const Field& field) override;
    /// @brief 原k/epsilon/omega数组及只读cell闭合；调用方不获得第二份storage。
    FDM::OnceStepTransportData onceStepTransport() override;

    /// @brief 查询有效动力粘度。
    /// @param field 主流场。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @param laminarMu 基础层流动力粘度。
    /// @return `laminarMu + mu_t`。
    double dynamicViscosity(const Field& field,
                            int i, int j, int k,
                            double laminarMu) const override;
    std::vector<std::string> stageScalarTargets() const override;
    FDM::TransportStageCoefficients stageCoefficients(const Field&,const FDM::TransportStageSamples&,int,int,int) const override;
    FDM::ScalarEquationData stageScalarData(const Field&,const std::string&) const override;
    void publishStageViscosity(const Field&,const std::vector<double>&) override;
    std::vector<std::string> distributedReadFields() const override;
    std::vector<std::string> distributedWriteFields() const override;
    int distributedHaloDepth() const override;
    /// @brief Publish non-owning views of the existing Manager arrays; no second storage.
    void registerDistributed(State::DistributedFieldRegistry& registry,
                             int blockId, Field& geometry);

    /// @brief 只读访问湍流标量场，用于诊断和VTK输出。
    /// @return 当前模型自管的 k / epsilon / omega / mu_t 标量缓存。
    const ScalarFields& scalarFields() const { return state_; }
    ScalarFields& scalarFields() { return state_; }

private:
    bool stageClosurePublished_=false;
    const FDM::IImmersedTurbulenceBoundary* wall_;
    FDM::TurbulenceConfig config_;
    ScalarFields state_;
    std::unique_ptr<IModel> model_;
};

/// @brief 从 IC 配置写入湍流标量初值。
/// @param flow 主流场，仅用于 set 索引查找。
/// @param settings 对应标量的初始条件列表。
/// @param state 湍流标量存储。
/// @param slot 目标标量槽位。
void applyScalarInitialConditions(const Field& flow,
                                  const std::vector<BCSetting<double>>& settings,
                                  ScalarFields& state,
                                  ScalarSlot slot);

/// @brief 保证标量为正值。
/// @param value 原始值。
/// @param floorValue 正值下界。
/// @return `max(value, floorValue)`。
double clampPositive(double value, double floorValue);

/// @brief 计算当前单元的局部应变率模。
/// @param flow 主流场。
/// @param i i方向索引。
/// @param j j方向索引。
/// @param k k方向索引。
/// @return `|S| = sqrt(2 Sij Sij)`。
double strainRateMagnitude(const Field& flow, int i, int j, int k);

/// @brief 估算 LES 所需的局部滤波尺度。
/// @param flow 主流场。
/// @param i i方向索引。
/// @param j j方向索引。
/// @param k k方向索引。
/// @param filterScale 额外尺度倍率。
/// @return 近似网格滤波尺度。
double characteristicFilterWidth(const Field& flow,
                                 int i, int j, int k,
                                 double filterScale);


} // namespace Turbulence
} // namespace SF
