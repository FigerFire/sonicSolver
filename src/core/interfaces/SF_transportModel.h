#pragma once

/// @file SF_transportModel.h
/// @brief 黏性、湍流及界面物性向求解器提供的统一输运闭合接口。

#include "core/field/SF_field.h"
#include "core/state/SF_valueTypes.h"

#include "SF_scalarEquationData.h"
#include <array>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {

/// @brief 求解器可见的有效输运闭合接口。
struct TransportStageCoefficients {
    double dynamicViscosity=0;
    std::array<double,2> diffusivity{},production{},sink{};
    double crossSource=0,blending1=0,blending2=0;
};
using TransportStageSamples=std::array<std::function<double(int,int,int)>,2>;
/// @brief 原位两标量输运的局部物理闭合；无dt、无状态写回或时间循环。
struct TransportCellClosure {
    double density=0;
    std::array<double,2> values{},source{},diffusionWeight{};
    double crossDiffusion=0;
};
/// @brief 借用原模型数组/边界定律；数值遍历、离散及dt写回属于solver。
struct OnceStepTransportData {
    std::array<std::vector<double>*,2> variables{};
    const std::vector<double>* eddyViscosity=nullptr;
    std::array<const std::vector<BCSetting<double>>*,2> boundary{};
    std::array<double,2> floors{};
    double molecularViscosity=0;
    bool molecularDiffusion=false;
    std::function<TransportCellClosure(const Field&,int,int,int)> cellClosure;
};
class ITransportModel {
public:
    virtual ~ITransportModel() = default;
    virtual void applyBoundary(const Field& field) = 0;
    /// @brief 刷新代数闭合，不积分输运变量。
    virtual void refreshClosure(const Field&) {
        throw std::runtime_error("Transport model has no algebraic closure publication port.");
    }
    /// @brief optional原位两标量端口；空端口不授权数值推进。
    virtual OnceStepTransportData onceStepTransport() { return {}; }
    virtual double dynamicViscosity(const Field& field,
                                    int i, int j, int k,
                                    double laminarMu) const = 0;

    /// @brief 只提供指定 stage 的闭合，不积分输运未知量。
    virtual std::vector<std::string> stageScalarTargets() const {return {};}
    virtual TransportStageCoefficients stageCoefficients(const Field&,const TransportStageSamples&,int,int,int) const {
        throw std::runtime_error("Transport model has no stage coefficient implementation.");
    }
    virtual ScalarEquationData stageScalarData(const Field&,const std::string&) const {
        throw std::runtime_error("Transport model has no stage scalar boundary port.");
    }
    virtual void publishStageViscosity(const Field&,const std::vector<double>&) {
        throw std::runtime_error("Transport model has no stage viscosity publication.");
    }
    /// @brief 输运/闭合前需要新鲜 halo 的模型自管字段。
    virtual std::vector<std::string> distributedReadFields() const {
        return {};
    }
    /// @brief 输运/闭合及边界处理会写入的模型自管字段。
    virtual std::vector<std::string> distributedWriteFields() const {
        return {};
    }
    /// @brief 上述模型 stencil 的实际 halo 深度；由模型而非注册点声明。
    virtual int distributedHaloDepth() const { return 0; }
};

} // namespace SF::FDM
