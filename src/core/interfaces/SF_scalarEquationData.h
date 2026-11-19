#pragma once
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include "core/state/SF_valueTypes.h"

namespace SF::FDM {
namespace Scalar {
enum class BoundMode { None, Diagnose, Project };
/// @brief 显式选择的有界处理；不由字段名或模型身份启用。
struct Bounds {bool enabled=false;double lower=0,upper=1;BoundMode mode=BoundMode::Diagnose;};
}
/// @brief 数值边界选择，与物理集合/条件分开。
struct ScalarBoundaryOptions {bool ilw=false;int order=0;Scalar::Bounds bounds;};
/// @brief 模型声明的物理初边值；不拥有数组、数学更新或时钟。
struct ScalarFieldInput {
    double initial=0;bool positiveInitial=false;
    std::vector<BCSetting<double>> initialConditions,boundaries;
};
using ScalarKnownValue = std::function<double(double,double,double,double)>;
enum class ScalarBoundaryKind { Unspecified, FixedValue, ZeroGradient };
struct ScalarBoundaryCondition {
    ScalarBoundaryKind kind=ScalarBoundaryKind::Unspecified;
    ScalarKnownValue value;
    /// @brief 显式三分量边界；标量端口继续使用 value，不广播未声明的 vector 数据。
    std::function<double(double,double,double,double,int)> componentValue;
};
/// Non-owning numerical data injection. Storage/clock remain in StateBundle.
/// Face order: x-/x+, y-/y+, z-/z+. Known functions receive x,y,z,stageTime.
/// @brief 实际网格集合上的常值边界；保留 cell 端点/ghost 闭合，不替代节点反射。
struct CellBoundaryData {
    ScalarBoundaryKind kind=ScalarBoundaryKind::Unspecified;
    std::vector<int> cells;
    std::vector<double> value;
};
struct ScalarEquationData {
    std::string target;
    std::array<bool,3> activeDimensions{true,true,true};
    std::array<ScalarBoundaryCondition,6> boundary;
    std::map<std::string,ScalarKnownValue> sources;
    std::vector<CellBoundaryData> cellBoundary;
    /// @brief 保持集合顺序的原标量边界，供所选metric/ILW实现消费。
    std::vector<BCSetting<double>> setBoundary;
};
/// @brief Caller keeps each immutable data object alive throughout execution.
struct ScalarInstanceData {
    std::string occurrence;
    std::weak_ptr<const ScalarEquationData> data;
};
}
