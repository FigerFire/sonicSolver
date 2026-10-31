#pragma once
/// WHICH occurrence selection is independent of equation math and execution order.
#include "core/system/SF_provenance.h"
#include "core/system/SF_dataContract.h"
#include "methods/math/SF_signature.h"
#include "core/config/types/SF_timeRecipe.h"
#include <optional>
#include <string>
#include <map>
#include <vector>
namespace SF::System {
/// @brief 数学算子引用与实现选择；地址来自原 AST，非另一份数学。
struct FormulaOperatorBinding {
    FormulaOperatorBinding()=default;
    FormulaOperatorBinding(std::string equation,std::string address,std::string implementation)
        :formula(std::move(equation)),occurrence(std::move(address)),provider(std::move(implementation)) {
        if(formula.empty() && !occurrence.empty())operation=Math::resolveOperator(occurrence);
    }
    std::optional<Math::Op> operation;
    std::string formula,occurrence,provider;
};
struct NumericalBinding {
    std::string equation;
    std::string method;
    std::vector<std::string> inputs;
    std::string occurrence;
    /// Selected numerical parameters; never inferred from the equation name.
    std::map<std::string,double> parameters;
    /// @brief 来源不参与数值匹配，仅解释显式选择/被覆盖默认。
    Provenance origin;
    /// @brief 局部覆盖先归并再编译；time 仅校验唯一 resolved recipe，不创建时间循环。
    std::vector<FormulaOperatorBinding> operators;
    std::optional<FDM::TimeRecipeId> time;
    bool inputsDeclared=false;
    /// @brief 用户 WHICH 选择覆盖贡献提供的默认；来源类别不决定优先级。
    bool explicitChoice=false;
};

}
