#pragma once

/// @file SF_stateRegistry.h
/// @brief STATE — base variables and their physical storage/closure contracts.
/// Call input versions do not choose equation destinations, order or numerical recipes.

#include "core/system/SF_formula.h"
#include "core/system/SF_dataContract.h"
#include <algorithm>
#include <optional>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace SF::System {
struct StateProducer;
/// @brief 跨轴身份显式区分命名空间，不根据名字或前缀猜测引用类别。
enum class ReferenceKind { Equation, Step };
/// @brief STATE 按 HOW 调用身份声明读取；唯一方程调用可用方程 ID 简写。
struct CallInputs {
    std::string call;
    std::vector<StateUse> reads;
    Provenance origin;
    ReferenceKind reference=ReferenceKind::Equation;
};
enum class VariableLocation {
    EulerianCell,
    EulerianFace,
    BodyConstraint,
    SurfaceConstraint,
    SolidGlobal
};

enum class OwnershipKind {
    EulerianGlobalDof,
    CanonicalFace,
    ConstraintGlobalDof,
    SolidGlobalDof
};

enum class ValueShape { Scalar, Vector, Tensor };
enum class StateRole { Primary, Transported, Algebraic, Multiplier, Derived };
enum class StorageBinding {
    PackedDistributed,
    NamedDistributed,
    /// Provider-owned physical arrays, exposed through non-owning distributed views.
    ProviderDistributed,
    TransientWorkspace,
    SpecializedExecutor
};

/// @brief 物化 workspace 的明确边界准备；不是物理 BC 或 MPI halo。
enum class WorkspaceBoundary { None, NearestInterior };

/// @brief A base symbol; Working/Correction/Stage are compiler-created views.
struct StateSymbol {
    std::string id;
    std::string name;
    VariableLocation location = VariableLocation::EulerianCell;
    int components = 1;
    OwnershipKind ownership = OwnershipKind::EulerianGlobalDof;
    ValueShape shape = ValueShape::Scalar;
    Math::Signature signature;
    FormulaExpr expression() const {return FormulaExpr::symbol(id,signature);}
    StateRole role = StateRole::Primary;
    StorageBinding storageBinding = StorageBinding::SpecializedExecutor;
    std::string storageKey;
    int componentOffset = 0;
    /// A derived constant has no writable runtime field; its value is frozen at composition.
    std::optional<double> constantValue;
    bool initializationRequired = true;
    bool boundaryRequired = true;
    bool restartEligible = true;
    bool outputEligible = true;
    bool runtimeStorageRequired = true;
    std::string nameSpace;
    Provenance origin;
    StateVersion availableVersion=StateVersion::Current;
    StateEvaluation evaluation=StateEvaluation::Direct;
    std::vector<std::string> dependencies;
    /// @brief 模块提供的不可变生产者；不持有物理数组、不选择更新顺序。
    std::shared_ptr<const StateProducer> producer;
    /// @brief 仅选定数值 workspace 请求；先求内部最终值，再准备空间读取区域。
    WorkspaceBoundary workspaceBoundary=WorkspaceBoundary::None;
    bool readOnly() const { return constantValue.has_value() || producer; }
};

/// Active case base variables, independently of the available builtin catalog.
/// Registration neither allocates runtime arrays nor schedules equations.
class StateRegistry {
public:
    void bindInputs(CallInputs value) {
        if(value.call.empty())throw std::runtime_error("Empty STATE call input binding.");
        if(std::any_of(callInputs_.begin(),callInputs_.end(),[&](const auto& previous){return previous.call==value.call && previous.reference==value.reference;}))
            throw std::runtime_error("Duplicate STATE call input binding: "+value.call+" from "+value.origin.source);
        for(const auto& read:value.reads) {
            if(read.symbol.empty() || read.symbol.find_first_of("*'")!=std::string::npos)
                throw std::runtime_error("Invalid STATE input symbol: "+read.symbol);
            if(std::count_if(value.reads.begin(),value.reads.end(),[&](const auto& other){return other.symbol==read.symbol;})!=1)
                throw std::runtime_error("Duplicate input read declaration: "+read.symbol+" at "+value.call);
        }
        callInputs_.push_back(std::move(value));
    }
    const std::vector<CallInputs>& callInputs() const {return callInputs_;}
    void add(StateSymbol symbol) {
        Math::requireKnown(symbol.signature);
        const int rank=symbol.shape==ValueShape::Scalar?0:symbol.shape==ValueShape::Vector?1:2;
        if(symbol.signature.tupleSize || symbol.signature.rank!=rank)throw std::runtime_error("STATE shape/signature mismatch: "+symbol.id);
        const auto extent=static_cast<long long>(symbol.signature.extent);
        const auto components=rank==0?1:rank==1?extent:
            symbol.signature.structure==Math::TensorStructure::Symmetric?extent*(extent+1)/2:extent*extent;
        if(symbol.components!=components)throw std::runtime_error("STATE component/signature mismatch: "+symbol.id);

        if (symbol.id.empty() || symbol.id.find_first_of("*'")!=std::string::npos
            || symbol.components<=0
            || (symbol.shape==ValueShape::Scalar && symbol.components!=1))
            throw std::runtime_error("Invalid base STATE symbol: "+symbol.id);
        if(symbol.workspaceBoundary!=WorkspaceBoundary::None &&
            (symbol.workspaceBoundary!=WorkspaceBoundary::NearestInterior
                || symbol.storageBinding!=StorageBinding::TransientWorkspace
                || symbol.evaluation!=StateEvaluation::Materialized
                || symbol.producer || symbol.constantValue))
            throw std::runtime_error("Invalid explicit workspace boundary preparation: "+symbol.id);
        if (contains(symbol.id)) throw std::runtime_error("Duplicate STATE symbol: "+symbol.id);
        symbols_.push_back(std::move(symbol));
    }
    bool contains(std::string_view id) const {
        return std::any_of(symbols_.begin(),symbols_.end(),[&](const auto& value) { return value.id==id; });
    }
    const StateSymbol& at(std::string_view id) const {
        const auto found=std::find_if(symbols_.begin(),symbols_.end(),[&](const auto& value) { return value.id==id; });
        if (found==symbols_.end()) throw std::runtime_error("Unknown base STATE symbol: "+std::string(id));
        return *found;
    }
    /// @brief 只核验已注册STATE声明的引用一致性，不重新推导AST数学。
    void validateSymbols(const FormulaExpr& expression) const {
        if(expression.kind()==FormulaExpr::Kind::Symbol && contains(expression.name()))
            Math::requireCompatible(expression.signature(),at(expression.name()).signature,"STATE reference "+expression.name());
        for(const auto& a:expression.arguments())validateSymbols(a);
    }
    const std::vector<StateSymbol>& symbols() const { return symbols_; }
    FormulaExpr expression(std::string_view id) const {return at(id).expression();}
    void selectSolution(std::string id) {
        if (!contains(id)) throw std::runtime_error("Solution STATE has no metadata: "+id);
        if (isSolution(id)) throw std::runtime_error("Duplicate solution STATE: "+id);
        solutionVariables_.push_back(std::move(id));
    }
    bool isSolution(std::string_view id) const {
        return std::find(solutionVariables_.begin(),solutionVariables_.end(),id)!=solutionVariables_.end();
    }
    const std::vector<std::string>& solutionVariables() const { return solutionVariables_; }
    void setSelectionOrigin(std::string value) { selectionOrigin_=std::move(value); }
    const std::string& selectionOrigin() const { return selectionOrigin_; }
    std::size_t size() const { return symbols_.size(); }
private:
    std::vector<CallInputs> callInputs_;
    std::vector<StateSymbol> symbols_;
    std::vector<std::string> solutionVariables_;
    std::string selectionOrigin_;
};

} // namespace SF::System
