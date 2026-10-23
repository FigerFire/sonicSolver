#pragma once

/// @file SF_stateViews.h
/// @brief Compiler output: lazily requested versions of base STATE symbols.
#include "core/system/SF_stateRegistry.h"
#include <functional>
#include <typeindex>

namespace SF { class Field; namespace State { struct DistributedFieldView; } }

namespace SF::System {
/// @brief Stage/OldTime are temporal lowering output, never user registrations.
enum class StateViewKind { Physical, Working, Correction, OldTime, Stage, Workspace };
/// @brief Only one storage authority supplies a view; aliases never own another Q.
enum class StateViewOwner { PhysicalState, CompilerWorkspace, NumericalProvider };

/// @brief 只绑定现有数组/闭合；返回的更新只能由显式 HOW operation 调用。
struct StateProducer {
    std::string id;
    FormulaExpr expression;
    std::vector<StateViewKind> versions{StateViewKind::Physical};
    /// Lazy 仅允许既有、声明过失效契约的只读 closure，不在此插入刷新/同步。
    std::function<std::function<void()>(Field&,
        const std::vector<State::DistributedFieldView>&,State::DistributedFieldView&)> bind;
};
inline bool producedBy(const StateSymbol& state,std::string_view id) {
    return state.producer && state.producer->id==id;
}

/// @brief No ordering data lives in the state-view allocation contract.
struct CompiledStateView {
    std::string symbol;
    StateViewKind kind=StateViewKind::Physical;
    int stage=-1;
    std::string storage;
    int componentOffset=0;
    int components=1;
    VariableLocation location=VariableLocation::EulerianCell;
    OwnershipKind ownership=OwnershipKind::EulerianGlobalDof;
    StateViewOwner owner=StateViewOwner::PhysicalState;
    /// @brief 引用原 STATE 生产者，不复制数学或选择第二种实现。
    std::shared_ptr<const StateProducer> producer;
    bool haloRequired=false;
};

/// @brief 非数学资源的编译需求；不伪装成 Field 或 STATE symbol。
struct ResourceViewRequirement {
    std::string id,owner;
    std::type_index type{typeid(void)};
    bool writable=false;
};

/// @brief 借用真实 owner 的资源；访问时核验类型、权限及 owner lifetime。
class ResourceView {
    ResourceViewRequirement requirement_;
    std::weak_ptr<void> owner_;
    void* object_=nullptr;
    std::function<void()> validate_;
public:
    const ResourceViewRequirement& requirement() const {return requirement_;}
    /// @brief 不刷新或分配资源；仅验证实际 owner 仍存活和已声明 availability。
    void requireAvailable() const {
        if (!object_ || owner_.expired()) throw std::runtime_error("Unavailable resource owner: "+requirement_.id);
        if (validate_) validate_();
    }
    template<class T,class Owner>
    static ResourceView bind(ResourceViewRequirement requirement,
            const std::shared_ptr<Owner>& owner,T& object,std::string_view ownerId,
            std::function<void()> validate={}) {
        if (!owner || requirement.id.empty() || requirement.owner!=ownerId || requirement.type!=typeid(T))
            throw std::runtime_error("Invalid resource view type/owner: "+requirement.id+" owner="+requirement.owner);
        ResourceView view;view.requirement_=std::move(requirement);view.owner_=owner;view.object_=&object;
        view.validate_=std::move(validate);return view;
    }
    template<class T> std::shared_ptr<const T> read() const { return access<const T>(false); }
    template<class T> std::shared_ptr<T> write() const { return access<T>(true); }
private:
    template<class T> std::shared_ptr<T> access(bool write) const {
        auto owner=owner_.lock();
        if (!owner || requirement_.type!=typeid(T) || (write && !requirement_.writable))
            throw std::runtime_error("Execution contract: expired/type/access resource view "+requirement_.id+" owner="+requirement_.owner);
        if (validate_) validate_();
        return std::shared_ptr<T>(std::move(owner),static_cast<T*>(object_));
    }
};

inline const char* toString(StateViewKind kind) {
    switch (kind) {
    case StateViewKind::Physical:return "Physical";
    case StateViewKind::Working:return "Working";
    case StateViewKind::Correction:return "Correction";
    case StateViewKind::OldTime:return "OldTime";
    case StateViewKind::Stage:return "Stage";
    case StateViewKind::Workspace:return "Workspace";
    }
    return "invalid";
}
} // namespace SF::System
