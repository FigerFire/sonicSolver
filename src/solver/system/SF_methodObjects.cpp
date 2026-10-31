/// @file SF_methodObjects.cpp
/// @brief 编译期数值实现注册；不拥有 HOW 编译或执行。
#include "SF_methodObjects.h"
#include <algorithm>
namespace SF::System {
void TemporalMethodRegistry::add(const ITemporalMethod& method) {
    if (std::any_of(methods_.begin(),methods_.end(),
        [&](const ITemporalMethod* existing) { return existing->id()==method.id(); }))
        throw std::runtime_error("Duplicate TemporalMethod registration.");
    methods_.push_back(&method);
}

const ITemporalMethod& TemporalMethodRegistry::at(FDM::TimeRecipeId id) const {
    for (const auto* method:methods_) if (method->id()==id) return *method;
    throw std::runtime_error("Selected TemporalMethod has no registered implementation.");
}

void ProviderRegistry::add(const IProvider& method) {
    if (std::any_of(methods_.begin(),methods_.end(),
        [&](const IProvider* existing) { return existing->id()==method.id(); }))
        throw std::runtime_error("Duplicate EquationMethod registration.");
    methods_.push_back(&method);
}

const IProvider& ProviderRegistry::at(std::string_view id) const {
    for (const auto* method:methods_) if (method->id()==id) return *method;
    throw std::runtime_error("EquationMethod '"+std::string(id)
        +"' has no registered implementation.");
}


} // namespace SF::System
