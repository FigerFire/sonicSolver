#pragma once

/// @file SF_sourceContribution.h
/// @brief Physical source models contribute terms to existing equations.

#include "core/config/SF_configTypes.h"

#include "core/system/SF_systemContribution.h"

namespace SF::Physics::SourceContribution {

/// @brief 源的单位属于具体目标等式，不能由全局主变量名字选择。
enum class Form { Acceleration, ForceDensity, EnergyDensity };
struct Target { std::string equation; Form form; };
std::vector<Target> targets(const System::EquationRegistry& equations);
void contribute(System::SystemContribution& system,const FDM::SourceConfig& config,
                const std::vector<Target>& targets);

} // namespace SF::Physics::SourceContribution
