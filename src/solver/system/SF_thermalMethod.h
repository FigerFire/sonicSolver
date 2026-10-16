#pragma once
#include "SF_methodObjects.h"
namespace SF::System {
/// @brief 现有标量 Upwind/Central 内核的独立温度离散契约；仅 Forward Euler。
void addThermalMethod(ProviderRegistry&);
inline constexpr const char* ThermalProvider="equation.thermal";
}
