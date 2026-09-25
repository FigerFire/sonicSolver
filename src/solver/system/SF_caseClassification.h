#pragma once

/// @file SF_caseClassification.h
/// @brief EXPLAIN-ONLY 案例分类事实。
///
/// OWNERSHIP: 这些值是 compile 阶段一次性写定的说明性结论，只回答
/// "这个 case 被理解成了什么"，不参与任何 execution 决策。
///
/// 它们刻意从编译产物（ExecutableEquationSystem / CompiledNumericalSystem /
/// CompiledSolvePlan / RuntimeRequirements）中分离出来：runtime consumer
/// 不应当、也无法通过 Program 读到这些字符串。只有 `explain`/`check`
/// 这类 inspection 表面才读取它们。

#include "SF_configTypes.h"
#include "SF_physicsTemplate.h"

#include <string>

namespace SF::System {

struct CaseClassification {
    PhysicsTemplateKind templateOrigin = PhysicsTemplateKind::SingleFluid;
    std::string densityBehavior;
    std::string thermodynamicCompressibility;
};

} // namespace SF::System
