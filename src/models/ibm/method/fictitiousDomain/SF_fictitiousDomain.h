#pragma once

/// @file SF_fictitiousDomain.h
/// @brief Fictitious-domain 方法族入口。
///
/// 包含显式自推进、FTS 自推进/给定速度以及全隐式给定速度/自推进分支。
/// FTS 与 fully implicit 是不同时间耦合，不提供含混的“implicit FTS”别名。

#include "method/SF_method.h"
