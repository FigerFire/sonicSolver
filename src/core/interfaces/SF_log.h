#pragma once

/// @file SF_log.h
/// @brief 与求解配置无关的轻量日志接口。

#include <cstdlib>
#include <iostream>

namespace SF {

/// @brief 当前进程是否允许输出普通日志。
inline bool logOutputEnabled = true;

/// @brief 启用或关闭当前进程的普通日志输出。
inline void setLogOutputEnabled(bool enabled) {
    logOutputEnabled = enabled;
}

/// @brief 判断当前 MPI 进程是否负责普通日志输出。
inline bool logFromThisProcess() {
    if (!logOutputEnabled) return false;

    const char* rankVariables[] = {
        "OMPI_COMM_WORLD_RANK",
        "PMI_RANK",
        "PMIX_RANK",
        "MV2_COMM_WORLD_RANK",
        "SLURM_PROCID"
    };
    for (const char* name : rankVariables) {
        const char* value = std::getenv(name);
        if (value != nullptr && value[0] != '\0') {
            return std::atoi(value) == 0;
        }
    }
    return true;
}

/// @brief 输出一条带 sonicSolver 前缀的日志。
template<typename T>
inline void broadcast(const char* key, const T& value) {
    if (!logFromThisProcess()) return;
    std::cout << "[SF] " << key << value << std::endl;
}

} // namespace SF
