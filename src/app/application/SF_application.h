#pragma once

/// @file SF_application.h
/// @brief case 加载、模块装配、运行与输出的应用入口。

#include <string>

namespace SF::Application {

/// @brief 已解析的命令行运行请求。
///
/// CLI 负责解析命令行并填充此结构；application 层不再感知原始 argv。
/// `argc` / `argv` 仅作为非拥有透传参数，供 MPI 初始化（MPI_Init）使用，
/// 其生命周期必须覆盖整个 run() 调用。
struct RunRequest {
    std::string casePath{"."};
    bool initialOutputOnly{false};
    int stepLimit{0};

    int argc{0};
    char** argv{nullptr};
};

/// @brief 运行一个由 native SonicFile case 驱动的任务。
///
/// 旧 OpenFOAM-like case 仍由 application compatibility adapter 读取；CLI
/// 的 `run`、默认 `sonicSolver CASE` 和 legacy flags 最终都进入此入口。
/// @return 进程兼容的退出码。
int run(const RunRequest& request);

} // namespace SF::Application
