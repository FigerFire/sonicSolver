/// @file SF_resultCommands.h
/// @brief result 目录清理与 ParaView 后处理命令接口。

#pragma once

#include <filesystem>

namespace SF::CLI {

/// @brief 清空指定 case 的 result 目录，但保留 result 目录本身。
/// @param caseDirectory case 目录；空路径表示当前工作目录。
/// @return 成功返回 0，路径或文件系统操作失败返回非零值。
int cleanResult(const std::filesystem::path& caseDirectory);

/// @brief 查找指定 case/result 中唯一的 PVD 文件并用 ParaView 打开。
/// @param caseDirectory case 目录；空路径表示当前工作目录。
/// @return 成功启动 ParaView 返回 0，否则返回非零值。
int openPostProcessing(const std::filesystem::path& caseDirectory);

} // namespace SF::CLI
