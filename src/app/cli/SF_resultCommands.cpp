/// @file SF_resultCommands.cpp
/// @brief 实现与求解流程解耦的 result 清理和 ParaView 启动命令。

#include "SF_resultCommands.h"

#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <spawn.h>
#include <stdexcept>
#include <string>
#include <system_error>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

extern char** environ;

namespace {

std::filesystem::path resultDirectory(const std::filesystem::path& caseDirectory)
{
    const auto workingDirectory = std::filesystem::current_path();
    const auto casePath = caseDirectory.empty()
        ? workingDirectory
        : (caseDirectory.is_absolute() ? caseDirectory : workingDirectory / caseDirectory);
    return (casePath / "result").lexically_normal();
}

std::vector<std::filesystem::path> pvdFiles(const std::filesystem::path& resultPath)
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(resultPath)) {
        const std::string name = entry.path().filename().string();
        // macOS 外置卷可能产生 `._name.pvd` AppleDouble 文件；隐藏文件不是
        // 求解器结果索引，不能参与唯一 PVD 判定。
        if (!name.empty() && name.front() != '.'
            && entry.is_regular_file() && entry.path().extension() == ".pvd") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

#ifdef __APPLE__
std::filesystem::path paraViewApplication()
{
    const char* environmentPath = std::getenv("PATH");
    const std::string searchPath = environmentPath ? environmentPath : "";
    std::size_t begin = 0;
    do {
        const auto end = searchPath.find(':', begin);
        const auto directory = searchPath.substr(begin, end - begin);
        const auto executable = std::filesystem::path(directory.empty() ? "." : directory)
            / "paraview";
        if (access(executable.c_str(), X_OK) == 0) {
            auto path = std::filesystem::canonical(executable);
            while (path != path.root_path()) {
                if (path.extension() == ".app") return path;
                path = path.parent_path();
            }
            throw std::runtime_error("PATH 中的 paraview 不属于 macOS .app: "
                + executable.string());
        }
        if (end == std::string::npos) break;
        begin = end + 1;
    } while (true);
    throw std::runtime_error("无法找到 paraview，请确认它已安装并位于 PATH");
}
#endif

int spawnViewer(const std::filesystem::path& pvdPath)
{
    const std::string file = pvdPath.string();
    pid_t processId = 0;

    std::vector<char> fileArgument(file.begin(), file.end());
    fileArgument.push_back('\0');
#ifdef __APPLE__
    // LaunchServices owns the GUI lifetime. Closing a macOS window need not
    // terminate the app, so the CLI waits only for the launch request, not -W.
    char command[] = "/usr/bin/open";
    char applicationOption[] = "-a";
    // Resolve the selected PATH executable to its bundle, including versioned
    // app names; a bundle identifier may not yet be registered with the OS.
    std::string application = paraViewApplication().string();
    char* arguments[] = {command, applicationOption, application.data(), fileArgument.data(), nullptr};
#else
    char command[] = "paraview";
    char* arguments[] = {command, fileArgument.data(), nullptr};
#endif
    const int status = posix_spawnp(&processId, command, nullptr, nullptr, arguments, environ);
    if (status != 0) {
        SF::broadcast("Fatal: 无法启动 ParaView: ", std::strerror(status));
        return 1;
    }

    // Reap the process we created and retain launch/viewer failure diagnostics.
    int exitStatus = 0;
    pid_t waited;
    do {
        waited = waitpid(processId, &exitStatus, 0);
    } while (waited == -1 && errno == EINTR);
    if (waited == -1) {
        SF::broadcast("Fatal: 无法获取 ParaView 启动进程状态: ", std::strerror(errno));
        return 1;
    }
    if (WIFSIGNALED(exitStatus)) {
        SF::broadcast("Fatal: ParaView 启动进程被信号终止: ", WTERMSIG(exitStatus));
        return 128 + WTERMSIG(exitStatus);
    }
    if (!WIFEXITED(exitStatus) || WEXITSTATUS(exitStatus) != 0) {
        const int code = WIFEXITED(exitStatus) ? WEXITSTATUS(exitStatus) : 1;
        SF::broadcast("Fatal: ParaView 启动进程返回错误: ", code);
        return code;
    }

    SF::broadcast("ParaView: ", pvdPath.string());
    return 0;
}

} // namespace

int SF::CLI::cleanResult(const std::filesystem::path& caseDirectory)
{
    const auto resultPath = resultDirectory(caseDirectory);
    std::error_code error;

    if (!std::filesystem::exists(resultPath, error)) {
        if (error) {
            SF::broadcast("Fatal: ", "无法检查 result 目录 " + resultPath.string()
                + ": " + error.message());
            return -1;
        }
        SF::broadcast("Result clean: 目录不存在，无需清理: ", resultPath.string());
        return 0;
    }
    if (!std::filesystem::is_directory(resultPath, error) || error) {
        SF::broadcast("Fatal: result 路径不是目录: ", resultPath.string());
        return -1;
    }
    const auto resultStatus = std::filesystem::symlink_status(resultPath, error);
    if (error) {
        SF::broadcast("Fatal: ", "无法检查 result 目录类型 "
            + resultPath.string() + ": " + error.message());
        return -1;
    }
    if (std::filesystem::is_symlink(resultStatus)) {
        SF::broadcast("Fatal: 拒绝清理符号链接 result 目录: ", resultPath.string());
        return -1;
    }

    std::uintmax_t removed = 0;
    for (const auto& entry : std::filesystem::directory_iterator(resultPath)) {
        removed += std::filesystem::remove_all(entry.path(), error);
        if (error) {
            SF::broadcast("Fatal: ", "清理失败 " + entry.path().string()
                + ": " + error.message());
            return -1;
        }
    }

    SF::broadcast("Result clean: ", resultPath.string()
        + " (removed entries: " + std::to_string(removed) + ")");
    return 0;
}

int SF::CLI::openPostProcessing(const std::filesystem::path& caseDirectory)
{
    const auto resultPath = resultDirectory(caseDirectory);
    if (!std::filesystem::is_directory(resultPath)) {
        SF::broadcast("Fatal: result 目录不存在: ", resultPath.string());
        return -1;
    }

    const auto files = pvdFiles(resultPath);
    if (files.empty()) {
        SF::broadcast("Fatal: result 目录中没有 PVD 文件: ", resultPath.string());
        return -1;
    }
    if (files.size() != 1) {
        SF::broadcast("Fatal: result 目录中存在多个 PVD 文件，请只保留一个后重试: ",
                      resultPath.string());
        for (const auto& file : files) {
            SF::broadcast("  ", file.filename().string());
        }
        return -1;
    }

    return spawnViewer(files.front());
}
