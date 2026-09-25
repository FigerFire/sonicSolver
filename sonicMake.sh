#!/usr/bin/env bash
set -euo pipefail

# sonicSolver build helper
#
# Usage:
#   ./sonicMake.sh          Clean, reconfigure, and build the whole project.
#   ./sonicMake.sh -solver  Configure for CLI only, then build sonicSolver.
#   ./sonicMake.sh -gui     Configure for GUI only, then build sonicGui.
#
# Notes:
# - The script uses CMake's cross-platform build command instead of raw make.
# - macOS uses clang++ by default. Linux uses the system/default C++ compiler,
#   unless the user provides CXX=... in the environment.
# - CLI builds install a sonicSolver symlink in ~/.local/bin by default.
#   They also install Bash/Zsh completion scripts under
#   ~/.local/share/sonicSolver/completions. Set SONIC_COMMAND_DIR or
#   SONIC_COMPLETION_DIR to customize the locations, or
#   SONIC_INSTALL_COMMAND=0 to disable both installs.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
GENERATOR="Unix Makefiles"
SONIC_COMMAND_INSTALLED=0
SONIC_COMPLETION_DIR_INSTALLED=""

usage() {
    sed -n '3,12p' "$0"
}

jobs_count() {
    local jobs=""
    if command -v sysctl >/dev/null 2>&1; then
        jobs="$(sysctl -n hw.ncpu 2>/dev/null || true)"
    fi
    if [ -z "$jobs" ] && command -v nproc >/dev/null 2>&1; then
        jobs="$(nproc 2>/dev/null || true)"
    fi
    if [ -z "$jobs" ] && command -v getconf >/dev/null 2>&1; then
        jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
    fi
    if [ -z "$jobs" ]; then
        jobs=1
    fi
    echo "$jobs"
}

is_macos() {
    [ "$(uname -s)" = "Darwin" ]
}

prepare_environment() {
    if is_macos; then
        # Avoid Xcode toolchain overrides that can make command-line CMake
        # select unexpected compilers.
        unset TOOLCHAINS
        unset METAL_DEVICE_COMMAND_LINE_TOOLS
    fi
}

ensure_compatible_build_dir() {
    if [ ! -f "${BUILD_DIR}/CMakeCache.txt" ]; then
        return
    fi
    if ! grep -q "CMAKE_GENERATOR:INTERNAL=${GENERATOR}" \
        "${BUILD_DIR}/CMakeCache.txt"; then
        echo "--- build 目录使用了其他 CMake generator，正在重建 build ---"
        rm -rf "$BUILD_DIR"
    fi
}

clean_for_full_build() {
    echo "--- 正在清理 build 目录并重新配置整个项目 ---"
    if is_macos && command -v dot_clean >/dev/null 2>&1; then
        # dot_clean is useful on macOS external drives, but it can be noisy
        # inside .git. Keep it non-fatal and quiet.
        dot_clean "$SCRIPT_DIR" >/dev/null 2>&1 || true
    fi
    rm -rf "$BUILD_DIR"
}

configure_project() {
    local build_cli="$1"
    local build_gui="$2"
    shift 2

    mkdir -p "$BUILD_DIR"
    local common_args=(-G "$GENERATOR")
    if is_macos && [ -z "${CXX:-}" ]; then
        common_args+=("-DCMAKE_CXX_COMPILER=clang++")
    fi

    echo "--- 正在配置 CMake: BUILD_CLI=${build_cli}, BUILD_GUI=${build_gui} ---"
    cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" \
        "${common_args[@]}" \
        -DBUILD_CLI="$build_cli" \
        -DBUILD_GUI="$build_gui" \
        "$@"
}

build_project() {
    local target="${1:-}"
    local jobs
    jobs="$(jobs_count)"

    if [ -n "$target" ]; then
        echo "--- 正在编译目标: ${target} (-j${jobs}) ---"
        cmake --build "$BUILD_DIR" --target "$target" -j "$jobs"
    else
        echo "--- 正在编译整个项目 (-j${jobs}) ---"
        cmake --build "$BUILD_DIR" -j "$jobs"
    fi
}

install_solver_command() {
    if [ "${SONIC_INSTALL_COMMAND:-1}" = "0" ]; then
        echo "--- 已按 SONIC_INSTALL_COMMAND=0 跳过终端命令安装 ---"
        return
    fi

    local command_dir="${SONIC_COMMAND_DIR:-${HOME}/.local/bin}"
    local command_path="${command_dir}/sonicSolver"
    local solver_path="${BUILD_DIR}/sonicSolver"

    if [ ! -x "$solver_path" ]; then
        echo "ERROR: 找不到已编译的 sonicSolver: $solver_path" >&2
        return 1
    fi
    if [ -e "$command_path" ] && [ ! -L "$command_path" ]; then
        echo "ERROR: 拒绝覆盖已有的非符号链接文件: $command_path" >&2
        echo "请设置 SONIC_COMMAND_DIR 指定其他安装目录。" >&2
        return 1
    fi

    mkdir -p "$command_dir"
    ln -sfn "$solver_path" "$command_path"
    SONIC_COMMAND_INSTALLED=1
    echo "--- 终端命令已安装: $command_path -> $solver_path ---"

    local completion_dir="${SONIC_COMPLETION_DIR:-${HOME}/.local/share/sonicSolver/completions}"
    mkdir -p "$completion_dir"
    cp "${SCRIPT_DIR}/scripts/sonicSolver-completion.bash" \
        "${completion_dir}/sonicSolver-completion.bash"
    cp "${SCRIPT_DIR}/scripts/_sonicSolver" \
        "${completion_dir}/_sonicSolver"
    SONIC_COMPLETION_DIR_INSTALLED="$completion_dir"
    echo "--- shell 自动补全已安装: $completion_dir ---"

    case ":${PATH}:" in
        *":${command_dir}:"*) ;;
        *)
            echo "提示：$command_dir 尚不在 PATH；请将其加入 shell 配置。"
            ;;
    esac
}

print_cli_commands() {
    local command="./build/sonicSolver"
    if [ "$SONIC_COMMAND_INSTALLED" -eq 1 ]; then
        command="sonicSolver"
    fi

    echo "可运行的 CLI 指令："
    echo "  ${command} --help"
    echo "  ${command} check [CASE]"
    echo "  ${command} explain [CASE]"
    echo "  ${command} doctor [CASE]"
    echo "  ${command} run [CASE]"
    echo "  ${command} models ibm --with mpi"
    echo "  ${command} why MODEL --with mpi"
    echo "  ${command} recipes"
    echo "  ${command} init [CASE] --recipe RECIPE"
    echo "  ${command} cleanCase [CASE]"
    echo "  ${command} postProcess [CASE]"
    echo "  ${command} -cleanResult [CASE]"
    echo "  ${command} -postProcessing [CASE]"

    if [ -n "$SONIC_COMPLETION_DIR_INSTALLED" ]; then
        echo "自动补全："
        echo "  Bash: source \"${SONIC_COMPLETION_DIR_INSTALLED}/sonicSolver-completion.bash\""
        echo "  Zsh:  fpath=(\"${SONIC_COMPLETION_DIR_INSTALLED}\" \$fpath); autoload -Uz compinit && compinit"
    fi
}

print_success() {
    local mode="$1"
    case "$mode" in
        all)
            echo "--- 编译完成：整个项目 ---"
            print_cli_commands
            ;;
        solver)
            echo "--- 编译完成：求解器 ---"
            print_cli_commands
            ;;
        gui)
            echo "--- 编译完成：GUI ---"
            if is_macos; then
                echo "可运行: ./build/sonicGui.app/Contents/MacOS/sonicGui"
            else
                echo "可运行: ./build/sonicGui"
            fi
            ;;
    esac
}

main() {
    local mode="all"
    if [ "$#" -gt 1 ]; then
        usage
        exit 2
    fi
    if [ "$#" -eq 1 ]; then
        case "$1" in
            -solver) mode="solver" ;;
            -gui) mode="gui" ;;
            -h|--help) usage; exit 0 ;;
            *)
                echo "未知参数: $1" >&2
                usage
                exit 2
                ;;
        esac
    fi

    prepare_environment
    cd "$SCRIPT_DIR"

    case "$mode" in
        all)
            clean_for_full_build
            configure_project ON ON
            build_project
            ;;
        solver)
            ensure_compatible_build_dir
            configure_project ON OFF
            build_project sonicSolver
            ;;
        gui)
            ensure_compatible_build_dir
            configure_project OFF ON
            build_project sonicGui
            ;;
    esac

    if [ "$mode" = "all" ] || [ "$mode" = "solver" ]; then
        install_solver_command
    fi

    print_success "$mode"
}

main "$@"
