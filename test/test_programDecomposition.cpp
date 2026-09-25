/// @file test_programDecomposition.cpp
/// @brief §18/§20/§35 回归：Program 是纯视图，且不携带任何运行期状态 authority。
///
/// 断言刻意做在"成员数量 + 成员类型 + 对象布局"这一层，因此无法靠放宽断言绕过：
///   - `Program` 必须恰好是四成员聚合，成员依次是四类编译产物；
///   - `Program` 只持有引用，sizeof 等于四个指针，不存在任何 owned state。

#include "solver/system/SF_resolvedSimulationSystem.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <type_traits>

namespace {

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "program decomposition test failed: " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

/// 结构化绑定要求成员数量精确匹配：新增第五个成员会直接编译失败，
/// 从而阻止 Program 重新变成携带运行期状态的对象。
void requireFourCompileProducts(const SF::System::Program& program) {
    const auto& [equations, numerics, solve, runtime] = program;
    using E = std::remove_cv_t<std::remove_reference_t<decltype(equations)>>;
    using N = std::remove_cv_t<std::remove_reference_t<decltype(numerics)>>;
    using S = std::remove_cv_t<std::remove_reference_t<decltype(solve)>>;
    using R = std::remove_cv_t<std::remove_reference_t<decltype(runtime)>>;
    static_assert(std::is_same_v<E,SF::System::ExecutableEquationSystem>,
        "Program member 0 must be the executable equation system");
    static_assert(std::is_same_v<N,SF::System::CompiledNumericalSystem>,
        "Program member 1 must be the compiled numerical system");
    static_assert(std::is_same_v<S,SF::System::CompiledSolvePlan>,
        "Program member 2 must be the compiled solve plan");
    static_assert(std::is_same_v<R,SF::System::RuntimeRequirements>,
        "Program member 3 must be the runtime requirements");
}

} // namespace

int main() {
    using SF::System::Program;

    // 纯引用视图：没有任何 owned state，也没有第二份 authority。
    static_assert(sizeof(Program) == 4 * sizeof(const void*),
                  "Program owns state: sizeof(Program) != four references");
    static_assert(std::is_trivially_copyable_v<Program>,
                  "Program stopped being a cheap view type");

    SF::System::ResolvedSimulationSystem resolved;
    requireFourCompileProducts(SF::System::programOf(resolved));

    const Program program = SF::System::programOf(resolved);
    require(&program.equations == &resolved.executableSystem,
            "Program copied the executable equation system");
    require(&program.numerics == &resolved.numericalSystem,
            "Program copied the compiled numerical system");
    require(&program.solve == &resolved.solvePlan,
            "Program copied the compiled solve plan");
    require(&program.runtime == &resolved.runtime,
            "Program copied the runtime requirements");

    const SF::System::CompilationResult compilation =
        SF::System::compilationOf(resolved);
    require(&compilation.raw == &resolved.rawSystem,
            "CompilationResult copied the raw equation system");
    require(&compilation.transformations == &resolved.transformations,
            "CompilationResult copied the transformation history");
    require(&compilation.program.equations == &resolved.executableSystem,
            "CompilationResult holds a detached program view");

    std::cout << "Program holds exactly the four compile products by reference\n";
    return 0;
}
