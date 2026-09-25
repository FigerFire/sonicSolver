/// @file SF_hypreBackend.cpp
/// @brief HYPRE ParCSR 线性代数后端实现。
///
/// Ownership: HYPRE 是 linear-algebra backend，因此实现位于
/// solver/linearAlgebra/backend/hypre。MPI infrastructure 只提供
/// communicator/ownership；基础设施层不得反向依赖 solver。

#include "solver/linearAlgebra/hypre/SF_hypre.h"
#include "SF_schurPreconditioner.h"

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SF::LinearAlgebra {

bool hypreBackendLinked() { return true; }

namespace {

void requireHypre(int code, const char* operation) {
    if (code == 0) return;
    char description[256] = {0};
    HYPRE_DescribeError(code, description);
    HYPRE_ClearAllErrors();
    throw std::runtime_error(
        std::string("HYPRE operation failed: ") + operation
        + " (code=" + std::to_string(code) + ", " + description + ").");
}

void initializeHypre() {
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (!mpiInitialized) {
        throw std::runtime_error(
            "HYPRE linear solve requires an initialized MPI environment.");
    }
    static bool initialized = false;
    if (!initialized) {
        requireHypre(HYPRE_Initialize(), "HYPRE_Initialize");
        initialized = true;
    }
}

struct Objects {
    HYPRE_IJMatrix matrix = nullptr;
    HYPRE_IJVector rhs = nullptr;
    HYPRE_IJVector solution = nullptr;
    HYPRE_Solver solver = nullptr;
    HYPRE_Solver preconditioner = nullptr;
    bool pcg = false;
    bool ilu = false;
    std::unique_ptr<SchurPreconditioner> schur;

    ~Objects() {
        if (solver) {
            if (pcg) HYPRE_ParCSRPCGDestroy(solver);
            else HYPRE_ParCSRFlexGMRESDestroy(solver);
        }
        if (preconditioner) {
            if (ilu) HYPRE_ILUDestroy(preconditioner);
            else HYPRE_BoomerAMGDestroy(preconditioner);
        }
        if (matrix) HYPRE_IJMatrixDestroy(matrix);
        if (rhs) HYPRE_IJVectorDestroy(rhs);
        if (solution) HYPRE_IJVectorDestroy(solution);
    }
};

} // namespace

struct HypreSolverSession::Impl {
    explicit Impl(FDM::LinearSolverConfig value) : config(std::move(value)) {}

    FDM::LinearSolverConfig config;
    std::unique_ptr<Objects> objects;
    std::int64_t firstRow = -1;
    std::int64_t lastRow = -1;
    std::int64_t globalSize = -1;
    std::vector<std::vector<std::int64_t>> pattern;
    ReuseStatistics statistics;
    bool forceRebuild = true;

    bool samePattern(const SparseSystem& system) const {
        if (!objects || forceRebuild || firstRow != system.firstRow
            || lastRow != system.lastRow || globalSize != system.globalSize
            || pattern.size() != system.rows.size()) return false;
        for (size_t row = 0; row < pattern.size(); ++row) {
            if (pattern[row] != system.rows[row].columns) return false;
        }
        const int interval = config.structureRebuildInterval;
        return interval <= 0 || statistics.solves == 0
            || statistics.solves % interval != 0;
    }

    void configureSolver() {
        objects->ilu =
            config.preconditioner == FDM::LinearPreconditioner::ILU;
        if(config.preconditioner==FDM::LinearPreconditioner::KKTBlockSchur) {
            objects->schur=std::make_unique<SchurPreconditioner>(config);
        } else if(objects->ilu) {
            requireHypre(HYPRE_ILUCreate(&objects->preconditioner),
                         "ILUCreate");
            requireHypre(HYPRE_ILUSetPrintLevel(
                             objects->preconditioner,0),"ILUSetPrintLevel");
            requireHypre(HYPRE_ILUSetMaxIter(
                             objects->preconditioner,1),"ILUSetMaxIter");
            requireHypre(HYPRE_ILUSetTol(
                             objects->preconditioner,0.0),"ILUSetTol");
            requireHypre(HYPRE_ILUSetType(
                             objects->preconditioner,config.iluType),
                         "ILUSetType");
            requireHypre(HYPRE_ILUSetLevelOfFill(
                             objects->preconditioner,config.iluLevelOfFill),
                         "ILUSetLevelOfFill");
        } else {
            requireHypre(HYPRE_BoomerAMGCreate(&objects->preconditioner),
                         "BoomerAMGCreate");
            requireHypre(HYPRE_BoomerAMGSetPrintLevel(
                             objects->preconditioner, 0),
                         "BoomerAMGSetPrintLevel");
            requireHypre(HYPRE_BoomerAMGSetMaxIter(
                             objects->preconditioner, 1),
                         "BoomerAMGSetMaxIter");
            requireHypre(HYPRE_BoomerAMGSetTol(
                             objects->preconditioner, 0.0),
                         "BoomerAMGSetTol");
            configureAMGCycle(objects->preconditioner,config);
        }

        objects->pcg = config.method == FDM::KrylovMethod::PCG;
        if (objects->pcg) {
            requireHypre(HYPRE_ParCSRPCGCreate(
                             MPI_COMM_WORLD, &objects->solver),
                         "ParCSRPCGCreate");
            requireHypre(HYPRE_PCGSetTol(
                             objects->solver, config.relativeTolerance),
                         "PCGSetTol");
            requireHypre(HYPRE_PCGSetAbsoluteTol(
                             objects->solver, config.absoluteTolerance),
                         "PCGSetAbsoluteTol");
            requireHypre(HYPRE_PCGSetMaxIter(
                             objects->solver, config.maxIterations),
                         "PCGSetMaxIter");
            requireHypre(HYPRE_PCGSetPrecond(
                             objects->solver,
                             reinterpret_cast<HYPRE_PtrToSolverFcn>(
                                 HYPRE_BoomerAMGSolve),
                             reinterpret_cast<HYPRE_PtrToSolverFcn>(
                                 HYPRE_BoomerAMGSetup),
                             objects->preconditioner),
                         "PCGSetPrecond");
            return;
        }
        requireHypre(HYPRE_ParCSRFlexGMRESCreate(
                         MPI_COMM_WORLD, &objects->solver),
                     "ParCSRFlexGMRESCreate");
        requireHypre(HYPRE_FlexGMRESSetTol(
                         objects->solver, config.relativeTolerance),
                     "FlexGMRESSetTol");
        requireHypre(HYPRE_FlexGMRESSetAbsoluteTol(
                         objects->solver, config.absoluteTolerance),
                     "FlexGMRESSetAbsoluteTol");
        requireHypre(HYPRE_FlexGMRESSetMaxIter(
                         objects->solver, config.maxIterations),
                     "FlexGMRESSetMaxIter");
        requireHypre(HYPRE_FlexGMRESSetKDim(
                         objects->solver, config.krylovDimension),
                     "FlexGMRESSetKDim");
        HYPRE_PtrToSolverFcn precondition =
            reinterpret_cast<HYPRE_PtrToSolverFcn>(
                objects->ilu ? HYPRE_ILUSolve : HYPRE_BoomerAMGSolve);
        HYPRE_PtrToSolverFcn setup =
            reinterpret_cast<HYPRE_PtrToSolverFcn>(
                objects->ilu ? HYPRE_ILUSetup : HYPRE_BoomerAMGSetup);
        HYPRE_Solver preconditionContext=objects->preconditioner;
        if(objects->schur) {
            precondition=reinterpret_cast<HYPRE_PtrToSolverFcn>(SchurPreconditioner::apply);
            setup=reinterpret_cast<HYPRE_PtrToSolverFcn>(SchurPreconditioner::setup);
            preconditionContext=reinterpret_cast<HYPRE_Solver>(objects->schur.get());
        }
        requireHypre(HYPRE_FlexGMRESSetPrecond(
                         objects->solver,precondition,setup,
                         preconditionContext),
                     "FlexGMRESSetPrecond");
    }

    void rebuild(const SparseSystem& system) {
        objects = std::make_unique<Objects>();
        firstRow = system.firstRow;
        lastRow = system.lastRow;
        globalSize = system.globalSize;
        pattern.clear();
        pattern.reserve(system.rows.size());
        for (const auto& row : system.rows) pattern.push_back(row.columns);

        const HYPRE_BigInt first = static_cast<HYPRE_BigInt>(firstRow);
        const HYPRE_BigInt last = static_cast<HYPRE_BigInt>(lastRow);
        requireHypre(HYPRE_IJMatrixCreate(
                         MPI_COMM_WORLD, first, last, first, last,
                         &objects->matrix),
                     "IJMatrixCreate");
        requireHypre(HYPRE_IJMatrixSetObjectType(
                         objects->matrix, HYPRE_PARCSR),
                     "IJMatrixSetObjectType");
        requireHypre(HYPRE_IJVectorCreate(
                         MPI_COMM_WORLD, first, last, &objects->rhs),
                     "IJVectorCreate(rhs)");
        requireHypre(HYPRE_IJVectorSetObjectType(
                         objects->rhs, HYPRE_PARCSR),
                     "IJVectorSetObjectType(rhs)");
        requireHypre(HYPRE_IJVectorCreate(
                         MPI_COMM_WORLD, first, last, &objects->solution),
                     "IJVectorCreate(solution)");
        requireHypre(HYPRE_IJVectorSetObjectType(
                         objects->solution, HYPRE_PARCSR),
                     "IJVectorSetObjectType(solution)");
        configureSolver();
        forceRebuild = false;
        ++statistics.structureRebuilds;
    }

    void updateMatrix(const SparseSystem& system) {
        requireHypre(HYPRE_IJMatrixInitialize(objects->matrix),
                     "IJMatrixInitialize(update)");
        for (const SparseRow& source : system.rows) {
            HYPRE_Int count = static_cast<HYPRE_Int>(source.columns.size());
            std::vector<HYPRE_BigInt> columns(source.columns.begin(),
                                              source.columns.end());
            const HYPRE_BigInt row =
                static_cast<HYPRE_BigInt>(source.globalRow);
            requireHypre(HYPRE_IJMatrixSetValues(
                             objects->matrix, 1, &count, &row,
                             columns.data(), source.values.data()),
                         "IJMatrixSetValues(update)");
        }
        requireHypre(HYPRE_IJMatrixAssemble(objects->matrix),
                     "IJMatrixAssemble(update)");
        ++statistics.coefficientUpdates;
    }

    std::vector<HYPRE_BigInt> updateVectors(const SparseSystem& system) {
        const HYPRE_BigInt first =
            static_cast<HYPRE_BigInt>(system.firstRow);
        const size_t localSize = system.rhs.size();
        std::vector<HYPRE_BigInt> rows(localSize);
        std::vector<double> guess(localSize, 0.0);
        for (size_t n = 0; n < localSize; ++n) {
            rows[n] = first + static_cast<HYPRE_BigInt>(n);
            if (!system.initialGuess.empty()) guess[n] = system.initialGuess[n];
        }
        requireHypre(HYPRE_IJVectorInitialize(objects->rhs),
                     "IJVectorInitialize(rhs update)");
        requireHypre(HYPRE_IJVectorInitialize(objects->solution),
                     "IJVectorInitialize(solution update)");
        requireHypre(HYPRE_IJVectorSetValues(
                         objects->rhs, static_cast<HYPRE_Int>(localSize),
                         rows.data(), system.rhs.data()),
                     "IJVectorSetValues(rhs update)");
        requireHypre(HYPRE_IJVectorSetValues(
                         objects->solution, static_cast<HYPRE_Int>(localSize),
                         rows.data(), guess.data()),
                     "IJVectorSetValues(solution update)");
        requireHypre(HYPRE_IJVectorAssemble(objects->rhs),
                     "IJVectorAssemble(rhs update)");
        requireHypre(HYPRE_IJVectorAssemble(objects->solution),
                     "IJVectorAssemble(solution update)");
        ++statistics.rhsUpdates;
        return rows;
    }

    SolveResult run(const SparseSystem& system) {
        initializeHypre();
        const bool zeroRightHandSide = std::all_of(
            system.rhs.begin(), system.rhs.end(),
            [](double value) { return value == 0.0; });
        const bool zeroInitialGuess = system.initialGuess.empty()
            || std::all_of(
                system.initialGuess.begin(), system.initialGuess.end(),
                [](double value) { return value == 0.0; });
        // 分布式系统中零 RHS 只能在所有 rank 同时确认后短路。否则某个
        // 约束行 rank 会提前返回，而其它 rank 进入 HYPRE 集合调用并死锁。
        int globallyZero = (zeroRightHandSide && zeroInitialGuess) ? 1 : 0;
        MPI_Allreduce(
            MPI_IN_PLACE, &globallyZero, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (globallyZero != 0) {
            // A*0=0 是精确解。部分 HYPRE/AMG 组合会把这个情况报告为
            // 0 iterations、residual=1 且 converged=false；在后端调用前
            // 解析识别精确零系统，既不放宽容差，也不替换算法。
            SolveResult result;
            result.solution.assign(system.rhs.size(), 0.0);
            result.iterations = 0;
            result.relativeResidual = 0.0;
            result.converged = true;
            ++statistics.solves;
            return result;
        }
        const bool ownsWholeSystem =
            system.firstRow == 0
            && system.lastRow + 1 == system.globalSize;
        if (ownsWholeSystem) {
            const size_t count = system.rhs.size();
            std::vector<double> guess(count, 0.0);
            if (!system.initialGuess.empty()) guess = system.initialGuess;
            double residualSquared = 0.0;
            double rhsSquared = 0.0;
            for (size_t local = 0; local < count; ++local) {
                double product = 0.0;
                const SparseRow& row = system.rows[local];
                for (size_t entry = 0;
                     entry < row.columns.size(); ++entry) {
                    product += row.values[entry]
                        * guess[(size_t)row.columns[entry]];
                }
                const double residual = product - system.rhs[local];
                residualSquared += residual * residual;
                rhsSquared += system.rhs[local] * system.rhs[local];
            }
            const double residualNorm = std::sqrt(residualSquared);
            const double rhsNorm = std::sqrt(rhsSquared);
            const bool absoluteConverged =
                residualNorm <= config.absoluteTolerance;
            const bool relativeConverged = rhsNorm > 0.0
                && residualNorm / rhsNorm <= config.relativeTolerance;
            if (absoluteConverged || relativeConverged) {
                // 串行全局矩阵下可直接计算真实初始残差；满足用户配置
                // 容差时返回初值，避免HYPRE对极小非零RHS报告伪失败。
                SolveResult result;
                result.solution = std::move(guess);
                result.iterations = 0;
                result.relativeResidual =
                    rhsNorm > 0.0 ? residualNorm / rhsNorm : 0.0;
                result.converged = true;
                ++statistics.solves;
                return result;
            }
        }
        const bool rebuilt = !samePattern(system);
        if (rebuilt) rebuild(system);
        updateMatrix(system);
        const auto rows = updateVectors(system);

        HYPRE_ParCSRMatrix matrix = nullptr;
        HYPRE_ParVector rhs = nullptr;
        HYPRE_ParVector solution = nullptr;
        requireHypre(HYPRE_IJMatrixGetObject(
                         objects->matrix,
                         reinterpret_cast<void**>(&matrix)),
                     "IJMatrixGetObject");
        requireHypre(HYPRE_IJVectorGetObject(
                         objects->rhs, reinterpret_cast<void**>(&rhs)),
                     "IJVectorGetObject(rhs)");
        requireHypre(HYPRE_IJVectorGetObject(
                         objects->solution,
                         reinterpret_cast<void**>(&solution)),
                     "IJVectorGetObject(solution)");

        const int refreshInterval =
            config.preconditionerRefreshInterval;
        const bool refreshPreconditioner = rebuilt
            || (refreshInterval > 0 && statistics.solves > 0
                && statistics.solves % refreshInterval == 0);
        if(objects->schur && refreshPreconditioner)objects->schur->update(system);
        if (objects->pcg) {
            if (refreshPreconditioner) {
                requireHypre(HYPRE_ParCSRPCGSetup(
                                 objects->solver, matrix, rhs, solution),
                             "ParCSRPCGSetup(refresh)");
                ++statistics.preconditionerRefreshes;
            }
            requireHypre(HYPRE_ParCSRPCGSolve(
                             objects->solver, matrix, rhs, solution),
                         "ParCSRPCGSolve");
        } else {
            if (refreshPreconditioner) {
                requireHypre(HYPRE_ParCSRFlexGMRESSetup(
                                 objects->solver, matrix, rhs, solution),
                             "ParCSRFlexGMRESSetup(refresh)");
                ++statistics.preconditionerRefreshes;
            }
            const int solveCode = HYPRE_ParCSRFlexGMRESSolve(
                objects->solver, matrix, rhs, solution);
            if (solveCode != 0 && objects->schur
                && !objects->schur->error().empty()) {
                char description[256] = {0};
                HYPRE_DescribeError(solveCode, description);
                HYPRE_ClearAllErrors();
                throw std::runtime_error(
                    "HYPRE KKT FGMRES preconditioner failed: "
                    + objects->schur->error() + " (code="
                    + std::to_string(solveCode) + ", "
                    + description + ").");
            }
            requireHypre(solveCode, "ParCSRFlexGMRESSolve");
        }

        SolveResult result;
        HYPRE_Int iterations = 0;
        HYPRE_Int converged = 0;
        HYPRE_Real residual = 0.0;
        if (objects->pcg) {
            requireHypre(HYPRE_PCGGetNumIterations(
                             objects->solver, &iterations),
                         "PCGGetNumIterations");
            requireHypre(HYPRE_PCGGetFinalRelativeResidualNorm(
                             objects->solver, &residual),
                         "PCGGetFinalRelativeResidualNorm");
            requireHypre(HYPRE_PCGGetConverged(objects->solver, &converged),
                         "PCGGetConverged");
        } else {
            requireHypre(HYPRE_FlexGMRESGetNumIterations(
                             objects->solver, &iterations),
                         "FlexGMRESGetNumIterations");
            requireHypre(HYPRE_FlexGMRESGetFinalRelativeResidualNorm(
                             objects->solver, &residual),
                         "FlexGMRESGetFinalRelativeResidualNorm");
            requireHypre(HYPRE_FlexGMRESGetConverged(
                             objects->solver, &converged),
                         "FlexGMRESGetConverged");
        }
        result.iterations = iterations;
        result.relativeResidual = residual;
        // HYPRE 对零初始残差可能返回 0 iterations 且 converged flag 未置位；
        // 此时报告的最终相对残差已满足同一用户容差，属于严格收敛判定。
        result.converged = converged != 0
            || (std::isfinite(residual)
                && residual <= config.relativeTolerance);
        result.solution.resize(system.rhs.size(), 0.0);
        requireHypre(HYPRE_IJVectorGetValues(
                         objects->solution,
                         static_cast<HYPRE_Int>(rows.size()), rows.data(),
                         result.solution.data()),
                     "IJVectorGetValues(solution)");
        ++statistics.solves;
        if (!result.converged || !std::isfinite(result.relativeResidual)) {
            throw std::runtime_error(
                "HYPRE linear solver did not reach the configured tolerance; "
                "iterations=" + std::to_string(result.iterations)
                + ", relativeResidual="
                + std::to_string(result.relativeResidual) + ".");
        }
        return result;
    }
};

HypreSolverSession::HypreSolverSession(FDM::LinearSolverConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

HypreSolverSession::~HypreSolverSession() = default;
HypreSolverSession::HypreSolverSession(HypreSolverSession&&) noexcept = default;
HypreSolverSession& HypreSolverSession::operator=(
    HypreSolverSession&&) noexcept = default;

SolveResult HypreSolverSession::solve(const SparseSystem& system) {
    return impl_->run(system);
}

void HypreSolverSession::invalidateStructure() {
    impl_->forceRebuild = true;
}

const ReuseStatistics& HypreSolverSession::statistics() const {
    return impl_->statistics;
}

} // namespace SF::LinearAlgebra
