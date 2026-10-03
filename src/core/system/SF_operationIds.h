#pragma once

/// @file SF_operationIds.h
/// @brief Stable built-in operation identities; custom OpIds remain open strings.

namespace SF::System::OpIds {

inline constexpr const char* PressurePrepare = "pressure.prepare";
inline constexpr const char* PressureStepBegin = "pressure.step.begin";
inline constexpr const char* PressureIterationBegin = "pressure.iteration.begin";
inline constexpr const char* MomentumAssemble = "momentum.assemble";
inline constexpr const char* MomentumSolve = "momentum.solve";
inline constexpr const char* PressureBoundaryPrepare = "pressure.boundary.prepare";
inline constexpr const char* PressureAssemble = "pressure.assemble";
inline constexpr const char* PressureSolve = "pressure.solve";
inline constexpr const char* PressureUpdatePrepare = "pressure.update.prepare";
inline constexpr const char* VelocityCorrect = "velocity.correct";
inline constexpr const char* FluxCorrect = "flux.correct";
inline constexpr const char* PressureCorrectionCommit = "pressure.correction.commit";
inline constexpr const char* PressureRelaxationApply = "pressure.relaxation.apply";
inline constexpr const char* PressureFluxConsistencyRestore = "pressure.flux.consistency.restore";
inline constexpr const char* PressureConvergenceEvaluate = "pressure.convergence.evaluate";
inline constexpr const char* PressureIterationEnd = "pressure.iteration.end";
inline constexpr const char* PressureStepCommit = "pressure.step.commit";
inline constexpr const char* MomentumPredictor = "momentum.predictor";
inline constexpr const char* PressureCorrection = "pressure.correction";
inline constexpr const char* IbmConstraintProject = "ibm.constraint.project";
inline constexpr const char* IbmKktSolve = "ibm.kkt.solve";
inline constexpr const char* TimeCommit = "time.commit";
inline constexpr const char* FlowStepPrepare = "flow.step.prepare";
inline constexpr const char* FlowDtCompute = "flow.dt.compute";
inline constexpr const char* FlowStepBegin = "flow.step.begin";
inline constexpr const char* TurbulenceAdvance = "turbulence.advance";
inline constexpr const char* TurbulenceClosureRefresh = "turbulence.closure.refresh";
inline constexpr const char* FlowStepCommit = "flow.step.commit";
inline constexpr const char* ExplicitStageExecute = "explicit.stage.execute";
inline constexpr const char* EeDtCompute = "ee.dt.compute";
inline constexpr const char* EeStepBegin = "ee.step.begin";
inline constexpr const char* EeInterphaseCompute = "ee.interphase.compute";
inline constexpr const char* EeSourcesAssemble = "ee.sources.assemble";
inline constexpr const char* EeTurbulencePrepare = "ee.turbulence.prepare";
inline constexpr const char* EeSourcesValidate = "ee.sources.validate";
inline constexpr const char* EeMomentumDiagonal = "ee.momentum.diagonal";
inline constexpr const char* EeMomentumFlux = "ee.momentum.flux";
inline constexpr const char* EeFaceFluxCanonical = "ee.faceFlux.canonical";
inline constexpr const char* EeContinuityAssemble = "ee.continuity.assemble";
inline constexpr const char* EeBoundaryPrepare = "ee.boundary.prepare";
inline constexpr const char* EeMomentumSolve = "ee.momentum.solve";
inline constexpr const char* EeInterphaseCorrect = "ee.interphase.correct";
inline constexpr const char* EeBoundaryAfterMomentum = "ee.boundary.afterMomentum";
inline constexpr const char* EeDiagonalSync = "ee.diagonal.sync";
inline constexpr const char* EeMomentumFluxAfter = "ee.momentum.flux.after";
inline constexpr const char* EeFaceFluxCanonicalAfter = "ee.faceFlux.canonical.after";
inline constexpr const char* EePressureSolve = "ee.pressure.solve";
inline constexpr const char* EePressurePublish = "ee.pressure.publish";
inline constexpr const char* EePressureSync = "ee.pressure.sync";
inline constexpr const char* EePhaseCorrect = "ee.phase.correct";
inline constexpr const char* EeFaceFluxCorrect = "ee.faceFlux.correct";
inline constexpr const char* EeBoundaryAfterPressure = "ee.boundary.afterPressure";
inline constexpr const char* EeFaceFluxCanonicalPressure = "ee.faceFlux.canonical.pressure";
inline constexpr const char* EeEnergySolve = "ee.energy.solve";
inline constexpr const char* EeTurbulenceSolve = "ee.turbulence.solve";
inline constexpr const char* EeBoundaryFinal = "ee.boundary.final";
inline constexpr const char* EeOuterValidate = "ee.outer.validate";
inline constexpr const char* EeStepCommit = "ee.step.commit";
inline constexpr const char* EeTimeCommit = "ee.time.commit";

} // namespace SF::System::OpIds
