#pragma once

/// @file SF_fields.h
/// @brief 暴露 VTK piece layout 与 scalar field composition helpers。
///
/// Data flow:
///   mesh / phase / interface / turbulence state
///       -> non-owning ResultWriter field descriptors
///       -> output writer
///
/// 本文件不写文件、不修改 physical state，也不参与 runtime scheduling。

#include "SF_resultWriter.h"

#include <vector>

namespace SF {
class MultiBlockMesh;
namespace State { class VariableRegistry; }
namespace Turbulence { class EquationSystem; }
namespace Physics::InterfaceModels { class Model; }
namespace Physics::Multiphase { class MultiPhaseModel; }
namespace Physics::PhaseSystems { class PhaseSystem; }
namespace Application::Output {

std::vector<ResultWriter::PieceLayout> buildMultiBlockVTKPieces(
    const MultiBlockMesh& mesh);
bool needsMultiFieldMPI(const MultiBlockMesh& mesh, int mpiSize);
std::vector<ResultWriter::ScalarField> multiPhaseVTKScalars(
    Physics::Multiphase::MultiPhaseModel& model,
    const State::VariableRegistry* registry = nullptr);
std::vector<ResultWriter::ScalarField> interfaceVTKScalars(
    Physics::InterfaceModels::Model& model,
    const State::VariableRegistry* registry = nullptr);
std::vector<ResultWriter::ScalarField> eulerianVTKScalars(
    Physics::PhaseSystems::PhaseSystem& system);
std::vector<ResultWriter::ScalarField> eulerianTurbulenceVTKScalars(
    const Turbulence::EquationSystem& turbulence);

} // namespace Application::Output
} // namespace SF
