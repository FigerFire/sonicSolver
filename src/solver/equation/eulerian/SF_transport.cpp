/// @file SF_transport.cpp
/// @brief 双欧拉面输运系数、梯度和通量的结构网格算子。

#include "solver/equation/eulerian/SF_transport.h"
#include "SF_canonicalFace.h"

#include <cmath>
#include <stdexcept>

namespace SF::EulerianEulerian {

LinearAlgebra::SparseSystem assembleTransportEquation(
        const Field& field,
        const LinearAlgebra::DistributedRowMap& rowMap,
        const ScalarField& unknown,
        const ScalarField& previousConserved,
        const ScalarField& capacity,
        const ScalarField& diffusivity,
        const ScalarField& explicitSource,
        const std::vector<double>& canonicalMassFlux,
        double timeStep,
        bool axisymmetric,
        int radialCoordinate,
        ScalarField& diagonal,
        const ScalarField* implicitSink) {
    LinearAlgebra::SparseSystem matrix;
    matrix.firstRow = rowMap.first;
    matrix.lastRow = rowMap.last;
    matrix.globalSize = rowMap.total;
    matrix.rows.reserve(rowMap.localCells.size());
    matrix.rhs.reserve(rowMap.localCells.size());
    matrix.initialGuess.reserve(rowMap.localCells.size());

    for (int cell : rowMap.localCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        const double cellCapacity = capacity.values()[(size_t)cell];
        if (!std::isfinite(cellCapacity) || cellCapacity <= 0.0) {
            throw std::runtime_error(
                "Eulerian transport capacity must remain positive.");
        }

        LinearAlgebra::SparseRow row;
        row.globalRow = rowMap.global[(size_t)cell];
        const double inverseVolume = Ops::inverseCellVolume(
            field, i, j, k, axisymmetric, radialCoordinate);
        double diagonalValue = cellCapacity / timeStep;
        if (implicitSink) {
            const double sink = implicitSink->values()[(size_t)cell];
            if (!std::isfinite(sink) || sink < 0.0) {
                throw std::runtime_error(
                    "Eulerian implicit transport sink is invalid.");
            }
            diagonalValue += sink;
        }
        double rhs = previousConserved.values()[(size_t)cell] / timeStep
                   + explicitSource.values()[(size_t)cell];

        auto addColumn = [&](std::int64_t column, double value) {
            for (size_t entry = 0; entry < row.columns.size(); ++entry) {
                if (row.columns[entry] == column) {
                    row.values[entry] += value;
                    return;
                }
            }
            row.columns.push_back(column);
            row.values.push_back(value);
        };

        for (int axis = 0; axis < 3; ++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign : {-1, 1}) {
                int di = 0, dj = 0, dk = 0;
                Ops::offset(axis, sign, di, dj, dk);
                const int ni = i + di;
                const int nj = j + dj;
                const int nk = k + dk;
                const int neighbour = field.getIdx(ni, nj, nk);
                const double spacing = Ops::spacing(
                    field, i, j, k, ni, nj, nk);
                const int lowerI = sign > 0 ? i : ni;
                const int lowerJ = sign > 0 ? j : nj;
                const int lowerK = sign > 0 ? k : nk;
                const auto area = Ops::faceCofactor(
                    field, axis, lowerI, lowerJ, lowerK,
                    axisymmetric, radialCoordinate);
                const double areaNorm = std::sqrt(
                    area[0]*area[0] + area[1]*area[1]
                    + area[2]*area[2]);
                const double gamma = 0.5
                    * (diffusivity.values()[(size_t)cell]
                       + diffusivity.values()[(size_t)neighbour]);
                const double diffusion =
                    inverseVolume * gamma * areaNorm / spacing;
                if (!std::isfinite(diffusion) || diffusion < 0.0) {
                    throw std::runtime_error(
                        "Eulerian diffusion coefficient is invalid.");
                }
                diagonalValue += diffusion;
                const auto global = rowMap.global[(size_t)neighbour];
                if (global >= 0) {
                    addColumn(global, -diffusion);
                } else {
                    rhs += diffusion
                        * unknown.values()[(size_t)neighbour];
                }

                const size_t face = PhaseFaceFlux::index(
                    field, axis, lowerI, lowerJ, lowerK);
                const double outwardMassFlux =
                    Numerics::CanonicalFace::outward(
                        canonicalMassFlux[face],
                        sign>0
                            ?Numerics::CanonicalFace::CellSide::LowerOwner
                            :Numerics::CanonicalFace::CellSide::UpperNeighbour);
                const double convection =
                    inverseVolume * outwardMassFlux;
                if (!std::isfinite(convection)) {
                    throw std::runtime_error(
                        "Eulerian canonical mass flux is non-finite.");
                }
                if (convection >= 0.0) {
                    diagonalValue += convection;
                } else if (global >= 0) {
                    addColumn(global, convection);
                } else {
                    rhs -= convection
                        * unknown.values()[(size_t)neighbour];
                }
            }
        }

        diagonal.values()[(size_t)cell] = diagonalValue;
        addColumn(row.globalRow, diagonalValue);
        matrix.rows.push_back(std::move(row));
        matrix.rhs.push_back(rhs);
        matrix.initialGuess.push_back(
            unknown.values()[(size_t)cell]);
    }
    return matrix;
}

} // namespace SF::EulerianEulerian
