#pragma once

/// @file SF_ilwPlan.h
/// @brief IBM-ILW ghost closure 的预计算数据；不属于通用 Field 状态定义。

#include <array>
#include <vector>

namespace SF {

struct IBMILWFitCandidate {
    int derivativeOrder = 0;
    int polynomialOrder = 0;
    double linearWeight = 0.0;
    double localSpacing = 1.0;
    std::vector<int> sampleCells;
    std::vector<std::vector<double>> basisRows;
    std::vector<std::vector<double>> projectionRows;
    std::vector<double> derivativeWeights;
};

struct IBMILWPointPlan {
    bool valid = false;
    bool useT2 = false;
    int maxOrder = 0;
    std::array<double, 3> wallPoint{0.0, 0.0, 0.0};
    std::array<double, 3> normal{1.0, 0.0, 0.0};
    std::array<double, 3> tangent1{0.0, 1.0, 0.0};
    std::array<double, 3> tangent2{0.0, 0.0, 1.0};
    std::array<double, 3> wallVelocity{0.0, 0.0, 0.0};
    double targetDistance = 0.0;
    double curvatureK11 = 0.0;
    double curvatureK12 = 0.0;
    double curvatureK22 = 0.0;
    IBMILWFitCandidate wallFit;
    std::vector<IBMILWFitCandidate> higherFits;
};

} // namespace SF
