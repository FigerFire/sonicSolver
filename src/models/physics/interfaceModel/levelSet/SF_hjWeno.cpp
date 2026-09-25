/// @file SF_hjWeno.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_hjWeno.h"

#include "core/mesh/SF_dimension.h"

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

std::string cellText(int i, int j, int k) {
    std::ostringstream os;
    os << "(" << i << "," << j << "," << k << ")";
    return os.str();
}

int axisLength(const Field& field, int axis) {
    if (axis == 0) return field.NX();
    if (axis == 1) return field.NY();
    return field.NZ();
}

int offsetI(int axis, int offset) { return axis == 0 ? offset : 0; }
int offsetJ(int axis, int offset) { return axis == 1 ? offset : 0; }
int offsetK(int axis, int offset) { return axis == 2 ? offset : 0; }

double phiAt(const LevelSetField& levelSet,
             const std::vector<double>& values,
             int i, int j, int k,
             int axis,
             int offset) {
    const int ii = i + offsetI(axis, offset);
    const int jj = j + offsetJ(axis, offset);
    const int kk = k + offsetK(axis, offset);
    const int id = levelSet.getIdx(ii, jj, kk);
    return values[(size_t)id];
}

double slopeAt(const LevelSetField& levelSet,
               const std::vector<double>& values,
               int i, int j, int k,
               int axis,
               int offset) {
    return phiAt(levelSet, values, i, j, k, axis, offset + 1)
         - phiAt(levelSet, values, i, j, k, axis, offset);
}

void validateWeightOptions(const HJWenoWeightOptions& options) {
    if (!std::isfinite(options.epsilon) || options.epsilon <= 0.0
        || !std::isfinite(options.power) || options.power <= 0.0) {
        throw std::runtime_error(
            "HJWeno weights require user-provided finite epsilon > 0 and power > 0.");
    }
}

double nonlinearWeight(double linearWeight,
                       double smoothness,
                       const HJWenoWeightOptions& options) {
    return linearWeight
        / std::pow(options.epsilon + smoothness, options.power);
}

double weno3Left(const double* v, const HJWenoWeightOptions& options) {
    const double p0 = -0.5 * v[0] + 1.5 * v[1];
    const double p1 =  0.5 * v[1] + 0.5 * v[2];
    const double b0 = (v[1] - v[0]) * (v[1] - v[0]);
    const double b1 = (v[2] - v[1]) * (v[2] - v[1]);
    const double a0 = nonlinearWeight(1.0 / 3.0, b0, options);
    const double a1 = nonlinearWeight(2.0 / 3.0, b1, options);
    const double sum = a0 + a1;
    return (a0 * p0 + a1 * p1) / sum;
}

double weno5Left(const double* v, const HJWenoWeightOptions& options) {
    const double p0 = ( 2.0 * v[0] - 7.0 * v[1] + 11.0 * v[2]) / 6.0;
    const double p1 = (-1.0 * v[1] + 5.0 * v[2] +  2.0 * v[3]) / 6.0;
    const double p2 = ( 2.0 * v[2] + 5.0 * v[3] -  1.0 * v[4]) / 6.0;

    const double b0 =
        13.0 / 12.0 * std::pow(v[0] - 2.0 * v[1] + v[2], 2)
        + 0.25 * std::pow(v[0] - 4.0 * v[1] + 3.0 * v[2], 2);
    const double b1 =
        13.0 / 12.0 * std::pow(v[1] - 2.0 * v[2] + v[3], 2)
        + 0.25 * std::pow(v[1] - v[3], 2);
    const double b2 =
        13.0 / 12.0 * std::pow(v[2] - 2.0 * v[3] + v[4], 2)
        + 0.25 * std::pow(3.0 * v[2] - 4.0 * v[3] + v[4], 2);

    const double a0 = nonlinearWeight(0.1, b0, options);
    const double a1 = nonlinearWeight(0.6, b1, options);
    const double a2 = nonlinearWeight(0.3, b2, options);
    const double sum = a0 + a1 + a2;
    return (a0 * p0 + a1 * p1 + a2 * p2) / sum;
}

double beta4(const double v[4], double a, double b) {
    const double d1 = v[1] - v[0];
    const double d2 = 0.5 * (v[2] - 2.0 * v[1] + v[0]);
    const double d3 = (v[3] - 3.0 * v[2] + 3.0 * v[1] - v[0]) / 6.0;

    const double A = d1 - d2 + 2.0 * d3;
    const double B = 2.0 * d2 - 6.0 * d3;
    const double C = 3.0 * d3;

    auto integral = [&](auto f) { return f(b) - f(a); };
    const double beta1 =
        A * A * integral([](double x) { return x; })
        + A * B * integral([](double x) { return x * x; })
        + (B * B + 2.0 * A * C)
              * integral([](double x) { return x * x * x / 3.0; })
        + B * C * integral([](double x) {
              const double x2 = x * x;
              return x2 * x2 / 2.0;
          })
        + C * C * integral([](double x) {
              const double x2 = x * x;
              return x2 * x2 * x / 5.0;
          });

    const double D = 2.0 * d2 - 6.0 * d3;
    const double E = 6.0 * d3;
    const double beta2 =
        D * D * integral([](double x) { return x; })
        + D * E * integral([](double x) { return x * x; })
        + E * E * integral([](double x) { return x * x * x / 3.0; });
    const double beta3 = 36.0 * d3 * d3 * integral([](double x) { return x; });
    return beta1 + beta2 + beta3;
}

double weno7Left(const double* v, const HJWenoWeightOptions& options) {
    const double p0 =
        (-3.0 * v[0] + 13.0 * v[1] - 23.0 * v[2] + 25.0 * v[3]) / 12.0;
    const double p1 =
        (        v[1] -  5.0 * v[2] + 13.0 * v[3] +  3.0 * v[4]) / 12.0;
    const double p2 =
        (      - v[2] +  7.0 * v[3] +  7.0 * v[4] -        v[5]) / 12.0;
    const double p3 =
        ( 3.0 * v[3] + 13.0 * v[4] -  5.0 * v[5] +        v[6]) / 12.0;

    double s0[4] = {v[0], v[1], v[2], v[3]};
    double s1[4] = {v[1], v[2], v[3], v[4]};
    double s2[4] = {v[2], v[3], v[4], v[5]};
    double s3[4] = {v[3], v[4], v[5], v[6]};
    const double b0 = beta4(s0, 2.5, 3.5);
    const double b1 = beta4(s1, 1.5, 2.5);
    const double b2 = beta4(s2, 0.5, 1.5);
    const double b3 = beta4(s3, -0.5, 0.5);

    const double a0 = nonlinearWeight(1.0 / 35.0, b0, options);
    const double a1 = nonlinearWeight(12.0 / 35.0, b1, options);
    const double a2 = nonlinearWeight(18.0 / 35.0, b2, options);
    const double a3 = nonlinearWeight(4.0 / 35.0, b3, options);
    const double sum = a0 + a1 + a2 + a3;
    return (a0 * p0 + a1 * p1 + a2 * p2 + a3 * p3) / sum;
}

double reconstructLeft(const double* values,
                       int order,
                       const HJWenoWeightOptions& options) {
    if (order == 3) return weno3Left(values, options);
    if (order == 5) return weno5Left(values, options);
    return weno7Left(values, options);
}

} // namespace

void HJWeno::validateOrder(int order, const char* context) {
    if (order == 3 || order == 5 || order == 7) return;
    throw std::runtime_error(
        std::string(context)
        + ": Level Set HJ-WENO order must be 3, 5, or 7; got "
        + std::to_string(order) + ".");
}

int HJWeno::requiredGhostLayers(int order) {
    validateOrder(order, "HJWeno::requiredGhostLayers");
    if (order == 3) return 2;
    if (order == 5) return 3;
    return 4;
}

void HJWeno::requireStencil(const Field& field,
                            int order,
                            const char* context) {
    validateOrder(order, context);
    const int required = requiredGhostLayers(order);
    if (field.NG() < required
        && (activeAxis(field, 0) || activeAxis(field, 1) || activeAxis(field, 2))) {
        throw std::runtime_error(
            std::string(context) + ": WENO" + std::to_string(order)
            + " requires at least " + std::to_string(required)
            + " ghost layers; Field.NG()=" + std::to_string(field.NG())
            + ". Increase mesh nGhost instead of silently downgrading order.");
    }
}

bool HJWeno::activeAxis(const Field& field, int axis) {
    return Math::isDirectionActiveIndex(axis) && axisLength(field, axis) > 1;
}

OneSidedDerivative HJWeno::derivative(const LevelSetField& levelSet,
                                      const std::vector<double>& values,
                                      int i, int j, int k,
                                      int axis,
                                      int order,
                                      const HJWenoWeightOptions& weights) {
    validateOrder(order, "HJWeno::derivative");
    validateWeightOptions(weights);
    if ((int)values.size() != levelSet.TotalSize()) {
        throw std::runtime_error(
            "HJWeno::derivative: phi snapshot size does not match LevelSetField.");
    }

    std::array<double, 7> minus{};
    std::array<double, 7> plus{};
    int count = 0;
    if (order == 3) {
        count = 3;
        const int mOffsets[3] = {-2, -1, 0};
        const int pOffsets[3] = {1, 0, -1};
        for (int n = 0; n < count; ++n) {
            minus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, mOffsets[n]);
            plus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, pOffsets[n]);
        }
    } else if (order == 5) {
        count = 5;
        const int mOffsets[5] = {-3, -2, -1, 0, 1};
        const int pOffsets[5] = {2, 1, 0, -1, -2};
        for (int n = 0; n < count; ++n) {
            minus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, mOffsets[n]);
            plus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, pOffsets[n]);
        }
    } else {
        count = 7;
        const int mOffsets[7] = {-4, -3, -2, -1, 0, 1, 2};
        const int pOffsets[7] = {3, 2, 1, 0, -1, -2, -3};
        for (int n = 0; n < count; ++n) {
            minus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, mOffsets[n]);
            plus[(size_t)n] = slopeAt(levelSet, values, i, j, k, axis, pOffsets[n]);
        }
    }

    OneSidedDerivative d;
    d.minus = reconstructLeft(minus.data(), order, weights);
    d.plus = reconstructLeft(plus.data(), order, weights);
    if (!std::isfinite(d.minus) || !std::isfinite(d.plus)) {
        throw std::runtime_error(
            "HJWeno::derivative: non-finite WENO derivative at "
            + cellText(i, j, k) + ".");
    }
    (void)count;
    return d;
}

double HJWeno::upwind(const OneSidedDerivative& d, double speed) {
    return speed >= 0.0 ? d.minus : d.plus;
}

double HJWeno::contravariantSpeed(const Field& field,
                                  int i, int j, int k,
                                  int axis,
                                  const Vector3& velocity) {
    if (axis == 0) {
        return velocity.x * field.XiX(i, j, k)
             + velocity.y * field.XiY(i, j, k)
             + velocity.z * field.XiZ(i, j, k);
    }
    if (axis == 1) {
        return velocity.x * field.EtX(i, j, k)
             + velocity.y * field.EtY(i, j, k)
             + velocity.z * field.EtZ(i, j, k);
    }
    return velocity.x * field.ZeX(i, j, k)
         + velocity.y * field.ZeY(i, j, k)
         + velocity.z * field.ZeZ(i, j, k);
}

Vector3 HJWeno::physicalGradient(const Field& field,
                                 int i, int j, int k,
                                 double dXi,
                                 double dEta,
                                 double dZeta) {
    return Vector3(
        field.XiX(i, j, k) * dXi
            + field.EtX(i, j, k) * dEta
            + field.ZeX(i, j, k) * dZeta,
        field.XiY(i, j, k) * dXi
            + field.EtY(i, j, k) * dEta
            + field.ZeY(i, j, k) * dZeta,
        field.XiZ(i, j, k) * dXi
            + field.EtZ(i, j, k) * dEta
            + field.ZeZ(i, j, k) * dZeta);
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
