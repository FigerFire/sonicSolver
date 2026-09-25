/// @file SF_heaviside.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_heaviside.h"

#include <cmath>
#include <string>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

constexpr double pi() {
    return 3.141592653589793238462643383279502884;
}

void requireEpsilon(double epsilon, const char* context) {
    if (!std::isfinite(epsilon) || epsilon <= 0.0) {
        throw std::runtime_error(
            std::string(context) + ": epsilon must be finite and > 0.");
    }
}

} // namespace

double Heaviside::sharp(double phi) {
    return phi >= 0.0 ? 1.0 : 0.0;
}

double Heaviside::regularized(double phi, double epsilon) {
    requireEpsilon(epsilon, "Heaviside::regularized");
    if (phi <= -epsilon) return 0.0;
    if (phi >= epsilon) return 1.0;
    return 0.5 * (1.0 + phi / epsilon
                  + std::sin(pi() * phi / epsilon) / pi());
}

double Heaviside::delta(double phi, double epsilon) {
    requireEpsilon(epsilon, "Heaviside::delta");
    if (std::abs(phi) >= epsilon) return 0.0;
    return 0.5 / epsilon * (1.0 + std::cos(pi() * phi / epsilon));
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
