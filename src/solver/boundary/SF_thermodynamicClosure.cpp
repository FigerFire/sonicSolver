/// @file SF_thermodynamicClosure.cpp
/// @brief 边界 primitive 状态到守恒能量/热力学状态的闭合。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#include "SF_thermodynamicClosure.h"

#include "SF_boundary.h"
#include "SF_utility.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace SF {
namespace Boundary {

namespace {

const char* thermalTypeName(ThermalBCType type) {
    switch (type) {
    case ThermalBCType::FixedTemperature: return "fixedTemperature";
    case ThermalBCType::ZeroGradient: return "zeroGradient";
    case ThermalBCType::Adiabatic: return "adiabatic";
    case ThermalBCType::HeatFlux: return "heatFlux";
    case ThermalBCType::Empty: return "empty";
    }
    return "unknown";
}

double requireLegacyDensity(const Field& field,
                            int i,
                            int j,
                            int k,
                            const char* operation) {
    const double rho = field(i, j, k, RHO);
    if (!std::isfinite(rho) || rho <= 0.0) {
        throw std::runtime_error(
            std::string(operation) + ": invalid density at ("
            + std::to_string(i) + "," + std::to_string(j) + ","
            + std::to_string(k) + "): rho=" + std::to_string(rho)
            + " kg/m3. Boundary closure cannot replace this state with a "
              "density floor.");
    }
    return rho;
}

void nearestThermalInsideSample(const Field& field,
                                int i,
                                int j,
                                int k,
                                int axis,
                                int& si,
                                int& sj,
                                int& sk) {
    const int ng = field.NG();
    si = i;
    sj = j;
    sk = k;
    if (axis == 0) {
        if (i == ng && field.NX() > 1) si = i + 1;
        else if (i == field.NX() + ng - 1 && field.NX() > 1) si = i - 1;
    } else if (axis == 1) {
        if (j == ng && field.NY() > 1) sj = j + 1;
        else if (j == field.NY() + ng - 1 && field.NY() > 1) sj = j - 1;
    } else if (axis == 2) {
        if (k == ng && field.NZ() > 1) sk = k + 1;
        else if (k == field.NZ() + ng - 1 && field.NZ() > 1) sk = k - 1;
    }
}

double pointDistance(const Field& field,
                     int i0,
                     int j0,
                     int k0,
                     int i1,
                     int j1,
                     int k1) {
    const double dx = field.X(i1, j1, k1) - field.X(i0, j0, k0);
    const double dy = field.Y(i1, j1, k1) - field.Y(i0, j0, k0);
    const double dz = field.Z(i1, j1, k1) - field.Z(i0, j0, k0);
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!std::isfinite(d) || d <= 1.0e-30) {
        std::cerr << "[SF FATAL] thermal boundary found invalid point spacing "
                  << "between (" << i0 << "," << j0 << "," << k0
                  << ") and (" << i1 << "," << j1 << "," << k1
                  << "): " << d << std::endl;
        std::exit(1);
    }
    return d;
}

double thermalConductivity(double dynamicViscosity, double prandtl) {
    if (!std::isfinite(dynamicViscosity) || dynamicViscosity <= 0.0) {
        std::cerr << "[SF FATAL] heatFlux boundary requires positive finite "
                  << "dynamic viscosity, got " << dynamicViscosity
                  << std::endl;
        std::exit(1);
    }
    if (!std::isfinite(prandtl) || prandtl <= 0.0) {
        std::cerr << "[SF FATAL] heatFlux boundary requires positive finite "
                  << "Prandtl number, got " << prandtl << std::endl;
        std::exit(1);
    }
    const double cp = DefaultIdealGasGamma * DefaultIdealGasConstant
        / (DefaultIdealGasGamma - 1.0);
    return dynamicViscosity * cp / prandtl;
}

} // namespace

double pressureAt(const Field& field, int i, int j, int k) {
    if (field.hasStateModel()) {
        return field.thermodynamicState(i, j, k).pressure;
    }
    const double rho = requireLegacyDensity(
        field, i, j, k, "pressureAt");
    const double ke = 0.5 * (field(i, j, k, RU) * field(i, j, k, RU)
                           + field(i, j, k, RV) * field(i, j, k, RV)
                           + field(i, j, k, RW) * field(i, j, k, RW))
                    / rho;
    const double pressure =
        (DefaultIdealGasGamma - 1.0) * (field(i, j, k, E) - ke);
    if (!std::isfinite(pressure) || pressure <= 0.0) {
        throw std::runtime_error(
            "pressureAt: invalid pressure at (" + std::to_string(i) + ","
            + std::to_string(j) + "," + std::to_string(k) + "): p="
            + std::to_string(pressure)
            + " Pa. Correct the conservative boundary state.");
    }
    return pressure;
}

double temperatureAt(const Field& field, int i, int j, int k) {
    if (field.hasStateModel()) {
        const double temperature = field.thermodynamicState(i, j, k).temperature;
        if (!std::isfinite(temperature) || temperature <= 0.0) {
            throw std::runtime_error(
                "FluidStateModel thermal boundary received invalid temperature.");
        }
        return temperature;
    }
    const double p = Numerics::requirePhysicalState("temperatureAt",
                                                    field, i, j, k);
    const double rho = field(i, j, k, RHO);
    const double T = p / (rho * DefaultIdealGasConstant);
    if (!std::isfinite(T) || T <= 0.0) {
        std::cerr << "[SF FATAL] invalid temperature at ("
                  << i << "," << j << "," << k << "): T=" << T
                  << std::endl;
        std::exit(1);
    }
    return T;
}

void setEnergyFromPressure(Field& field, int i, int j, int k,
                           double pressure) {
    if (!std::isfinite(pressure) || pressure <= 0.0) {
        throw std::runtime_error(
            "setEnergyFromPressure: boundary pressure must be finite and "
            "positive at (" + std::to_string(i) + ","
            + std::to_string(j) + "," + std::to_string(k) + "): p="
            + std::to_string(pressure) + " Pa.");
    }
    if (field.hasStateModel()) {
        std::vector<double> q((size_t)field.NVar(), 0.0);
        for (int v = 0; v < field.NVar(); ++v) q[(size_t)v] = field(i,j,k,v);
        field(i,j,k,field.stateModel()->energyIndex()) =
            field.stateModel()->totalEnergyFromPressure(
                q.data(), field.NVar(), pressure);
        return;
    }
    const double rho = requireLegacyDensity(
        field, i, j, k, "setEnergyFromPressure");
    const double ke = 0.5 * (field(i, j, k, RU) * field(i, j, k, RU)
                           + field(i, j, k, RV) * field(i, j, k, RV)
                           + field(i, j, k, RW) * field(i, j, k, RW))
                    / rho;
    field(i, j, k, E) = pressure / (DefaultIdealGasGamma - 1.0) + ke;
}

void setEnergyFromTemperature(Field& field, int i, int j, int k,
                              double temperature) {
    if (!std::isfinite(temperature) || temperature <= 0.0) {
        std::cerr << "[SF FATAL] thermal boundary produced invalid temperature at ("
                  << i << "," << j << "," << k << "): T="
                  << temperature << std::endl;
        std::exit(1);
    }
    if (field.hasStateModel()) {
        std::vector<double> q((size_t)field.NVar(), 0.0);
        for (int v = 0; v < field.NVar(); ++v) q[(size_t)v] = field(i,j,k,v);
        field(i,j,k,field.stateModel()->energyIndex()) =
            field.stateModel()->totalEnergyFromTemperature(
                q.data(), field.NVar(), temperature);
        return;
    }
    const double rho = field(i, j, k, RHO);
    if (!std::isfinite(rho) || rho <= 0.0) {
        std::cerr << "[SF FATAL] thermal boundary cannot rebuild energy with "
                  << "invalid rho at (" << i << "," << j << "," << k
                  << "): rho=" << rho << std::endl;
        std::exit(1);
    }
    const double ke = 0.5 * (field(i, j, k, RU) * field(i, j, k, RU)
                           + field(i, j, k, RV) * field(i, j, k, RV)
                           + field(i, j, k, RW) * field(i, j, k, RW))
                    / rho;
    const double pressure = rho * DefaultIdealGasConstant * temperature;
    field(i, j, k, E) =
        pressure / (DefaultIdealGasGamma - 1.0) + ke;
}

void updateEnergyFromPressure(Field& field,
                              const std::vector<BCSetting<double>>& pSettings,
                              bool useILW,
                              int ilwOrder) {
    for (const auto& bc : pSettings) {
        forBoundarySettingPoints(field, bc, [&](int i, int j, int k, int axis) {
            if (useILW && bc.type != EMPTY) {
                using namespace ILW::BoundaryClosure;
                const auto constraint = bc.type == FIXED_VALUE ? ScalarConstraint::Value
                    : bc.type == SYMMETRY ? ScalarConstraint::Even : ScalarConstraint::NormalGradient;
                writeScalarTaylor(field,i,j,k,axis,ilwOrder,constraint,
                    bc.type == FIXED_VALUE ? bc.value : 0.,PressureGetter{},
                    [&](int a,int b,int c,double p) { setEnergyFromPressure(field,a,b,c,p); },
                    "pressure");
                return;
            }
            switch (bc.type) {
                case FIXED_VALUE: {


                    auto applyFixed = [&](int gi, int gj, int gk,
                                          int ri, int rj, int rk) {
                        int refI = ri, refJ = rj, refK = rk;
                        nearestBoundaryInteriorCell(field, ri, rj, rk,
                                                    refI, refJ, refK);
                        const double pGhost =
                            2.0 * bc.value - pressureAt(field, refI, refJ,
                                                        refK);
                        setEnergyFromPressure(field, gi, gj, gk, pGhost);
                    };

                    setEnergyFromPressure(field, i, j, k, bc.value);
                    forBoundaryGhostsAlongAxis(field, i, j, k, axis, applyFixed);
                    break;
                }
                case ZERO_GRADIENT:
                case SYMMETRY: {


                    auto applyZeroGradient = [&](int gi, int gj, int gk,
                                                 int ri, int rj, int rk) {
                        int refI = ri, refJ = rj, refK = rk;
                        nearestBoundaryInteriorCell(field, ri, rj, rk,
                                                    refI, refJ, refK);
                        setEnergyFromPressure(field, gi, gj, gk,
                                              pressureAt(field, refI, refJ,
                                                         refK));
                    };

                    const int ng = field.NG();
                    int srcI = i, srcJ = j, srcK = k;
                    if (axis == 0) {
                        if (i == ng && field.NX() > 1) srcI = i + 1;
                        else if (i == field.NX() + ng - 1 && field.NX() > 1)
                            srcI = i - 1;
                    } else if (axis == 1) {
                        if (j == ng && field.NY() > 1) srcJ = j + 1;
                        else if (j == field.NY() + ng - 1 && field.NY() > 1)
                            srcJ = j - 1;
                    } else if (axis == 2) {
                        if (k == ng && field.NZ() > 2) srcK = k + 1;
                        else if (k == field.NZ() + ng - 1 && field.NZ() > 2)
                            srcK = k - 1;
                    }
                    setEnergyFromPressure(field, i, j, k,
                                          pressureAt(field, srcI, srcJ,
                                                     srcK));
                    forBoundaryGhostsAlongAxis(field, i, j, k, axis, applyZeroGradient);
                    break;
                }
                case EMPTY: {
                    int di = 0, dj = 0, dk = 0;
                    if (axis == 0) di = 1;
                    else if (axis == 1) dj = 1;
                    else dk = 1;
                    if (axis == 0 && i == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i - b * di, j, k,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    if (axis == 0 && i == field.NX() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i + b * di, j, k,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    if (axis == 1 && j == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i, j - b * dj, k,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    if (axis == 1 && j == field.NY() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i, j + b * dj, k,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    if (axis == 2 && k == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i, j, k - b * dk,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    if (axis == 2 && k == field.NZ() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b) {
                            setEnergyFromPressure(field, i, j, k + b * dk,
                                                  pressureAt(field, i, j, k));
                        }
                    }
                    break;
                }
                default:
                    break;
            }
        });
    }
}

void updateEnergyFromThermalBoundary(
    Field& field,
    const std::vector<ThermalBCSetting>& thermalSettings,
    double dynamicViscosity,
    double prandtl,
    bool useILW,
    int ilwOrder) {
    if (thermalSettings.empty()) return;


    for (const auto& bc : thermalSettings) {
        const auto& allSets = field.getAllSets();
        auto setIt = allSets.find(bc.name);
        if (setIt == allSets.end()) {
            std::cerr << "[SF FATAL] thermal boundary set '" << bc.name
                      << "' is not a mesh set." << std::endl;
            std::exit(1);
        }
        // The global set identity survives decomposition, but a patch can
        // have no local members. There is then no thermal face to reconstruct.
        if (setIt->second.empty()) continue;

        const int axis = (bc.type == ThermalBCType::Empty)
            ? Geometry::boundaryAxisForSet(field, bc.name)
            : Geometry::activeBoundaryAxisForSet(field, bc.name);

        if ((bc.type == ThermalBCType::FixedTemperature
             || bc.type == ThermalBCType::HeatFlux)
            && !std::isfinite(bc.value)) {
            std::cerr << "[SF FATAL] thermal boundary '" << bc.name
                      << "' has non-finite " << thermalTypeName(bc.type)
                      << " value: " << bc.value << std::endl;
            std::exit(1);
        }

        for (int idx : setIt->second) {
            int i = 0, j = 0, k = 0;
            field.getIJK(idx, i, j, k);
            if (!Geometry::isPhysicalPoint(field, i, j, k)) continue;

            if (useILW && bc.type != ThermalBCType::Empty) {
                using namespace ILW::BoundaryClosure;
                ScalarConstraint constraint=ScalarConstraint::NormalGradient;
                double prescribed=0.;
                if (bc.type == ThermalBCType::FixedTemperature) {
                    constraint=ScalarConstraint::Value;
                    prescribed=bc.value;
                    if (prescribed<=0.)
                        throw std::runtime_error("ILW fixedTemperature requires positive temperature.");
                } else if (bc.type == ThermalBCType::HeatFlux) {
                    int si=i,sj=j,sk=k;
                    nearestThermalInsideSample(field,i,j,k,axis,si,sj,sk);
                    const double kappa=field.hasStateModel()
                        ? field.thermodynamicState(si,sj,sk).thermalConductivity
                        : thermalConductivity(dynamicViscosity,prandtl);
                    if (!std::isfinite(kappa) || kappa<=0.)
                        throw std::runtime_error("ILW heatFlux requires positive finite thermal conductivity.");
                    // Preserve the existing heatFlux convention: k*dT/dn=value.
                    prescribed=bc.value/kappa;
                }
                auto temperature=[](const Field& f,int a,int b,int c) { return temperatureAt(f,a,b,c); };
                writeScalarTaylor(field,i,j,k,axis,ilwOrder,constraint,prescribed,temperature,
                    [&](int a,int b,int c,double t) { setEnergyFromTemperature(field,a,b,c,t); },"thermal");
                continue;
            }
            switch (bc.type) {
                case ThermalBCType::FixedTemperature: {
                    if (bc.value <= 0.0) {
                        std::cerr << "[SF FATAL] fixedTemperature boundary '"
                                  << bc.name << "' must be > 0, got "
                                  << bc.value << std::endl;
                        std::exit(1);
                    }
                    setEnergyFromTemperature(field, i, j, k, bc.value);
                    auto applyFixed = [&](int gi, int gj, int gk,
                                          int ri, int rj, int rk) {
                        const double realT = temperatureAt(field, ri, rj, rk);
                        const double ghostT = 2.0 * bc.value - realT;
                        setEnergyFromTemperature(field, gi, gj, gk, ghostT);
                    };
                    forBoundaryGhostsAlongAxis(
                        field, i, j, k, axis, applyFixed);
                    break;
                }
                case ThermalBCType::ZeroGradient:
                case ThermalBCType::Adiabatic: {
                    int si = i, sj = j, sk = k;
                    nearestThermalInsideSample(field, i, j, k, axis,
                                               si, sj, sk);
                    setEnergyFromTemperature(field, i, j, k,
                                             temperatureAt(field, si, sj, sk));
                    auto copyTemperature = [&](int gi, int gj, int gk,
                                               int ri, int rj, int rk) {
                        setEnergyFromTemperature(
                            field, gi, gj, gk,
                            temperatureAt(field, ri, rj, rk));
                    };
                    forBoundaryGhostsAlongAxis(
                        field, i, j, k, axis, copyTemperature);
                    break;
                }
                case ThermalBCType::HeatFlux: {
                    int si = i, sj = j, sk = k;
                    nearestThermalInsideSample(field, i, j, k, axis,
                                               si, sj, sk);
                    const double inwardDistance =
                        pointDistance(field, si, sj, sk, i, j, k);
                    const double kappa = field.hasStateModel()
                        ? field.thermodynamicState(si, sj, sk).thermalConductivity
                        : thermalConductivity(dynamicViscosity, prandtl);
                    if (!std::isfinite(kappa) || kappa <= 0.0) {
                        throw std::runtime_error(
                            "FluidStateModel heatFlux boundary requires explicit positive thermal conductivity.");
                    }
                    const double boundaryT =
                        temperatureAt(field, si, sj, sk)
                        + bc.value * inwardDistance / kappa;
                    setEnergyFromTemperature(field, i, j, k, boundaryT);

                    auto applyHeatFlux = [&](int gi, int gj, int gk,
                                             int, int, int) {
                        const double outwardDistance =
                            pointDistance(field, i, j, k, gi, gj, gk);
                        const double ghostT =
                            boundaryT + bc.value * outwardDistance / kappa;
                        setEnergyFromTemperature(field, gi, gj, gk, ghostT);
                    };
                    forBoundaryGhostsAlongAxis(
                        field, i, j, k, axis, applyHeatFlux);
                    break;
                }
                case ThermalBCType::Empty: {
                    auto copyTemperature = [&](int gi, int gj, int gk) {
                        setEnergyFromTemperature(
                            field, gi, gj, gk,
                            temperatureAt(field, i, j, k));
                    };
                    if (axis == 0 && i == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i - b, j, k);
                    }
                    if (axis == 0 && i == field.NX() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i + b, j, k);
                    }
                    if (axis == 1 && j == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i, j - b, k);
                    }
                    if (axis == 1 && j == field.NY() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i, j + b, k);
                    }
                    if (axis == 2 && k == field.NG()) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i, j, k - b);
                    }
                    if (axis == 2 && k == field.NZ() + field.NG() - 1) {
                        for (int b = 1; b <= field.NG(); ++b)
                            copyTemperature(i, j, k + b);
                    }
                    break;
                }
            }
        }
        // Intersecting patches use the explicitly ordered thermal laws. The
        // next patch must see this patch's energy publication, not cached T.
        field.invalidateThermodynamicCache();
    }
}

} // namespace Boundary
} // namespace SF
