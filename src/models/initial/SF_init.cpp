/// @file SF_init.cpp
/// @brief 初值设置与 legacy 配置迁移实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_init.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace SF{
    namespace Init {

        // 辅助工具：从 rho, u, v, w, p 计算总能 E
        // 公式：E = P/(gamma-1) + 0.5 * rho * |V|^2
        inline double calcEnergy(double rho, double u, double v, double w,
                                 double p, double gamma) {
            return p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v + w*w);
        }

        inline double calcEnergyFromTemperature(double rho,
                                                double u,
                                                double v,
                                                double w,
                                                double temperature,
                                                double gamma,
                                                double gasConstant) {
            if (!std::isfinite(rho) || rho <= 0.0
                || !std::isfinite(temperature) || temperature <= 0.0) {
                std::cerr << "[SF FATAL] invalid temperature initial state: rho="
                          << rho << ", T=" << temperature << std::endl;
                std::exit(1);
            }
            const double p = rho * gasConstant * temperature;
            return calcEnergy(rho, u, v, w, p, gamma);
        }

    template <typename T>
    void applyIC(SF::Field& field, const std::vector<BCSetting<T>>& icSettings, int vIdx) {
        for (const auto& ic : icSettings) {
            const auto& indices = field.getSet(ic.name);
            for (int idx : indices) {
                int i, j, k;
                field.getIJK(idx, i, j, k);
                field.setVal<T>(i, j, k, vIdx, ic.value);
            }
        }
    }

    inline double requirePositiveDensity(const SF::Field& field,
                                         int i,
                                         int j,
                                         int k,
                                         const char* context) {
        const double rho = field(i, j, k, RHO);
        if (!std::isfinite(rho) || rho <= 0.0) {
            std::ostringstream oss;
            oss << context << ": invalid rho before primitive velocity "
                << "conversion at (" << i << "," << j << "," << k
                << "), rho=" << rho;
            throw std::runtime_error(oss.str());
        }
        return rho;
    }

    void applyVelocityIC(SF::Field& field,
                         const std::vector<BCSetting<Vector3>>& icSettings) {
        for (const auto& ic : icSettings) {
            const auto& indices = field.getSet(ic.name);
            for (int idx : indices) {
                int i, j, k;
                field.getIJK(idx, i, j, k);
                const double rho =
                    requirePositiveDensity(field, i, j, k, "Init::applyVelocityIC");
                field(i, j, k, RU) = rho * ic.value.x;
                field(i, j, k, RV) = rho * ic.value.y;
                field(i, j, k, RW) = rho * ic.value.z;
            }
        }
    }

        void setupAll(SF::Field& field,
                      const FDM::InitialConditionConfig& config,
                      double idealGasGamma,
                      double idealGasConstant) {
            applyIC(field, config.density, SF::RHO);
            applyVelocityIC(field, config.velocity);

            // 将 p_IC (压力) 转换为总能 E 并写入场
            for (const auto& ic : config.pressure) {
                const auto& indices = field.getSet(ic.name);
                for (int idx : indices) {
                    int i, j, k;
                    field.getIJK(idx, i, j, k);
                    double rho = requirePositiveDensity(
                        field, i, j, k, "Init::setupAll p_IC");
                    double u   = field(i, j, k, RU) / rho;
                    double v   = field(i, j, k, RV) / rho;
                    double w   = field(i, j, k, RW) / rho;
                    double E_val = calcEnergy(
                        rho, u, v, w, ic.value, idealGasGamma);
                    field(i, j, k, E) = E_val;
                }
            }

            for (const auto& ic : config.temperature) {
                const auto& indices = field.getSet(ic.name);
                for (int idx : indices) {
                    int i, j, k;
                    field.getIJK(idx, i, j, k);
                    double rho = field(i, j, k, RHO);
                    if (!std::isfinite(rho) || rho <= 0.0) {
                        std::cerr << "[SF FATAL] invalid rho before T_IC at ("
                                  << i << "," << j << "," << k << "): "
                                  << rho << std::endl;
                        std::exit(1);
                    }
                    double u   = field(i, j, k, RU) / rho;
                    double v   = field(i, j, k, RV) / rho;
                    double w   = field(i, j, k, RW) / rho;
                    field(i, j, k, E) =
                        calcEnergyFromTemperature(
                            rho, u, v, w, ic.value,
                            idealGasGamma, idealGasConstant);
                }
            }

            // restart 的显式守恒能量最后写入，禁止再由 p/T 初值覆盖。
            applyIC(field, config.energy, SF::E);
        }

    } // namespace Init
} // namespace SF
