/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_heaviside.h
/// @brief Level Set正则化Heaviside和Delta函数。

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief 界面指示函数工具。
class Heaviside {
public:
    /// @brief 锐界面Heaviside, phi>=0为1, 否则为0。
    static double sharp(double phi);

    /// @brief Sussman/Osher常用余弦正则化Heaviside。
    /// @param phi Level Set值。
    /// @param epsilon 界面半厚度, 必须为正。
    static double regularized(double phi, double epsilon);

    /// @brief 正则化Dirac delta。
    /// @param phi Level Set值。
    /// @param epsilon 界面半厚度, 必须为正。
    static double delta(double phi, double epsilon);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
