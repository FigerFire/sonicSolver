/// @file SF_riemann.cpp
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include <cmath>
#include "SF_riemann.h"
#include "SF_utility.h"

namespace SF {
namespace Riemann {

void splitStegerWarming(const double* Q, const double* n,
                        double gamma, double* f_pos, double* f_neg) {
    double r    = Q[0];
    double p  = Numerics::requirePhysicalState(
        "splitStegerWarming", Q[0], Q[1], Q[2], Q[3], Q[4], gamma);
    double invR = 1.0 / r;
    double u    = Q[1] * invR;
    double v    = Q[2] * invR;
    double w    = Q[3] * invR;
    double e    = Q[4];

    double q2 = u*u + v*v + w*w;
    double a  = std::sqrt(gamma * p * invR);
    double un = u * n[0] + v * n[1] + w * n[2];

    double lam[5] = {un, un, un, un + a, un - a};
    double lp[5], lm[5];
    const double eps_sw = 0.01;
    for (int i = 0; i < 5; ++i) {
        if (std::abs(lam[i]) < 2.0 * eps_sw) {
            lp[i] = 0.5 * (lam[i]*lam[i] / (4.0 * eps_sw) + lam[i] + eps_sw);
            lm[i] = lam[i] - lp[i];
        } else {
            lp[i] = 0.5 * (lam[i] + std::abs(lam[i]));
            lm[i] = 0.5 * (lam[i] - std::abs(lam[i]));
        }
    }

    double common_p = r / (2.0 * gamma);
    double W1p = 2.0 * (gamma - 1.0) * lp[0] + lp[3] + lp[4];
    double W2p = lp[3] - lp[4];
    double W3p = lp[3] + lp[4];

    f_pos[0] = common_p * W1p;
    f_pos[1] = common_p * (u * W1p + a * n[0] * W2p);
    f_pos[2] = common_p * (v * W1p + a * n[1] * W2p);
    f_pos[3] = common_p * (w * W1p + a * n[2] * W2p);
    f_pos[4] = common_p * (
        0.5 * q2 * W1p + a * un * W2p
        + (a*a / (gamma - 1.0)) * W3p);

    double W1m = 2.0 * (gamma - 1.0) * lm[0] + lm[3] + lm[4];
    double W2m = lm[3] - lm[4];
    double W3m = lm[3] + lm[4];

    f_neg[0] = common_p * W1m;
    f_neg[1] = common_p * (u * W1m + a * n[0] * W2m);
    f_neg[2] = common_p * (v * W1m + a * n[1] * W2m);
    f_neg[3] = common_p * (w * W1m + a * n[2] * W2m);
    f_neg[4] = common_p * (
        0.5 * q2 * W1m + a * un * W2m
        + (a*a / (gamma - 1.0)) * W3m);
}

void buildCharacteristicMatrix(const double* Q_avg, const double* n,
                               double gamma, double* L, double* R) {
    double r    = Q_avg[0];
    double p  = Numerics::requirePhysicalState("buildCharacteristicMatrix",
                                               Q_avg[0], Q_avg[1], Q_avg[2],
                                               Q_avg[3], Q_avg[4], gamma);
    double invR = 1.0 / r;
    double u    = Q_avg[1] * invR;
    double v    = Q_avg[2] * invR;
    double w    = Q_avg[3] * invR;
    double e    = Q_avg[4];

    double q2 = u*u + v*v + w*w;
    double c  = std::sqrt(gamma * p * invR);
    double H  = (e + p) * invR;

    double n1 = n[0], n2 = n[1], n3 = n[2];
    double l1, l2, l3, m1, m2, m3;

    if (std::abs(n3) <= std::abs(n2)) {
        double ss = std::sqrt(n1*n1 + n2*n2);
        l1 = -n2/ss; l2 = n1/ss; l3 = 0.0;
    } else {
        double ss = std::sqrt(n1*n1 + n3*n3);
        l1 = -n3/ss; l2 = 0.0; l3 = n1/ss;
    }
    m1 = n2*l3 - n3*l2;
    m2 = n3*l1 - n1*l3;
    m3 = n1*l2 - n2*l1;

    double un = u*n1 + v*n2 + w*n3;
    double ul = u*l1 + v*l2 + w*l3;
    double um = u*m1 + v*m2 + w*m3;
    double b1 = (gamma - 1.0) / (c*c);
    double b2 = 0.5 * b1;
    double x1 = 1.0 / (2.0 * c);

    auto setL = [&](int r, int c, double val) { L[r*5 + c] = val; };
    setL(0,0, 1.0 - b2*q2);    setL(0,1, b1*u);        setL(0,2, b1*v);        setL(0,3, b1*w);        setL(0,4, -b1);
    setL(1,0, -ul);             setL(1,1, l1);          setL(1,2, l2);          setL(1,3, l3);          setL(1,4, 0.0);
    setL(2,0, -um);             setL(2,1, m1);          setL(2,2, m2);          setL(2,3, m3);          setL(2,4, 0.0);
    // The acoustic left eigenvectors contain (gamma-1)|u|^2/(4c^2).
    // b2 already equals (gamma-1)/(2c^2), so the kinetic contribution
    // is 0.5*b2*q2.  Using b2*q2 makes L and R non-inverse and changes a
    // constant conservative state during characteristic reconstruction.
    setL(3,0, 0.5*b2*q2 + x1*un); setL(3,1,-x1*n1 - b2*u); setL(3,2,-x1*n2 - b2*v); setL(3,3,-x1*n3 - b2*w); setL(3,4, b2);
    setL(4,0, 0.5*b2*q2 - x1*un); setL(4,1, x1*n1 - b2*u); setL(4,2, x1*n2 - b2*v); setL(4,3, x1*n3 - b2*w); setL(4,4, b2);

    auto setR = [&](int r, int c, double val) { R[r*5 + c] = val; };
    setR(0,0, 1.0);     setR(0,1, 0.0); setR(0,2, 0.0); setR(0,3, 1.0);           setR(0,4, 1.0);
    setR(1,0, u);       setR(1,1, l1);   setR(1,2, m1);   setR(1,3, u - c*n1);     setR(1,4, u + c*n1);
    setR(2,0, v);       setR(2,1, l2);   setR(2,2, m2);   setR(2,3, v - c*n2);     setR(2,4, v + c*n2);
    setR(3,0, w);       setR(3,1, l3);   setR(3,2, m3);   setR(3,3, w - c*n3);     setR(3,4, w + c*n3);
    setR(4,0, 0.5*q2);  setR(4,1, ul);   setR(4,2, um);   setR(4,3, H - c*un);     setR(4,4, H + c*un);
}

void splitLaxFriedrichs(const double* Q, const double* n,
                        double gamma, double* f_pos, double* f_neg) {
    double r    = Q[0];
    double p  = Numerics::requirePhysicalState(
        "splitLaxFriedrichs", Q[0], Q[1], Q[2], Q[3], Q[4], gamma);
    double invR = 1.0 / r;
    double u    = Q[1] * invR;
    double v    = Q[2] * invR;
    double w    = Q[3] * invR;
    double e    = Q[4];

    double q2 = u*u + v*v + w*w;
    double a  = std::sqrt(gamma * p * invR);

    double un = u * n[0] + v * n[1] + w * n[2];
    double lam_max = std::abs(un) + a;                          // 最大特征值

    // 物理通量 F_phys(Q) · n̂
    double F[5];
    F[0] = r * un;                            // ρ * u_n
    F[1] = Q[1] * un + p * n[0];             // ρu*u_n + p*n_x
    F[2] = Q[2] * un + p * n[1];             // ρv*u_n + p*n_y
    F[3] = Q[3] * un + p * n[2];             // ρw*u_n + p*n_z
    F[4] = (e + p) * un;                      // (E+p)*u_n

    for (int i = 0; i < 5; ++i) {
        f_pos[i] = 0.5 * (F[i] + lam_max * Q[i]);
        f_neg[i] = 0.5 * (F[i] - lam_max * Q[i]);
    }
}

} // namespace Riemann
} // namespace SF
