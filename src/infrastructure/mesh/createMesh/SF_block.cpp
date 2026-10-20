/// @file SF_block.cpp
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_block.h"
#include <cmath>

namespace SF
{
    Vector3 calculateTFI(double xi, double eta, double zeta,
                        const Vector3 v[8], Edge* edges[12]) {
        
        // --- Step 1: 基础三线性插值 (顶点贡献) ---
        Vector3 P_vertex(0, 0, 0);
        for (int kk = 0; kk < 2; ++kk) {
            double wz = kk ? zeta : (1.0 - zeta);
            for (int jj = 0; jj < 2; ++jj) {
                double wy = jj ? eta : (1.0 - eta);
                for (int ii = 0; ii < 2; ++ii) {
                    double wx = ii ? xi : (1.0 - xi);
                    
                    // OpenFOAM 映射表: map[kk][jj][ii]
                    static const int map[2][2][2] = {{{0, 1}, {3, 2}}, {{4, 5}, {7, 6}}};
                    int vi = map[kk][jj][ii];
                    P_vertex = P_vertex + (wx * wy * wz) * v[vi];
                }
            }
        }

        // --- Step 2: 3D TFI = P_ξ + P_η + P_ζ - P_ξη - P_ηζ - P_ξζ + T ---
        Vector3 P(0, 0, 0);

        // P_ξ: 单向插值沿 ξ 边
        for (int e = 0; e < 4; ++e) {
            double w = ((e&1) ? eta : (1-eta)) * ((e&2) ? zeta : (1-zeta));
            int v0 = (e==0)?0:((e==1)?3:((e==2)?4:7));
            int v1 = (e==0)?1:((e==1)?2:((e==2)?5:6));
            Vector3 ep = edges[e] ? edges[e]->evaluate(xi) : ((1-xi)*v[v0] + xi*v[v1]);
            P = P + w * ep;
        }
        // P_η: 单向插值沿 η 边
        for (int e = 0; e < 4; ++e) {
            double w = ((e&1) ? xi : (1-xi)) * ((e&2) ? zeta : (1-zeta));
            int v0 = (e==0)?0:((e==1)?1:((e==2)?4:5));
            int v1 = (e==0)?3:((e==1)?2:((e==2)?7:6));
            Vector3 ep = edges[4+e] ? edges[4+e]->evaluate(eta) : ((1-eta)*v[v0] + eta*v[v1]);
            P = P + w * ep;
        }
        // P_ζ: 单向插值沿 ζ 边
        for (int e = 0; e < 4; ++e) {
            double w = ((e&1) ? xi : (1-xi)) * ((e&2) ? eta : (1-eta));
            int v0 = (e==0)?0:((e==1)?1:((e==2)?3:2));
            int v1 = (e==0)?4:((e==1)?5:((e==2)?7:6));
            Vector3 ep = edges[8+e] ? edges[8+e]->evaluate(zeta) : ((1-zeta)*v[v0] + zeta*v[v1]);
            P = P + w * ep;
        }

        // - P_ξη: 双线性面 k=0/k=1, 在 ζ 混合
        {
            auto bilin_xy = [&](int kface) {
                int i0=kface*4, i1=i0+1, i2=i0+2, i3=i0+3;
                return (1-xi)*(1-eta)*v[i0] + xi*(1-eta)*v[i1] + (1-xi)*eta*v[i3] + xi*eta*v[i2];
            };
            P = P - ((1-zeta)*bilin_xy(0) + zeta*bilin_xy(1));
        }
        // - P_ηζ: 双线性面 i=0/i=1, 在 ξ 混合
        {
            auto bilin_yz = [&](int iface) {
                if (iface==0) return (1-eta)*(1-zeta)*v[0] + eta*(1-zeta)*v[3] + (1-eta)*zeta*v[4] + eta*zeta*v[7];
                else          return (1-eta)*(1-zeta)*v[1] + eta*(1-zeta)*v[2] + (1-eta)*zeta*v[5] + eta*zeta*v[6];
            };
            P = P - ((1-xi)*bilin_yz(0) + xi*bilin_yz(1));
        }
        // - P_ξζ: 双线性面 j=0/j=1, 在 η 混合
        {
            auto bilin_xz = [&](int jface) {
                if (jface==0) return (1-xi)*(1-zeta)*v[0] + xi*(1-zeta)*v[1] + (1-xi)*zeta*v[4] + xi*zeta*v[5];
                else          return (1-xi)*(1-zeta)*v[3] + xi*(1-zeta)*v[2] + (1-xi)*zeta*v[7] + xi*zeta*v[6];
            };
            P = P - ((1-eta)*bilin_xz(0) + eta*bilin_xz(1));
        }

        // + T (三线性顶点项 — 即 P_vertex)
        P = P + P_vertex;

        return P;
    }
}
