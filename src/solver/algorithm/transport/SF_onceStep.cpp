#include "SF_transportBinding.h"
#include "solver/boundary/SF_scalarBoundary.h"
#include "methods/numerics/structured/SF_structured.h"
#include <algorithm>
#include <cmath>
namespace SF::SolverAlgorithm {
namespace {
double eddyDiffusion(const Field& flow,
                                const std::vector<double>& viscosity,
                                const std::vector<double>& variable,
                                double invSigma,
                                int i, int j, int k) {
    double result = 0.0;
    double Jcell = 1.0 / std::max(std::abs(flow.Jac(i, j, k)), 1.0e-30);

    auto gammaAt = [&](int ii, int jj, int kk) -> double {
        return std::max(viscosity[flow.getIdx(ii, jj, kk)], 0.0) * invSigma;
    };

    auto phiAt = [&](int ii, int jj, int kk) -> double {
        return variable[flow.getIdx(ii, jj, kk)];
    };

    // ── XI direction ──
    if (Math::isDirectionActive(Math::XI)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i - 1, j, k, Math::XI, cfL);
        Math::faceMetrics(flow, i,     j, k, Math::XI, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i - 1, j, k) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i + 1, j, k));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i - 1, j, k));
        double fluxR = gamR * JgR * (phiAt(i + 1, j, k) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    // ── ETA direction ──
    if (Math::isDirectionActive(Math::ETA)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i, j - 1, k, Math::ETA, cfL);
        Math::faceMetrics(flow, i, j,     k, Math::ETA, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i, j - 1, k) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i, j + 1, k));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i, j - 1, k));
        double fluxR = gamR * JgR * (phiAt(i, j + 1, k) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    // ── ZETA direction ──
    if (Math::isDirectionActive(Math::ZETA)) {
        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i, j, k - 1, Math::ZETA, cfL);
        Math::faceMetrics(flow, i, j, k,     Math::ZETA, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                     / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                     / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i, j, k - 1) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i, j, k + 1));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i, j, k - 1));
        double fluxR = gamR * JgR * (phiAt(i, j, k + 1) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    }

    return result;
}

double molecularDiffusion(const Field& flow,
                          const std::vector<double>& viscosity,
                          const std::vector<double>& variable,
                          double sigma,
                          double laminarMu,
                          int i, int j, int k) {
    double result = 0.0;
    double Jcell = 1.0 / std::max(std::abs(flow.Jac(i, j, k)), 1.0e-30);

    auto gammaAt = [&](int ii, int jj, int kk) -> double {
        return std::max(laminarMu, 0.0)
             + sigma * std::max(viscosity[flow.getIdx(ii, jj, kk)], 0.0);
    };

    auto phiAt = [&](int ii, int jj, int kk) -> double {
        return variable[flow.getIdx(ii, jj, kk)];
    };

    auto addDirection = [&](Math::Dir d) {
        if (!Math::isDirectionActive(d)) return;

        int di, dj, dk;
        Math::dirOffset(d, di, dj, dk);

        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i - di, j - dj, k - dk, d, cfL);
        Math::faceMetrics(flow, i,      j,      k,      d, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                   / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                   / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i - di, j - dj, k - dk) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i + di, j + dj, k + dk));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i - di, j - dj, k - dk));
        double fluxR = gamR * JgR * (phiAt(i + di, j + dj, k + dk) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    };

    addDirection(Math::XI);
    addDirection(Math::ETA);
    addDirection(Math::ZETA);

    return result;
}

}
void applyTransportBoundary(const Field& field,FDM::ITransportModel& model) {
    auto data=model.onceStepTransport();
    if (data.variables[0]) {
        for(int n=0;n<2;++n) {
            if(!data.variables[n] || !data.boundary[n] || data.variables[n]->size()!=static_cast<size_t>(field.TotalSize()))
                throw std::runtime_error("Once-step transport boundary backing is incomplete.");
            Boundary::applyTransportScalarBoundary(field,*data.boundary[n],*data.variables[n]);
        }
    }
    model.applyBoundary(field);
}
void initializeTransport(const Field& field,FDM::ITransportModel& model) {
    applyTransportBoundary(field,model);model.refreshClosure(field);
}
void updateOnceStepTransport(const Field& field,FDM::ITransportModel& model,double dt) {
    auto data=model.onceStepTransport();
    if(!data.variables[0]) { model.refreshClosure(field);return; }
    if(!std::isfinite(dt) || dt<=0 || !data.variables[1] || !data.eddyViscosity || !data.cellClosure
        || data.variables[0]->size()!=static_cast<size_t>(field.TotalSize())
        || data.variables[1]->size()!=data.variables[0]->size() || data.eddyViscosity->size()!=data.variables[0]->size())
        throw std::runtime_error("Once-step transport requires valid dt, original arrays and a cell closure.");
    const int ng=field.NG();
    // The legacy sweep reads live neighbours updated earlier in this exact k/j/i order.
    // mu_t is frozen through the sweep; publishing a snapshot RHS would change the method.
    for(int k=ng;k<ng+field.NZ();++k)for(int j=ng;j<ng+field.NY();++j)for(int i=ng;i<ng+field.NX();++i) {
        if(field.CellFlag(i,j,k)!=FLUID_CELL)continue;
        auto cell=data.cellClosure(field,i,j,k);
        double diffusion[2];
        for(int n=0;n<2;++n) diffusion[n]=data.molecularDiffusion
            ? molecularDiffusion(field,*data.eddyViscosity,*data.variables[n],cell.diffusionWeight[n],data.molecularViscosity,i,j,k)
            : eddyDiffusion(field,*data.eddyViscosity,*data.variables[n],cell.diffusionWeight[n],i,j,k);
        cell.source[0]+=diffusion[0]/cell.density;
        if(data.molecularDiffusion) cell.source[1]+=(diffusion[1]+cell.crossDiffusion)/cell.density;
        else cell.source[1]+=diffusion[1]/cell.density;
        const int index=field.getIdx(i,j,k);
        (*data.variables[0])[index]=std::max(cell.values[0]+dt*cell.source[0],data.floors[0]);
        (*data.variables[1])[index]=std::max(cell.values[1]+dt*cell.source[1],data.floors[1]);
    }
    applyTransportBoundary(field,model);model.refreshClosure(field);
}
}
