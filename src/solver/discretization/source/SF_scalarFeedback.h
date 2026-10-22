#pragma once
#include "core/interfaces/SF_termKernel.h"

#include "methods/numerics/structured/SF_structured.h"
#include <cmath>
namespace SF::Discretization {
/// Linear body-force density and its mechanical work use the same Stage
/// tracer. Neither callback mutates physical STATE or controls stage order.
inline System::StageConservativeSource linearScalarFeedback(std::string target,double kappa,int axis,bool work) {
        System::StageConservativeSource source;source.reads={std::move(target)};
        if(work){source.reads.push_back("rho");source.reads.push_back("rhoU");}
        source.evaluate=[work,kappa,axis](Field& field,Residual& residual,const auto& reads) {
            Math::forFluidInterior(field,[&](int i,int j,int k) {
                const int cell=field.getIdx(i,j,k);double value=kappa*reads[0](cell,0);
                if(work) {
                    const double rho=reads[1](cell,0);
                    if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Scalar work has invalid Stage density.");
                    value*=reads[2](cell,axis)/rho;
                }
                if(!std::isfinite(value))throw std::runtime_error("Scalar feedback source is nonfinite.");
                residual.source(i,j,k,work?4:1+axis)+=value;
            });
        };return source;
}

}
