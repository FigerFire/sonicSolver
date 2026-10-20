#pragma once
#include "core/interfaces/SF_termKernel.h"
#include "core/system/SF_formula.h"
#include "methods/numerics/structured/SF_structured.h"
#include <cmath>
namespace SF::Discretization {
/// Linear body-force density and its mechanical work use the same Stage
/// tracer. Neither callback mutates physical STATE or controls stage order.
inline System::StageConservativeSource linearScalarFeedback(const System::FormulaExpr& expression,bool work) {
        const auto& args=expression.arguments;using K=System::FormulaExpr::Kind;
        if (args.size()!=4 || args[1].kind!=K::Symbol || args[2].kind!=K::Constant || !std::isfinite(args[2].constant)
            || args[3].kind!=K::Constant || args[3].constant<0 || args[3].constant>2
            || args[3].constant!=std::floor(args[3].constant))
            throw std::runtime_error("Scalar feedback requires Stage symbol, finite coefficient and axis 0/1/2.");
        const double kappa=args[2].constant;const int axis=static_cast<int>(args[3].constant);
        System::StageConservativeSource source;source.reads={args[1].name};
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
