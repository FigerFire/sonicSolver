#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::System {
namespace ScalarOps {
inline constexpr const char* Dt="scalar.dt.fixed";
inline constexpr const char* Begin="scalar.snapshot";
inline constexpr const char* Stage="scalar.stage";
inline constexpr const char* Commit="scalar.publish";
inline constexpr const char* Provider="equation.scalar-central2";
inline constexpr const char* TransportProvider="equation.scalar-transport";
inline constexpr const char* StageTransportProvider="equation.stage-scalar-transport";
inline constexpr const char* SharedFluxProvider="equation.scalar-shared-flux";
}
/// Immutable AST lowering output. No physical storage or temporal coefficients.
struct CompiledScalarEquation {
    std::string target;
    double diffusivity=0;
    bool conservativeTransport=false;
    bool sharedMassFlux=false;
    bool primitiveStorage=false,modelCoefficients=false;
    std::vector<std::string> modelTargets;
    std::string diffusivityState;
    std::vector<std::string> sources;
    std::vector<std::string> stageReads;
    FormulaValueKernel rhs;
};
void addScalarMethod(ProviderRegistry& registry);
}
