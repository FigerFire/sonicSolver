#include "solver/system/SF_providerCatalog.h"
#include "solver/system/SF_termProviderCatalog.h"
#include "models/physics/SF_sourceContribution.h"
#include "models/physics/gravity/SF_accelerationProvider.h"
#include "models/physics/mrf/SF_frameProvider.h"
#include "models/physics/heat/SF_wallFluxProvider.h"
#include "SF_field.h"
#include "core/residual/SF_residual.h"
#include "models/physics/mrf/SF_rotation.h"

#include <stdexcept>
#include <string>

using namespace SF::System;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    ExecutableOperation operation;
    operation.operation="pressure.solve";
    operation.requirements={OperationCapability::PressureLinearSolve};
    ProviderCatalog catalog;
    catalog.add({"serial",{OperationCapability::PressureLinearSolve}});
    auto result=catalog.resolve(operation,"serial");
    require(result.status==BindingStatus::Resolved && result.provider=="serial",
            "selected provider must resolve");
    result=catalog.resolve(operation,"distributed");
    require(result.status==BindingStatus::Unsupported,
            "missing selected provider must not fall back");
    catalog.add({"distributed",{OperationCapability::PressureLinearSolve}});
    require(catalog.resolve(operation,"serial").provider=="serial"
            && catalog.resolve(operation,"distributed").provider=="distributed",
            "registering another implementation changed WHICH selection");
    bool duplicateProviderRejected=false;
    try { catalog.add({"serial",{OperationCapability::PressureLinearSolve}}); }
    catch (const std::invalid_argument&) { duplicateProviderRejected=true; }
    require(duplicateProviderRejected,"duplicate implementation identity was accepted");
    operation.requirements={OperationCapability::ImmersedConstraint};
    require(catalog.resolve(operation,"serial").status==BindingStatus::Unsupported,
            "selected provider lacking a capability must not fall back");
    require(catalog.resolve(operation,{}).status==BindingStatus::Unsupported,
            "catalog inferred a provider without a WHICH binding");

    using Expr=FormulaExpr;
    const Equation momentum{"E_USER_FORCE",Expr::symbol("U"),
                           Expr::symbol("gravity"),{}};
    const auto gravity=Expr::op("source",{Expr::symbol("gravity")},
                                "momentum.gravity");
    const TermMatchContext termContext{
        momentum,gravity,"momentum.gravity","U","PressureMomentum",nullptr};
    auto sourceCatalog=TermProviderCatalog::builtIn();
    SF::ZoneVectorSetting acceleration;
    acceleration.value=SF::Vector3(0.01,0.0,0.0);
    for (auto descriptor:SF::Physics::Gravity::termProviders({acceleration}))
        sourceCatalog.addSource(std::move(descriptor));
    const auto builtIn=sourceCatalog.resolve(termContext);
    require(builtIn.status==BindingStatus::Resolved
            && builtIn.id=="source.gravity.primitive"
            && builtIn.primitiveSource,
            "custom Momentum role did not resolve primitive gravity");
    TermProviderCatalog lifetimeCatalog=TermProviderCatalog::builtIn();
    {
        SF::FDM::SourceConfig transient;
        transient.gravity.push_back(acceleration);
        SF::RotatingSetting rotation;
        rotation.axis=SF::Vector3(0.0,0.0,1.0);
        rotation.omega=2.0;
        transient.rotating.push_back(rotation);
        for (auto descriptor:SF::Physics::Gravity::termProviders(
                 transient.gravity)) lifetimeCatalog.addSource(std::move(descriptor));
        for (auto descriptor:SF::Physics::MRF::termProviders(
                 transient.rotating)) lifetimeCatalog.addSource(std::move(descriptor));
    }
    SF::Field geometry;
    geometry.setup(1,1,1,0);
    geometry.setBoundarySets({{"all",{0}}});
    const auto compiledGravity=lifetimeCatalog.resolve(termContext);
    require(compiledGravity.compiledDataAvailable
            && compiledGravity.owner=="models.physics.gravity"
            && compiledGravity.primitiveSource(
                geometry,0,0,0,SF::Vector3()).x==0.01,
            "gravity binding retained raw config or lost compiled acceleration");
    const auto mrfTerm=Expr::op("source",{Expr::symbol("MRF")},
                                "momentum.mrf");
    const auto compiledMRF=lifetimeCatalog.resolve(
        {momentum,mrfTerm,"momentum.mrf","U","PressureMomentum",nullptr});
    require(compiledMRF.compiledDataAvailable
            && compiledMRF.owner=="models.physics.mrf"
            && compiledMRF.primitiveSource(
                geometry,0,0,0,SF::Vector3(0.0,1.0,0.0)).x==-4.0,
            "MRF binding retained raw config or changed acceleration");
    SF::Field conservative;
    conservative.setup(1,1,1,0);
    conservative.setBoundarySets({{"all",{0}}});
    conservative(0,0,0,SF::RHO)=1.0;
    SF::Residual residual;
    residual.setup(1,1,1,5);
    const auto conservativeGravity=lifetimeCatalog.resolve(
        {momentum,gravity,"momentum.gravity","rhoU","ConservativeResidual",nullptr});
    require(conservativeGravity.compiledDataAvailable
            && conservativeGravity.conservativeSource,
            "gravity conservative binding lacks compiled data");
    conservativeGravity.conservativeSource(conservative,residual);
    require(residual.source(0,0,0,SF::RU)==0.01,
            "compiled gravity conservative kernel changed momentum source");
    residual.clearSource();
    conservative(0,0,0,SF::RV)=1.0;
    const auto conservativeMRF=lifetimeCatalog.resolve(
        {momentum,mrfTerm,"momentum.mrf","rhoU","ConservativeResidual",nullptr});
    require(conservativeMRF.compiledDataAvailable
            && conservativeMRF.conservativeSource,
            "MRF conservative binding lacks compiled data");
    conservativeMRF.conservativeSource(conservative,residual);
    require(residual.source(0,0,0,SF::RU)==-4.0,
            "compiled MRF conservative kernel changed Coriolis source");
    auto wallHeat=SF::Physics::WallHeat::termProvider({});
    const Equation energy{"energy",Expr::symbol("rhoE"),
                         Expr::symbol("wallHeat"),{}};
    const auto heatTerm=Expr::op("source",{Expr::symbol("wallHeat")},
                                 "energy.wallHeat");
    require(wallHeat.match({energy,heatTerm,"energy.wallHeat","rhoE",
                            "ConservativeResidual",nullptr})
            && !wallHeat.match({momentum,heatTerm,"energy.wallHeat","U",
                                 "PressureMomentum",nullptr}),
            "wallHeat provider matched pressure-only Momentum");
    require(sourceCatalog.resolve({momentum,gravity,"momentum.gravity",
                                   "p","PressureMomentum",nullptr}).status
                ==BindingStatus::Unsupported,
            "primitive source bound the wrong HOW Output");
    const auto unknown=Expr::op("source",{Expr::symbol("customForce")},
                                "momentum.customForce");
    require(sourceCatalog.resolve(
                {momentum,unknown,"momentum.customForce","U",
                 "PressureMomentum",nullptr}).status==BindingStatus::Unsupported,
            "unknown source must not fall back to a built-in source provider");
    const auto convection=Expr::op("div",{Expr::symbol("momentumFlux")},
                                   "momentum.convection");
    const auto wrongRecipe=SF::FDM::TermRecipe::diffusionRecipe(
        SF::FDM::TermRecipeId::Central4Explicit,
        SF::FDM::ViscousScheme::Central4);
    require(TermProviderCatalog::builtIn().resolve(
                {momentum,convection,"momentum.convection","U",
                 "PressureMomentum",&wrongRecipe}).status==BindingStatus::Unsupported,
            "unsupported convection recipe must not bind a diffusion provider");
    const auto conservativeRecipe=SF::FDM::builtInConvectionRecipe(
        SF::FDM::ConvectionScheme::WENO7,SF::FDM::FluxSplitter::Rusanov);
    const auto wrongFlux=Expr::op("div",{Expr::symbol("unimplementedFlux")});
    require(TermProviderCatalog::builtIn().resolve(
                {momentum,wrongFlux,"momentum.convection","rhoU",
                 "ConservativePressureMomentum",&conservativeRecipe}).status==BindingStatus::Unsupported,
            "conservative pressure silently reused its fused kernel for another flux operand");
    const auto accept=[](const TermMatchContext&) { return true; };
    TermProviderCatalog first;
    first.add({"zeta",accept,{}});
    first.add({"alpha",accept,{}});
    TermProviderCatalog reversed;
    reversed.add({"alpha",accept,{}});
    reversed.add({"zeta",accept,{}});
    const auto ambiguous=first.resolve(termContext);
    require(ambiguous.status==BindingStatus::Invalid
            && ambiguous.reason==reversed.resolve(termContext).reason
            && ambiguous.reason.find("alpha, zeta")!=std::string::npos,
            "term provider ambiguity must be sorted independently of registration");

    const std::vector<FormulaOperatorBinding> overrides{
        {"","source","zeta"},{momentum.id,"","zeta"},
        {momentum.id,"momentum.gravity","alpha"}};
    require(first.resolve(termContext,overrides).id=="alpha"
            && reversed.resolve(termContext,overrides).id=="alpha",
            "Occurrence override did not beat formula/global defaults.");
    require(first.resolve(termContext,{{"","source","alpha"},
                                      {momentum.id,"","zeta"}}).id=="zeta",
            "Equation default did not beat global operator default.");
    require(first.resolve(termContext,{{momentum.id,"momentum.gravity","missing"}})
                .status==BindingStatus::Unsupported,
            "Explicit unavailable provider fell back to another provider.");
    bool duplicateRejected=false;
    try {
        (void)first.resolve(termContext,{{momentum.id,"momentum.gravity","alpha"},
                                        {momentum.id,"momentum.gravity","zeta"}});
    } catch (const std::runtime_error&) { duplicateRejected=true; }
    require(duplicateRejected,"Duplicate occurrence bindings were accepted.");

    SystemContribution contribution;
    {
        SF::FDM::SourceConfig sources;
        sources.enabled.push_back(SF::FDM::SourceKind::Gravity);
        SF::Physics::SourceContribution::contribute(contribution,sources);
    }
    const SystemContribution copied=contribution;
    require(copied.mathematicalExtensions.size()==1
            && copied.mathematicalExtensions.front().families==std::vector<std::string>{"momentum"},
            "model contribution must extend Momentum mathematics");

    const SF::Vector3 origin;
    const auto coriolis=SF::Source::Rotating::rotatingFrameAcceleration(
        origin,SF::Vector3(0.0,1.0,0.0),origin,
        SF::Vector3(0.0,0.0,-2.0),origin);
    require(coriolis.x==-4.0 && coriolis.y==0.0 && coriolis.z==0.0,
            "MRF primitive acceleration changed Coriolis sign or magnitude");
}
