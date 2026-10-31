#include "solver/system/SF_lowering.h"
#include "core/system/SF_operationIds.h"
/// @file SF_providerResolver.cpp
/// @brief Freeze method-declared provider ownership and validate its capabilities.

#include "SF_providerResolver.h"

#include "SF_couplingStatus.h"
#include "SF_executionContract.h"
#include "SF_providerCatalog.h"
#include "SF_solvePlan.h"
#include "SF_scalarMethod.h"
#include "SF_eulerianAssembly.h"
#include "SF_pressureCoupling.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <set>
#include <stdexcept>

namespace SF::System {
namespace {

bool hasCapability(const std::vector<OperationCapability>& values,
                   OperationCapability wanted) {
    return std::find(values.begin(),values.end(),wanted) != values.end();
}

std::string validateNativeProvider(const std::string& selected,
        const ExecutableEquationSystem& equations,const CompiledNumericalSystem& numerics,
        const CompiledSolvePlan& plan,const ExecutionCapabilitySignature& signature,
        bool additionalContributions) {
    if (selected==ScalarOps::Provider || selected==ScalarOps::TransportProvider || selected==ScalarOps::SharedFluxProvider || selected==ScalarOps::CapacityProvider) {
        const bool temporal=!plan.compiledProgram.temporalParticipants.empty();
        if(!temporal) {
            // Static calls own only their declared leaves, not the surrounding
            // pressure program's clock, correction operations or loop topology.
            std::vector<std::string> owned;
            for(const auto& call:plan.compiledProgram.steps)if(call.backendProvider==selected)
                for(const auto& op:call.operations)owned.push_back(op.operation);
            std::set<std::string> found;
            bool clock=false;
            const auto local=[&](const auto& self,const SolvePlanNode& node)->bool {
                clock|=node.operation==OpIds::TimeCommit;
                if(node.provider==selected && !node.operation.empty()) {
                    if(std::find(owned.begin(),owned.end(),node.operation)==owned.end())return false;
                    found.insert(node.operation);
                }
                for(const auto& child:node.children)if(!self(self,child))return false;
                return true;
            };
            if(!local(local,plan.root) || std::any_of(owned.begin(),owned.end(),[&](const auto& op){return !found.count(op);}))
                return "Local scalar method has a missing or undeclared numerical operation";
            if(clock && !numerics.pressureMomentum)
                return "Static scalar calls have no implemented surrounding physical publication/clock contract";
            return {};
        }
        const auto reason=Lowering::validateTemporalPlan(plan,numerics.time.recipe.stageCount());
        if(!reason.empty())return reason;
        if (!plan.compiledProgram.hasTemporalRoot)
            return "scalar provider requires a compiled synchronous temporal group";
        int commits=0;
        std::vector<std::string> required{TemporalOps::Dt,TemporalOps::Open,TemporalOps::Ready,
            TemporalOps::RhsReady,TemporalOps::Close,TemporalOps::PublishReady,OpIds::TimeCommit};
        // A common plan can contain ordinary calls from other selected providers.
        // Validate their frozen operations too; provider availability is checked separately.
        for(const auto& call:plan.compiledProgram.steps)if(!call.temporalResidual)for(const auto& operation:call.operations)
            if(std::find(required.begin(),required.end(),operation.operation)==required.end())required.push_back(operation.operation);
        for (const auto& p:plan.compiledProgram.temporalParticipants) {
            if(!p.coefficientPrepare.empty())required.push_back(p.coefficientPrepare);
            if(!p.coefficientPublish.empty())required.push_back(p.coefficientPublish);
            for (const auto& op:p.preparation) required.push_back(op);
            for (const auto& op:{p.snapshot,p.prepareStage,p.rhs,p.advance,p.publish}) required.push_back(op);
        }
        std::vector<std::string> found;
        const auto check=[&](const auto& self,const SolvePlanNode& node)->bool {
            if (node.operation==OpIds::TimeCommit) ++commits;
            if (!node.operation.empty()) {
                if (std::find(required.begin(),required.end(),node.operation)==required.end()
                    || std::find(found.begin(),found.end(),node.operation)!=found.end()) return false;
                found.push_back(node.operation);
            }
            for (const auto& child:node.children) if (!self(self,child)) return false;
            return true;
        };
        if (!check(check,plan.root) || commits!=1 || found.size()!=required.size())
            return "scalar plan contains flow-only/wrong-target operations or an invalid physical commit";
        if (!std::isfinite(numerics.dt.maxDeltaT) || numerics.dt.maxDeltaT<=0)
            return "scalar fixed dt must be finite and positive";
    } else if (selected=="equation.thermal") {
        if(numerics.time.recipe.id()!=FDM::TimeRecipeId::ForwardEuler || plan.compiledProgram.temporalParticipants.empty()
            || signature.pressureConstraint)
            return "Unsupported thermal execution: Forward Euler synchronous scalar/flow group required; implicit/pressure thermal treatment is not implemented";
    } else if (selected=="flow.eulerian-pressure" || selected=="flow.eulerian-turbulence") {
        if (!hasConstraint(equations,"C_SHARED_PRESSURE") || !hasConstraint(equations,"C_VOLUME_FRACTION")
            || numerics.phaseTransport.convection!=FDM::PhaseConvectionScheme::Upwind
            || plan.compiledProgram.hasTemporalRoot)
            return "native Eulerian provider requires shared-pressure/volume closures, current Upwind backend and native non-temporal schedule";
        try {validateEulerianHow(plan.compiledProgram);}
        catch(const std::runtime_error& error){return error.what();}
    } else if (selected=="flow.homogeneous") {
        return "Unsupported native homogeneous execution: the current reconstructed convection kernel requires five-variable PerfectGas and cannot assemble generic-EOS component fluxes. No Rusanov fallback is permitted.";
    } else if (selected=="ibm.constraint") {
        const bool kkt=std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),
            [](const auto& call) { return call.backendProvider=="ibm.constraint"
                && std::any_of(call.operations.begin(),call.operations.end(),
                    [](const auto& op) {return op.operation==OpIds::IbmKktSolve;}); });
        if (kkt && (plan.compiledProgram.hasTemporalRoot || !signature.momentumPredictor))
            return "native immersed KKT requires an implemented pressure predictor/block schedule; explicit density stages are not a supported KKT predictor";
        if (!kkt) {
            if (!plan.compiledProgram.hasTemporalRoot)
                return "native immersed correction requires its implemented complete explicit predictor; pressure/implicit coupling is unavailable";
            const auto first=std::find_if(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),
                [](const auto& call) {return call.backendProvider=="ibm.constraint";});
            if (first==plan.compiledProgram.steps.end()
                || std::none_of(plan.compiledProgram.steps.begin(),first,[](const auto& call) {return call.temporalResidual;})
                || std::any_of(first,plan.compiledProgram.steps.end(),[](const auto& call) {return call.temporalResidual;}))
                return "native immersed correction must follow the complete physical predictor and precede commit; no pre-predictor or per-stage implementation exists";
        }
    } else if (selected=="flow.turbulence" || selected=="flow.turbulence-closure") {
        if (signature.pressureConstraint || !plan.compiledProgram.hasTemporalRoot)
            return "single-fluid transported turbulence requires explicit flow and the implemented unconstrained predictor treatment";
    } else if (selected=="flow.conservative" && std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),
            [](const auto& call) {return call.equationMethod=="MixtureBalance";})) {
        return "Unsupported native mixture execution: tracked phase-mass kernel exists, but compatible flow thermodynamic binding and native closure publication are not implemented. No fallback is permitted.";
    } else if (selected=="flow.conservative" && signature.pressureConstraint) {
        if (additionalContributions || !signature.conservativeState
            || !signature.momentumPredictor || !signature.pressureCorrection
            || signature.outerCorrectors!=1 || signature.pressureCorrectors<=0
            || signature.nonOrthogonalCorrectors!=0 || numerics.time.recipe.stageCount()!=1
            || equations.constraints.size()!=1)
            return "selected conservative pressure provider requires a single-stage single-fluid predictor, one outer pass and no non-orthogonal correction";
        for (const auto& layout:std::vector<std::pair<std::string,int>>{{"rho",0},{"rhoU",1},{"rhoE",4}}) {
            if (!equations.state.contains(layout.first)) return "conservative pressure requires STATE "+layout.first;
            const auto& symbol=equations.state.at(layout.first);
            if (symbol.storageBinding!=StorageBinding::PackedDistributed || symbol.storageKey!="conservative"
                || symbol.componentOffset!=layout.second || symbol.components!=(layout.first=="rhoU" ? 3 : 1))
                return "conservative pressure has incompatible physical storage: "+layout.first;
        }
        if (!equations.state.contains("U") || !equations.state.contains("p")
            || !producedBy(equations.state.at("U"),"closure.thermodynamicVelocity")
            || !producedBy(equations.state.at("p"),"closure.thermodynamic.p"))
            return "conservative pressure requires velocity and pressure views of the EOS-backed physical state";

    }
    return {};
}

} // namespace



std::string validateSolutionProviderContract(const StateRegistry& state,
        const EquationCompositionConfig& composition) {
    if (!composition.declared) return {};
    const auto& eos=composition.thermoDynamics.equationOfState;
    if (state.isSolution("U") && eos!="rhoConst")
        return "Unsupported provider contract: solution U requires variable-density primitive momentum/continuity and pressure-response operators; current PressureMomentum provider requires rhoConst. STATE and HOW were not rewritten.";
    if (state.isSolution("rho") && eos=="rhoConst")
        return "Unsupported provider contract: selected writable rho with rhoConst requires a closure-consistent density evolution and conservative energy provider. STATE and HOW were not rewritten.";
    if (composition.thermoDynamics.thermo=="janaf")
        return "Unsupported provider contract: native single-fluid energy requires a JANAF caloric provider; the current fused provider implements hConst.";
    if (composition.thermoDynamics.transport=="sutherland")
        return "Unsupported provider contract: native single-fluid coefficients require a Sutherland transport provider; the current provider implements const transport.";
    if (state.isSolution("U") && composition.thermoDynamics.transport.empty())
        return "Unsupported provider contract: PressureMomentum requires a registered transport closure.";
    if (state.isSolution("U")
        && std::find(composition.equations.begin(),composition.equations.end(),"Energy")!=composition.equations.end())
        return "Unsupported provider contract: PressureMomentum currently has no primitive Energy/Enthalpy equation executor.";
    return {};
}

ExecutionCapabilitySignature compileExecutionCapabilities(
        const ExecutableEquationSystem& system,
        const CompiledSolvePlan* plan) {
    ExecutionCapabilitySignature result;
    result.pressureConstraint = hasConstraint(system,"C_INCOMPRESSIBILITY");
    result.constantDensity = std::any_of(
        system.state.symbols().begin(),system.state.symbols().end(),[](const StateSymbol& item) {
            return item.id == "rho" && item.storageKey == "rhoConst";
        });
    result.conservativeState = std::any_of(
        system.state.symbols().begin(),system.state.symbols().end(),[](const StateSymbol& item) {
            return item.storageBinding == StorageBinding::PackedDistributed
                && item.storageKey == "conservative";
        });
    result.momentumPredictor = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& item) {
            return item.operatorBinding == OpIds::MomentumPredictor;
        });
    result.pressureCorrection = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& item) {
            return item.operatorBinding == OpIds::PressureCorrection;
        });
    if (plan && !plan->compiledProgram.steps.empty()) {
        for (const auto& call:plan->compiledProgram.steps) {
            result.momentumPredictor |= call.equationMethod=="PressureMomentum"
                || call.equationMethod=="ConservativePressureMomentum";
            result.pressureCorrection |= call.equationMethod=="PressureCorrection"
                || call.equationMethod=="ConservativePressureCorrection";
        }
    }
    if (plan && !plan->compiledProgram.steps.empty()) {
        result.outerCorrectors=1;
        result.pressureCorrectors=1;
        result.nonOrthogonalCorrectors=0;
        const auto collect=[&](const auto& self,const SolvePlanNode& node)->void {
            if (node.kind==PlanNodeKind::Loop) {
                if ((node.id=="outer" || node.id=="EE.outer")) result.outerCorrectors=node.repetitions;
                if ((node.id=="pressure" || node.id=="EE.pressure")) result.pressureCorrectors=node.repetitions;
                if ((node.id=="nonOrthogonal" || node.id=="EE.nonOrthogonal")) result.nonOrthogonalCorrectors=node.repetitions-1;
            }
            for (const auto& child:node.children) self(self,child);
        };
        collect(collect,plan->root);
    }
    return result;
}

std::vector<ResolvedOperationBinding> compileOperationBindings(
        const ExecutableEquationSystem& equations,
        const CompiledNumericalSystem& numerics,
        CompiledSolvePlan& plan,
        bool additionalContributions) {
    const auto signature=compileExecutionCapabilities(equations,&plan);
    const ProviderCatalog catalog=ProviderCatalog::builtIn();
    std::vector<ResolvedOperationBinding> bindings;
    // These operation identities were selected by actual numerical calls/HOW,
    // not by the legacy model catalog. Runtime must still bind every leaf.
    std::set<std::pair<OpId,std::string>> direct;
    const auto add=[&](const OpId& op,const std::string& owner){if(!op.empty())direct.emplace(op,owner);};
    for(const auto& call:plan.compiledProgram.steps)if(!call.legacyProvider) {
        add(call.backendOperation,call.backendProvider);
        for(const auto* actions:{&call.stepInputs,&call.iterationInputs,&call.iterationOutputs,&call.groupInputs})
            for(const auto& action:*actions)add(action.operation,action.provider);
        add(call.physicalPublication.operation,call.physicalPublication.provider);
        if(!call.physicalPublication.operation.empty())add(OpIds::TimeCommit,call.physicalPublication.provider);
        const auto fragment=[&](const auto& self,const SolvePlanNode& n)->void {
            add(n.operation,n.provider);for(const auto& child:n.children)self(self,child);
        };fragment(fragment,call.fragment);
        for(const auto& op:call.temporalPreparation)add(op,call.backendProvider);
        for(const auto& op:{call.temporalStepSize,call.temporalSnapshot,call.temporalStagePrepare,call.temporalRhs,
                call.temporalAdvance,call.temporalPublish})add(op,call.backendProvider);
    }
    if(plan.compiledProgram.hasTemporalRoot) {
        for(const auto* op:{TemporalOps::Dt,TemporalOps::Open,TemporalOps::Ready,TemporalOps::RhsReady,
                TemporalOps::Close,TemporalOps::PublishReady,OpIds::TimeCommit})add(op,TemporalOps::Provider);
    }
    const auto bind=[&](const auto& self,SolvePlanNode& leaf)->void {
        if (!leaf.operation.empty()) {
            ResolvedOperationBinding binding;
            binding.operation=leaf.operation;
            const auto declaration=std::find_if(equations.operations.begin(),equations.operations.end(),
                [&](const auto& item) { return item.operation==leaf.operation; });
            if (declaration==equations.operations.end()) {
                binding.reason=leaf.unsupportedReason.empty()
                    ? "no ExecutableOperation declares '"+leaf.operation+"'" : leaf.unsupportedReason;
            } else {
                if(direct.count({leaf.operation,leaf.provider})) {
                    binding.provider=leaf.provider;binding.status=BindingStatus::Resolved;
                } else binding=catalog.resolve(*declaration,leaf.provider);
                // STATE can reject a selected implementation; it cannot choose another.
                if (binding.status==BindingStatus::Resolved) {
                    std::vector<std::string> context;
                    std::vector<CapabilityBinding> boundContext;
                    for (const auto& boundary:equations.boundaryClosures) {
                        context.insert(context.end(),boundary.capabilities.begin(),boundary.capabilities.end());
                        boundContext.insert(boundContext.end(),boundary.boundCapabilities.begin(),boundary.boundCapabilities.end());
                    }
                    std::string requirement;
                    for (const auto& call:plan.compiledProgram.steps) if (call.backendProvider==binding.provider) {
                        requirement=missingCapability(plan.compiledProgram,context,call,boundContext);
                        if (!requirement.empty()) break;
                    }
                    const bool legacy=std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),
                        [&](const auto& call){return call.legacyProvider && call.backendProvider==binding.provider;});
                    const auto reason=!requirement.empty()?requirement:legacy ? validateNativeProvider(binding.provider,equations,numerics,plan,
                        signature,additionalContributions) : std::string{};
                    if (!reason.empty()) {
                        binding.status=BindingStatus::Unsupported;
                        binding.provider.clear();
                        binding.reason=reason;
                    }
                }
            }
            const auto existing=std::find_if(bindings.begin(),bindings.end(),
                [&](const auto& item) { return item.operation==binding.operation; });
            if (existing==bindings.end()) bindings.push_back(std::move(binding));
            else if (existing->provider!=binding.provider || existing->status!=binding.status) {
                existing->status=BindingStatus::Invalid;
                existing->provider.clear();
                existing->reason="conflicting compiled providers for '"+leaf.operation+"'";
            }
        }
        for (auto& child:leaf.children) self(self,child);
    };
    bind(bind,plan.root);
    return bindings;
}

RuntimeReport reportOperationBindings(
        const CompiledSolvePlan& plan,
        const std::vector<ResolvedOperationBinding>& bindings) {
    RuntimeReport report;
    report.requiredOperations = SolvePlanner::requiredOperations(plan);
    for (const OpId& id : report.requiredOperations) {
        const auto found = std::find_if(
            bindings.begin(),bindings.end(),
            [&](const ResolvedOperationBinding& binding) {
                return binding.operation == id;
            });
        if (found == bindings.end()) {
            throw std::runtime_error(
                "Provider resolver omitted plan operation '"+id+"'.");
        }
        if (found->status != BindingStatus::Resolved) {
            report.missingOperations.push_back(id);
            if (report.reason.empty()) report.reason = found->reason;
            if (found->status == BindingStatus::Invalid)
                report.status = RuntimeStatus::Invalid;
        }
    }
    if (report.status != RuntimeStatus::Invalid)
        report.status = report.missingOperations.empty()
            ? RuntimeStatus::Runnable : RuntimeStatus::Unsupported;
    if (!report.missingOperations.empty()) {
        report.reason += "; unresolved operations [";
        for (std::size_t i=0;i<report.missingOperations.size();++i) {
            if (i) report.reason += ", ";
            report.reason += report.missingOperations[i];
        }
        report.reason += "]";
    }
    return report;
}

} // namespace SF::System
