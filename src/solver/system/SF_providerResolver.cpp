#include "core/system/SF_operationIds.h"
/// @file SF_providerResolver.cpp
/// @brief Freeze method-declared provider ownership and validate its capabilities.

#include "SF_providerResolver.h"

#include "SF_couplingStatus.h"
#include "SF_legacyNumerics.h"
#include "SF_providerCatalog.h"
#include "SF_solvePlan.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::System {
namespace {

bool hasCapability(const std::vector<OperationCapability>& values,
                   OperationCapability wanted) {
    return std::find(values.begin(),values.end(),wanted) != values.end();
}

bool pressureOperatorRequirementsSatisfied(
        const ExecutableEquationSystem& equations,
        const ExecutionCapabilitySignature& signature,
        const std::vector<LegacyExecutionPolicy>& policies,
        const CompiledNumericalSystem& numerics,
        bool additionalContributions) {
    if (additionalContributions || !policies.empty()
        || numerics.time.recipe.topology()!=FDM::TimeTopology::ExplicitStages
        || numerics.time.recipe.stageCount()!=1
        || !signature.pressureConstraint || !signature.constantDensity
        || signature.conservativeState
        || !signature.momentumPredictor || !signature.pressureCorrection
        || signature.auxiliarySchedule || signature.outerCorrectors<=0
        || signature.pressureCorrectors<=0
        || signature.nonOrthogonalCorrectors!=0
        || !numerics.pressureOperator
        || numerics.pressureFaceCoupling!=PressureFaceCoupling::RhieChow)
        return false;
    const auto density=std::find_if(equations.state.symbols().begin(),
        equations.state.symbols().end(),[](const StateSymbol& unknown) {
            return unknown.id=="rho" && unknown.constantValue.has_value();
        });
    if (density==equations.state.symbols().end()
        || !std::isfinite(*density->constantValue)
        || *density->constantValue<=0.0) return false;
    if (!equations.registry.contains("momentum")
        || std::count_if(equations.constraints.begin(),equations.constraints.end(),
            [](const ConstraintDescriptor& constraint) { return constraint.id=="C_INCOMPRESSIBILITY"; })!=1
        || equations.constraints.size()!=1) return false;
    const auto sourceBinding=std::find_if(numerics.operators.begin(),
        numerics.operators.end(),[](const CompiledSpatialBinding& binding) {
            return binding.equationMethod=="PressureMomentum"
                && binding.output=="U";
        });
    if (sourceBinding==numerics.operators.end()) return false;
    const Equation* formula=nullptr;
    try { formula=&equations.registry.at(sourceBinding->formulaId); }
    catch (const std::runtime_error&) { return false; }
    int transient=0,divergence=0,diffusion=0;
    const auto checkLeft=[&](const auto& self,const FormulaExpr& node)->bool {
        if (node.kind==FormulaExpr::Kind::Add)
            return self(self,node.arguments[0]) && self(self,node.arguments[1]);
        if (node.kind!=FormulaExpr::Kind::Operator) return false;
        if (node.name=="ddt" && node.arguments.size()==1
            && node.arguments[0].kind==FormulaExpr::Kind::Symbol
            && node.arguments[0].name=="U") { ++transient; return true; }
        if (node.name=="div" && node.arguments.size()==1
            && node.arguments[0].kind==FormulaExpr::Kind::Symbol
            && node.arguments[0].name=="momentumFlux") {
            ++divergence; return true;
        }
        if (node.name=="diffusion" && node.arguments.size()==2
            && node.arguments[0].kind==FormulaExpr::Kind::Symbol
            && node.arguments[0].name=="nu"
            && node.arguments[1].kind==FormulaExpr::Kind::Symbol
            && node.arguments[1].name=="U") { ++diffusion; return true; }
        return false;
    };
    if (!checkLeft(checkLeft,formula->lhs)
        || transient!=1 || divergence!=1 || diffusion>1) return false;
    const auto checkRight=[&](const auto& self,const FormulaExpr& node)->bool {
        if (node.kind==FormulaExpr::Kind::Add)
            return self(self,node.arguments[0]) && self(self,node.arguments[1]);
        if (node.kind==FormulaExpr::Kind::Constant) return node.constant==0.0;
        const auto symbol=node.kind==FormulaExpr::Kind::Symbol ? node.name
            : node.kind==FormulaExpr::Kind::Operator
                && node.name=="source" && node.arguments.size()==1
                && node.arguments[0].kind==FormulaExpr::Kind::Symbol
                ? node.arguments[0].name : std::string{};
        if (symbol.empty()) return false;
        return std::any_of(numerics.operators.begin(),numerics.operators.end(),
            [&](const CompiledSpatialBinding& binding) {
                return binding.formulaId==formula->id && binding.output=="U"
                    && binding.operatorName=="source"
                    && binding.side==CompiledSpatialBinding::Side::Right
                    && binding.primary==symbol && binding.primitiveSource;
            });
    };
    if (!checkRight(checkRight,formula->rhs)) return false;
    bool convection=false;
    for (const CompiledSpatialBinding& term:numerics.operators) {
        if (term.formulaId!=formula->id || term.output!="U") continue;
        if (term.operatorName=="div") {
            if (term.primary!="momentumFlux" || !term.recipe
                || term.recipe->id()!=FDM::TermRecipeId::PrimitiveUpwind1)
                return false;
            convection=true;
        } else if (term.operatorName=="diffusion") {
            if (term.primary!="U" || term.secondary!="nu" || !term.recipe
                || term.recipe->id()!=FDM::TermRecipeId::Central2Explicit)
                return false;
        } else if (term.operatorName=="source") {
            if (term.side!=CompiledSpatialBinding::Side::Right
                || !term.primitiveSource) return false;
        }
    }
    return convection;
}

std::string validateNativeProvider(const std::string& selected,
        const ExecutableEquationSystem& equations,const CompiledNumericalSystem& numerics,
        const CompiledSolvePlan& plan,const ExecutionCapabilitySignature& signature,
        const std::vector<LegacyExecutionPolicy>& policies,bool additionalContributions) {
    if (selected=="flow.eulerian-pressure" || selected=="flow.eulerian-turbulence") {
        if (!hasConstraint(equations,"C_SHARED_PRESSURE") || !hasConstraint(equations,"C_VOLUME_FRACTION")
            || numerics.phaseTransport.convection!=FDM::PhaseConvectionScheme::Upwind
            || plan.compiledProgram.hasTemporalRoot || !policies.empty())
            return "native Eulerian provider requires shared-pressure/volume closures, current Upwind backend and no extra legacy constraint schedule";
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
        if (!equations.boundaryClosures.empty() || std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),
            [](const auto& call) {return call.backendProvider=="ibm.constraint";}))
            return "native turbulence + IBM requires unimplemented stage/velocity/wall-distance/boundary coupling contracts";
        if (signature.pressureConstraint || !plan.compiledProgram.hasTemporalRoot
            || !equations.constraints.empty() || !policies.empty())
            return "single-fluid transported turbulence requires explicit flow, no pressure/IBM constraints";
    } else if (selected=="flow.pressure-operators" || selected=="flow.rhie-chow") {
        if (!pressureOperatorRequirementsSatisfied(equations,signature,policies,numerics,additionalContributions))
            return "selected native pressure provider requires a supported constant-density pressure-multiplier state, equation, recipe and fixed-time capability";
        for (const auto& layout:std::vector<std::pair<std::string,int>>{{"U",3},{"p",1},{"phi",1}}) {
            if (!equations.state.contains(layout.first)) return "selected pressure provider requires STATE "+layout.first;
            const auto& symbol=equations.state.at(layout.first);
            if (symbol.components!=layout.second || symbol.derivation!=StateDerivation::None || symbol.constantValue)
                return "selected pressure provider has incompatible STATE contract: "+layout.first;
        }
    } else if (selected=="flow.conservative" && signature.pressureConstraint) {
        if (additionalContributions || !policies.empty() || !signature.conservativeState
            || !signature.momentumPredictor || !signature.pressureCorrection || signature.auxiliarySchedule
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
            || equations.state.at("U").derivation!=StateDerivation::Velocity
            || equations.state.at("p").derivation!=StateDerivation::Pressure)
            return "conservative pressure requires velocity and pressure views of the EOS-backed physical state";
    } else if (selected=="flow.conservative" && plan.compiledProgram.hasTemporalRoot) {
        const std::vector<std::string> layout{"rho","rhoU","rhoE"};
        std::vector<CompiledEquationCall> calls;
        for (const auto& call:plan.compiledProgram.steps) if (call.temporalResidual) calls.push_back(call);
        if (calls.size()!=layout.size()) return "selected fused conservative provider requires exactly rho/rhoU/rhoE occurrences";
        for (std::size_t i=0;i<layout.size();++i) {
            if (!calls[i].temporalResidual || calls[i].backendProvider!=selected
                || calls[i].source.target.symbol!=layout[i])
                return "selected fused conservative provider has incompatible occurrence/target layout";
            const auto& symbol=equations.state.at(layout[i]);
            const int components=i==1 ? 3 : 1;
            const int offset=i==0 ? 0 : i==1 ? 1 : 4;
            if (symbol.storageBinding!=StorageBinding::PackedDistributed || symbol.storageKey!="conservative"
                || symbol.components!=components || symbol.componentOffset!=offset)
                return "selected fused conservative provider has incompatible physical storage: "+layout[i];
        }
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
        const std::vector<LegacyExecutionPolicy>& policies,
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
    result.auxiliarySchedule = std::any_of(
        system.legacyEquations.begin(),system.legacyEquations.end(),
        [](const EquationDescriptor& item) {
            return item.id != "continuity" && item.id != "momentum"
                && item.id != "energy" && item.id != "momentum"
                && item.category == EquationCategory::PhysicalEquation;
        });
    const auto pressure = std::find_if(
        policies.begin(),policies.end(),[](const LegacyExecutionPolicy& item) {
            return item.id == kPressureScheduleId;
        });
    if (pressure != policies.end()) {
        result.outerCorrectors = pressure->repeatCount;
        result.pressureCorrectors = pressure->nestedRepeatCount;
        result.nonOrthogonalCorrectors = pressure->innerRepeatCount-1;
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
        const std::vector<LegacyExecutionPolicy>& policies,
        bool additionalContributions) {
    const auto signature=compileExecutionCapabilities(equations,policies,&plan);
    const ProviderCatalog catalog=ProviderCatalog::builtIn();
    std::vector<ResolvedOperationBinding> bindings;
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
                if (leaf.legacyAdapter && leaf.provider.empty())
                    leaf.provider=Legacy::selectOperationProvider(*declaration,equations,numerics,policies);
                binding=catalog.resolve(*declaration,leaf.provider);
                // STATE can reject a selected implementation; it cannot choose another.
                if (binding.status==BindingStatus::Resolved && !leaf.legacyAdapter) {
                    const auto reason=validateNativeProvider(binding.provider,equations,numerics,plan,
                        signature,policies,additionalContributions);
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
