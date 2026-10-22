#include "core/system/SF_operationIds.h"
/// @file SF_pressureCoupling.cpp
/// @brief pressure-constraint coupling preset 的匹配与贡献实现。

#include "SF_pressureCoupling.h"
#include "SF_equationContribution.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {
namespace {

/// @brief 约束声明的 pressure 未知量是否存在（role 不限：Derived closure 与
///        Multiplier 都是合法压力表示）。
const ConstraintDescriptor* pressureConstraintOf(
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs=nullptr) {
    for (const ConstraintDescriptor& constraint : raw.constraints) {
        if (inputs ? constraint.id==inputs->constraint : (constraint.id == "C_INCOMPRESSIBILITY"
            || constraint.id == "C_SHARED_PRESSURE")) {
            return &constraint;
        }
    }
    return nullptr;
}

bool rawHasUnknown(const RawEquationSystem& raw, const std::string& id) {
    return std::any_of(
        raw.state.symbols().begin(),raw.state.symbols().end(),
        [&](const StateSymbol& unknown) { return unknown.id == id; });
}

bool rawHasEquation(const RawEquationSystem& raw,const std::string& id) {
    return raw.registry.contains(id);
}

bool hasMomentumEquation(const RawEquationSystem& raw) {
    return std::any_of(raw.registry.entries().begin(),raw.registry.entries().end(),
        [&](const Equation& equation) {
            const auto vectorDerivative=[&](const auto& self,const FormulaExpr& expression)->bool {
                if (expression.kind()==FormulaExpr::Kind::Operator && expression.operation()==Math::Op::Ddt
                    && expression.arguments().size()==1 && expression.arguments()[0].kind()==FormulaExpr::Kind::Symbol) {
                    const auto& id=expression.arguments()[0].name();
                    return raw.state.contains(id) && raw.state.at(id).components>1;
                }
                return std::any_of(expression.arguments().begin(),expression.arguments().end(),[&](const auto& child){return self(self,child);});
            };
            return vectorDerivative(vectorDerivative,equation.lhs);
        });
}

} // namespace

std::vector<Equation> conservativePressureRelations() {
    using Expr=FormulaExpr;
    const auto symbol=[](const char* id) { return PressureMathematics::symbol(id,true); };
    const Provenance origin{OriginKind::Generated,"conservative pressure correction"};
    // Corrector::assemble uses reaction + face-averaged inverse-density diffusion.
    const auto compliance=Expr::multiply(symbol("rho"),Expr::multiply(
        Expr::multiply(symbol("soundSpeed"),symbol("soundSpeed")),
        Expr::multiply(symbol("dt"),symbol("dt"))));
    std::vector<Equation> equations{
        {"pSimple",Expr::subtract(Expr::divide(symbol("pPrime"),compliance),
            Expr::op(Math::Op::Div,{Expr::multiply(
                PressureMathematics::op("faceAverage",{Expr::divide(Expr::constantValue(1.0),symbol("rho"))}),
                Expr::op(Math::Op::Grad,{symbol("pPrime")}))})),
            Expr::negate(Expr::divide(Expr::op(Math::Op::Div,{symbol("U")}),symbol("dt"))),origin},
        {"correctU",symbol("rhoU"),Expr::subtract(symbol("rhoU"),
            Expr::multiply(Expr::multiply(symbol("momentumRelaxation"),symbol("dt")),
                Expr::op(Math::Op::Grad,{symbol("pPrime")}))),origin},
        {"correctP",symbol("p"),Expr::add(symbol("p"),
            Expr::multiply(symbol("pressureRelaxation"),symbol("pPrime"))),origin},
        {"correctFluxp",symbol("fluxValidity"),
            PressureMathematics::op("invalidateDerivedFlux",{symbol("rhoU")}),origin},
        {"publishPressure",symbol("rhoE"),PressureMathematics::op("totalEnergyFromPressure",
            {symbol("rho"),symbol("rhoU"),symbol("p")}),origin},
        {"relaxIterate",symbol("iterate"),PressureMathematics::op("relax",
            {symbol("iterate"),symbol("laggedIterate")}),origin},
        {"restoreFlux",symbol("fluxValidity"),PressureMathematics::op("consistentFlux",
            {symbol("rhoU"),symbol("p")}),origin},
        {"checkConvergence",symbol("converged"),PressureMathematics::op("residualConvergence",
            {symbol("iterate"),symbol("laggedIterate")}),origin}};
    return equations;
}

const char* toString(CouplingStatus status) {
    switch (status) {
        case CouplingStatus::Active: return "active";
        case CouplingStatus::Inactive: return "inactive";
        case CouplingStatus::Invalid: return "invalid";
        case CouplingStatus::Unsupported: return "unsupported";
    }
    return "inactive";
}

CouplingPresetRequest couplingRequestFrom(
        const FDM::PressureCouplingConfig& coupling,
        bool explicitlyRegistered) {
    CouplingPresetRequest request;
    request.presetKind = coupling.preset;
    request.preset = FDM::toString(coupling.preset);
    request.outerCorrectors = coupling.outerCorrectors;
    request.pressureCorrectors = coupling.pressureCorrectors;
    request.nonOrthogonalCorrectors = coupling.nonOrthogonalCorrectors;
    request.explicitlyRegistered = explicitlyRegistered;
    request.origin = {OriginKind::BuiltinPreset,"coupling preset"};
    return request;
}

PressureCouplingInputs pressureCouplingInputs(EquationCall predictor,PressurePredictorKind representation) {
    const bool conservative=representation==PressurePredictorKind::ConservativeMomentum;
    PressureCouplingInputs input;
    input.predictor=std::move(predictor);input.representation=representation;
    input.predictorTarget={input.predictor.target.symbol,conservative ? TargetKind::Physical : TargetKind::Working};
    input.predictorOccurrence="predictor";input.constraint="C_INCOMPRESSIBILITY";
    input.correction={"pSimple",{"p",TargetKind::Correction}};
    input.velocityUpdate={"correctU",{input.predictor.target.symbol}};
    input.pressureUpdate={"correctP",{"p",conservative ? TargetKind::Working : TargetKind::Physical}};
    input.fluxUpdate={"correctFluxp",{conservative ? "fluxValidity" : "phi",conservative ? TargetKind::Workspace : TargetKind::Physical}};
    if (conservative) {
        input.fusedBalances={{"continuity",{"rho"}},{"energy",{"rhoE"}}};
        input.replacedCommit="physicalStep.commit";
        input.publication={"publishPressure",{"rhoE"}};
    }
    input.relaxation={"relaxIterate",{"iterate",TargetKind::Workspace}};
    input.fluxConsistency={"restoreFlux",{conservative ? "fluxValidity" : "phi",TargetKind::Workspace}};
    input.convergence={"checkConvergence",{"converged",TargetKind::Workspace}};
    return input;
}

void applyPressureExecution(ExecutionProgram& program,
                            const CouplingPresetRequest& request,
                            const PressureCouplingInputs& input) {
    auto& entries=program.root.children;
    const auto matches=[](const ExecutionScope& node,const EquationCall& call) {
        return node.kind==ExecutionKind::EquationCall && node.step.equation==call.equation
            && node.step.target.symbol==call.target.symbol && node.step.target.kind==call.target.kind
            && node.step.occurrence==call.occurrence;
    };
    if (input.predictorOccurrence.empty() || input.predictorTarget.symbol!=input.predictor.target.symbol)
        throw std::runtime_error("Pressure method requires an explicit predictor occurrence and matching target.");
    const auto count=std::count_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,input.predictor);});
    if (count!=1) throw std::runtime_error("Pressure HOW requires one referenced predictor: "+input.predictor.equation+" / "+input.predictor.occurrence);
    for (const auto& balance:input.fusedBalances)
        if (std::count_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,balance);})!=1)
            throw std::runtime_error("Pressure HOW requires one referenced fused balance: "+balance.equation);
    const auto found=std::find_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,input.predictor);});
    ExecutionScope predictor=*found;
    entries.erase(found);
    predictor.step.target=input.predictorTarget;
    predictor.step.occurrence=input.predictorOccurrence;
    predictor.order=10;
    predictor.origin={OriginKind::Generated,request.preset+" transforms "+input.predictor.equation+" occurrence="+input.predictor.occurrence+" -> "+input.predictor.target.symbol+" @ 10; source="+predictor.origin.source};
    entries.erase(std::remove_if(entries.begin(),entries.end(),[&](const auto& node) {
        return (!input.replacedCommit.empty() && node.kind==ExecutionKind::Commit && node.id==input.replacedCommit)
            || std::any_of(input.fusedBalances.begin(),input.fusedBalances.end(),[&](const auto& call){return matches(node,call);});
    }),entries.end());
    const auto call=[&](const EquationCall& reference,Order order) {
        if (reference.equation.empty() || reference.target.symbol.empty())
            throw std::runtime_error("Pressure contribution has an empty equation/target reference.");
        ExecutionScope node;
        node.kind=ExecutionKind::EquationCall;
        node.order=order;node.step=reference;
        node.origin={OriginKind::BuiltinPreset,request.preset};
        return node;
    };
    const auto loop=[&](std::string id,int count,Order order) {
        ExecutionScope node;
        node.kind=ExecutionKind::Loop;
        node.id=std::move(id);
        node.repetitions=count;
        node.order=order;
        node.origin={OriginKind::BuiltinPreset,request.preset};
        return node;
    };
    const auto correction=[&] {
        ExecutionScope body;
        body.id="pressureCorrection";
        auto nonOrthogonal=loop("nonOrthogonal",request.nonOrthogonalCorrectors+1,10);
        nonOrthogonal.children.push_back(call(input.correction,10));
        body.children.push_back(std::move(nonOrthogonal));
        body.children.push_back(call(input.velocityUpdate,20));
        body.children.push_back(call(input.pressureUpdate,30));
        body.children.push_back(call(input.fluxUpdate,40));
        if (!input.publication.equation.empty())
            body.children.push_back(call(input.publication,50));
        body.origin={OriginKind::BuiltinPreset,request.preset};
        return body;
    };
    if (request.presetKind==FDM::PressureCouplingPreset::PISO) {
        entries.push_back(std::move(predictor));
        auto pressure=loop("pressure",request.pressureCorrectors,20);
        pressure.children.push_back(correction());
        entries.push_back(std::move(pressure));
    } else {
        auto outer=loop("outer",request.outerCorrectors,20);
        outer.terminationSignal=kPressureOuterConvergedSignal;
        outer.children.push_back(std::move(predictor));
        if (request.presetKind==FDM::PressureCouplingPreset::SIMPLE) {
            auto body=correction();body.order=20;
            outer.children.push_back(std::move(body));
        } else {
            auto pressure=loop("pressure",request.pressureCorrectors,20);
            pressure.children.push_back(correction());
            outer.children.push_back(std::move(pressure));
        }
        outer.children.push_back(call(input.relaxation,30));
        outer.children.push_back(call(input.fluxConsistency,40));
        outer.children.push_back(call(input.convergence,50));
        entries.push_back(std::move(outer));
    }
    ExecutionScope commit;
    commit.kind=ExecutionKind::Commit;
    commit.id="physicalStep.commit";
    commit.order=1000000;
    commit.origin={OriginKind::BuiltinPreset,request.preset};
    entries.push_back(std::move(commit));

}

std::vector<NumericalBinding> pressureNumerics(const CouplingPresetRequest& request,
        const PressureCouplingInputs& input) {
    const bool conservative=input.representation==PressurePredictorKind::ConservativeMomentum;
    const auto binding=[&](const EquationCall& call,const char* method) {
        return NumericalBinding{call.equation,method,{},call.occurrence};
    };
    auto predictor=binding(input.predictor,conservative ? "ConservativePressureMomentum" : "PressureMomentum");
    predictor.occurrence=input.predictorOccurrence;
    for (const auto& balance:input.fusedBalances) predictor.inputs.push_back(balance.equation);
    std::vector<NumericalBinding> bindings{predictor,
        binding(input.correction,conservative ? "ConservativePressureCorrection" : "PressureCorrection"),
        binding(input.velocityUpdate,conservative ? "ConservativeVelocityCorrection" : "VelocityCorrection"),
        binding(input.pressureUpdate,conservative ? "ConservativePressureUpdate" : "PressureUpdate"),
        binding(input.fluxUpdate,conservative ? "ConservativeFluxCorrection" : "FluxCorrection")};
    if (!input.publication.equation.empty()) bindings.push_back(binding(input.publication,"ConservativePressurePublication"));
    if (request.presetKind!=FDM::PressureCouplingPreset::PISO) {
        bindings.push_back(binding(input.relaxation,conservative ? "ConservativeFixedTimeRelaxation" : "FixedTimeRelaxation"));
        bindings.push_back(binding(input.fluxConsistency,conservative ? "ConservativeFluxConsistency" : "FluxConsistency"));
        bindings.push_back(binding(input.convergence,conservative ? "ConservativeResidualConvergence" : "ResidualConvergence"));
    }
    return bindings;
}

CouplingReport matchPressureCoupling(
        const CouplingPresetRequest& request,
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs) {
    CouplingReport report;
    report.id = "coupling."+request.preset;
    report.preset = request.preset;
    report.requirements = {
        "momentum equation",
        "incompressibility / pressure-multiplier constraint",
        "pressure multiplier unknown"};

    if (inputs ? !rawHasEquation(raw,inputs->predictor.equation) : !hasMomentumEquation(raw)) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "momentum equation not present in the resolved equation system";
        return report;
    }
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw,inputs);
    if (!pressureConstraint) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "incompressibility / pressure-multiplier constraint not present";
        return report;
    }
    const std::string multiplier = pressureConstraint->multiplierUnknown;
    const bool pressureBound = multiplier.empty()
        ? rawHasUnknown(raw,"p")
        : rawHasUnknown(raw,multiplier);
    if (!pressureBound) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "the pressure constraint is not bound to a declared pressure "
            "unknown";
        return report;
    }
    // pressure-constraint formulation 需要生成算法压力方程；它由 transformer
    // 产生（pSimple / E_SHARED_PRESSURE），不是用户输入。
    const bool shared =
        pressureConstraint->id == "C_SHARED_PRESSURE";
    const auto pressureEquation=inputs ? inputs->correction.equation : (shared ? "E_SHARED_PRESSURE" : "pSimple");
    const bool hasPressureEquation=rawHasEquation(raw,pressureEquation);
    if (!hasPressureEquation) {
        report.derivedEquations = {
            pressureEquation};
    }
    if (request.outerCorrectors <= 0 || request.pressureCorrectors <= 0
        || request.nonOrthogonalCorrectors < 0) {
        report.status = CouplingStatus::Invalid;
        report.reason =
            "registered coupling preset has invalid corrector counts";
        return report;
    }
    if (request.presetKind == FDM::PressureCouplingPreset::PISO
        && request.outerCorrectors != 1) {
        report.status = CouplingStatus::Invalid;
        report.reason = "PISO does not define an outer fixed-point iteration; "
            "use PIMPLE for outerCorrectors > 1.";
        return report;
    }
    if (request.presetKind == FDM::PressureCouplingPreset::SIMPLE
        && request.pressureCorrectors != 1) {
        report.status = CouplingStatus::Invalid;
        report.reason = "SIMPLE uses one pressure correction per outer "
            "iteration; use PIMPLE for multiple inner pressure corrections.";
        return report;
    }
    // Final operations are reported from the compiled owned plan.
    report.status = CouplingStatus::Active;
    report.reason = "resolved equation/constraint system satisfies the "
                    "pressure-coupling requirements";
    return report;
}

CouplingReport contributePressureCoupling(
        SystemCompositionBuilder& system,
        const CouplingPresetRequest& request,
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs) {
    CouplingReport report = matchPressureCoupling(request,raw,inputs);
    if (report.status != CouplingStatus::Active) return report;

    // 1. formulation / transformation contribution。
    //    Equation source 可能已经声明了同一 formulation（例如 Eulerian 模板
    //    本身要求 shared-pressure formulation），此时不重复请求。
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw,inputs);
    const bool shared = pressureConstraint
        && pressureConstraint->id == "C_SHARED_PRESSURE";
    if (shared) {
        if (!system.requestsTransformation("sharedPressureConstraint")) {
            system.requestTransformation({
                "sharedPressureConstraint",
                request.preset+" shared-pressure transformation",
                110,true,request.origin});
        }
    } else if (!system.requestsTransformation("pressureConstraint")) {
        system.requestTransformation({
            "pressureConstraint",
            request.preset+" pressure transformation",
            100,true,request.origin});
    }

    // Coupling counts are authored only in native HOW.
    return report;
}

} // namespace SF::System
