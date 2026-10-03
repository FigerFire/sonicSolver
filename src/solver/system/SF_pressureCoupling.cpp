#include "core/system/SF_operationIds.h"
/// @file SF_pressureCoupling.cpp
/// @brief pressure-constraint coupling preset 的匹配与贡献实现。

#include "SF_pressureCoupling.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {
namespace {

bool rawHasConstraint(const RawEquationSystem& raw, const std::string& id) {
    return std::any_of(
        raw.constraints.begin(),raw.constraints.end(),
        [&](const ConstraintDescriptor& constraint) {
            return constraint.id == id;
        });
}

/// @brief 约束声明的 pressure 未知量是否存在（role 不限：Derived closure 与
///        Multiplier 都是合法压力表示）。
const ConstraintDescriptor* pressureConstraintOf(
        const RawEquationSystem& raw) {
    for (const ConstraintDescriptor& constraint : raw.constraints) {
        if (constraint.id == "C_INCOMPRESSIBILITY"
            || constraint.id == "C_SHARED_PRESSURE") {
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
        [](const Equation& equation) { return equation.id.rfind("momentum",0)==0; });
}

} // namespace

std::vector<Equation> conservativePressureRelations() {
    using Expr=FormulaExpr;
    const auto symbol=[](const char* id) { return Expr::symbol(id); };
    const Provenance origin{OriginKind::Generated,"conservative pressure correction"};
    // Corrector::assemble uses reaction + face-averaged inverse-density diffusion.
    const auto compliance=Expr::multiply(symbol("rho"),Expr::multiply(
        Expr::multiply(symbol("soundSpeed"),symbol("soundSpeed")),
        Expr::multiply(symbol("dt"),symbol("dt"))));
    std::vector<Equation> equations{
        {"pSimple",Expr::subtract(Expr::divide(symbol("pPrime"),compliance),
            Expr::op("div",{Expr::multiply(
                Expr::op("faceAverage",{Expr::divide(Expr::constantValue(1.0),symbol("rho"))}),
                Expr::op("grad",{symbol("pPrime")}))})),
            Expr::negate(Expr::divide(Expr::op("div",{symbol("U")}),symbol("dt"))),origin},
        {"correctU",symbol("rhoU"),Expr::subtract(symbol("rhoU"),
            Expr::multiply(Expr::multiply(symbol("momentumRelaxation"),symbol("dt")),
                Expr::op("grad",{symbol("pPrime")}))),origin},
        {"correctP",symbol("p"),Expr::add(symbol("p"),
            Expr::multiply(symbol("pressureRelaxation"),symbol("pPrime"))),origin},
        {"correctFluxp",symbol("fluxValidity"),
            Expr::op("invalidateDerivedFlux",{symbol("rhoU")}),origin},
        {"publishPressure",symbol("rhoE"),Expr::op("totalEnergyFromPressure",
            {symbol("rho"),symbol("rhoU"),symbol("p")}),origin},
        {"relaxIterate",symbol("iterate"),Expr::op("relax",
            {symbol("iterate"),symbol("laggedIterate")}),origin},
        {"restoreFlux",symbol("fluxValidity"),Expr::op("consistentFlux",
            {symbol("rhoU"),symbol("p")}),origin},
        {"checkConvergence",symbol("converged"),Expr::op("residualConvergence",
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

void applyPressureExecution(ExecutionProgram& program,
                            const CouplingPresetRequest& request,
                            std::string_view momentumTarget) {
    const bool conservative=momentumTarget=="rhoU";
    auto& entries=program.root.children;
    const auto found=std::find_if(entries.begin(),entries.end(),
        [&](const ExecutionScope& node) {
            return node.kind==ExecutionKind::EquationCall
                && node.step.equation=="momentum"
                && node.step.target.symbol==momentumTarget && node.step.target.kind==TargetKind::Physical
                && (node.id=="momentum.default" || node.step.occurrence.empty());
        });
    if (found==entries.end())
        throw std::runtime_error("Pressure HOW transformation requires a default momentum occurrence.");
    ExecutionScope predictor=*found;
    entries.erase(found);
    predictor.step.target.kind=conservative ? TargetKind::Physical : TargetKind::Working;
    predictor.step.occurrence="predictor";
    predictor.order=10;
    predictor.origin={OriginKind::Generated,request.preset+" transforms default momentum -> "+std::string(momentumTarget)+" @ 10"};
    if (conservative) {
        entries.erase(std::remove_if(entries.begin(),entries.end(),[](const auto& node) {
            return node.origin.kind==OriginKind::BuiltinPreset && node.origin.source=="NavierStokes"
                && (node.kind==ExecutionKind::Commit
                    || (node.kind==ExecutionKind::EquationCall && node.step.occurrence.empty()
                        && (node.step.equation=="continuity" || node.step.equation=="energy")));
        }),entries.end());
    }
    const auto call=[&](const char* equation,const char* symbol,Order order,TargetKind kind) {
        ExecutionScope node;
        node.kind=ExecutionKind::EquationCall;
        node.order=order;
        node.step={equation,{symbol,kind}};
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
        nonOrthogonal.children.push_back(call("pSimple","p",10,
            TargetKind::Correction));
        body.children.push_back(std::move(nonOrthogonal));
        body.children.push_back(call("correctU",conservative ? "rhoU" : "U",20,TargetKind::Physical));
        body.children.push_back(call("correctP","p",30,conservative ? TargetKind::Working : TargetKind::Physical));
        body.children.push_back(call("correctFluxp",conservative ? "fluxValidity" : "phi",40,
            conservative ? TargetKind::Workspace : TargetKind::Physical));
        if (conservative)
            body.children.push_back(call("publishPressure","rhoE",50,TargetKind::Physical));
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
        outer.children.push_back(call("relaxIterate","iterate",30,
            TargetKind::Workspace));
        outer.children.push_back(call("restoreFlux",conservative ? "fluxValidity" : "phi",40,
            TargetKind::Workspace));
        outer.children.push_back(call("checkConvergence","converged",50,
            TargetKind::Workspace));
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
        std::string_view momentumTarget) {
    if (momentumTarget=="rhoU") {
        std::vector<NumericalBinding> bindings{
            {"momentum","ConservativePressureMomentum",{"continuity","energy"},"predictor"},
            {"pSimple","ConservativePressureCorrection"},
            {"correctU","ConservativeVelocityCorrection"},
            {"correctP","ConservativePressureUpdate"},
            {"correctFluxp","ConservativeFluxCorrection"},
            {"publishPressure","ConservativePressurePublication"}};
        if (request.presetKind!=FDM::PressureCouplingPreset::PISO) {
            bindings.insert(bindings.end(),{{"relaxIterate","ConservativeFixedTimeRelaxation"},
                {"restoreFlux","ConservativeFluxConsistency"},
                {"checkConvergence","ConservativeResidualConvergence"}});
        }
        return bindings;
    }
    std::vector<NumericalBinding> bindings{
        {"momentum","PressureMomentum",{},"predictor"},
        {"pSimple","PressureCorrection"},{"correctU","VelocityCorrection"},
        {"correctP","PressureUpdate"},{"correctFluxp","FluxCorrection"}};
    if (request.presetKind!=FDM::PressureCouplingPreset::PISO) {
        bindings.insert(bindings.end(),{{"relaxIterate","FixedTimeRelaxation"},
            {"restoreFlux","FluxConsistency"},{"checkConvergence","ResidualConvergence"}});
    }
    return bindings;
}

/// @brief PISO/分离式压力修正的执行片段。
///
/// 片段只表达顺序与重复；每个 Leaf 引用 formulation 的 OperationStage。
/// 具体 OpId 由 executable operation authority 解析，因此这里不再维护任何
/// operation/provider 名单。


namespace Legacy {


} // namespace Legacy

namespace {

void collectMissingStages(const PlanFragmentNode& node,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing) {
    if (node.kind == PlanFragmentNode::Kind::Leaf) {
        // Exact EquationCall resolution happens after EquationMethod compilation.
        if (!node.step.empty()) return;
        if (!node.operation.empty()) {
            if (std::any_of(executable.operations.begin(),executable.operations.end(),
                [&](const ExecutableOperation& item) {
                    return item.operation == node.operation;
                })) return;
            if (std::find(missing->begin(),missing->end(),node.operation)
                == missing->end()) missing->push_back(node.operation);
            return;
        }
        // 与 planner 的 lowering 规则一致：带具名缺失标记的 leaf 永远
        // unresolved；否则按 stage 解析由 formulation 声明的 operation。
        if (node.missingOperation.empty()
            && findExecutableOperation(executable,node.stage)) {
            return;
        }
        const std::string id = node.missingOperation.empty()
            ? node.id : node.missingOperation;
        if (std::find(missing->begin(),missing->end(),id) == missing->end()) {
            missing->push_back(id);
        }
        return;
    }
    for (const PlanFragmentNode& child : node.children) {
        collectMissingStages(child,executable,missing);
    }
}

} // namespace

bool couplingPlanResolved(const LegacyPlanFragment& fragment,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing) {
    std::vector<std::string> found;
    for (const PlanFragmentNode& node : fragment.nodes) {
        collectMissingStages(node,executable,&found);
    }
    if (missing) *missing = found;
    return found.empty();
}



CouplingReport matchPressureCoupling(
        const CouplingPresetRequest& request,
        const RawEquationSystem& raw) {
    CouplingReport report;
    report.id = "coupling."+request.preset;
    report.preset = request.preset;
    report.requirements = {
        "momentum equation",
        "incompressibility / pressure-multiplier constraint",
        "pressure multiplier unknown"};

    if (!hasMomentumEquation(raw)) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "momentum equation not present in the resolved equation system";
        return report;
    }
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw);
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
    const bool hasPressureEquation =
        rawHasEquation(raw,"pSimple") || rawHasEquation(raw,"E_SHARED_PRESSURE");
    if (!hasPressureEquation) {
        report.derivedEquations = {
            shared ? "E_SHARED_PRESSURE" : "pSimple"};
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
        const RawEquationSystem& raw) {
    CouplingReport report = matchPressureCoupling(request,raw);
    if (report.status != CouplingStatus::Active) return report;

    // 1. formulation / transformation contribution。
    //    Equation source 可能已经声明了同一 formulation（例如 Eulerian 模板
    //    本身要求 shared-pressure formulation），此时不重复请求。
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw);
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
