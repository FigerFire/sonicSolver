/// @file SF_pressureCoupling.cpp
/// @brief pressure-constraint coupling preset 的匹配与贡献实现。

#include "SF_pressureCoupling.h"

#include <algorithm>

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
        raw.unknowns.begin(),raw.unknowns.end(),
        [&](const UnknownDescriptor& unknown) { return unknown.id == id; });
}

bool rawHasEquation(const RawEquationSystem& raw, const std::string& id) {
    return std::any_of(
        raw.equations.begin(),raw.equations.end(),
        [&](const EquationDescriptor& equation) {
            return equation.id == id;
        });
}

/// @brief 是否存在动量方程（单流体 E_MOMENTUM 或逐相 E_MOMENTUM.<phase>）。
bool hasMomentumEquation(const RawEquationSystem& raw) {
    return std::any_of(
        raw.equations.begin(),raw.equations.end(),
        [](const EquationDescriptor& equation) {
            return equation.id.compare(0,10,"E_MOMENTUM") == 0;
        });
}

} // namespace

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

/// @brief PISO/分离式压力修正的执行片段。
///
/// 片段只表达顺序与重复；每个 Leaf 引用 formulation 的 OperationStage。
/// 具体 OpId 由 executable operation authority 解析，因此这里不再维护任何
/// operation/provider 名单。
PlanFragment couplingPlanFragment(const CouplingPresetRequest& request) {
    const auto leaf = [](PlanNodeKind kind, std::string id, std::string name,
                         OperationStage stage, std::string equation = {}) {
        PlanFragmentNode node;
        node.kind = PlanFragmentNode::Kind::Leaf;
        node.id = std::move(id);
        node.name = std::move(name);
        node.leafKind = kind;
        node.stage = stage;
        node.equation = std::move(equation);
        return node;
    };
    const auto loop = [](std::string id, std::string name, int repetitions) {
        PlanFragmentNode node;
        node.kind = PlanFragmentNode::Kind::Loop;
        node.id = std::move(id);
        node.name = std::move(name);
        node.nodeKind = PlanNodeKind::Loop;
        node.repetitions = repetitions;
        return node;
    };
    const auto sequence = [](std::string id, std::string name) {
        PlanFragmentNode node;
        node.kind = PlanFragmentNode::Kind::Sequence;
        node.id = std::move(id);
        node.name = std::move(name);
        node.nodeKind = PlanNodeKind::Sequence;
        return node;
    };

    PlanFragment fragment;
    fragment.id = request.preset+" pressure schedule";
    fragment.consumedPolicies = {kPressureScheduleId};
    fragment.priority = 20;
    if (request.presetKind != FDM::PressureCouplingPreset::PISO) {
        // SIMPLE/PIMPLE fixed point 需要 dedicated predictor provider；当前
        // 没有实现，因此片段引用一个具名缺失 operation，planner 保留真实的
        // control flow，provider resolver 报告 Unsupported。
        const auto missing = [&](PlanNodeKind kind, std::string id,
                                 std::string name) {
            PlanFragmentNode node = leaf(kind,std::move(id),std::move(name),
                                        OperationStage::Prepare);
            node.missingOperation = "pressure.schedule."+node.id;
            node.unsupportedReason =
                "The selected SIMPLE/PIMPLE fixed-point schedule has no "
                "dedicated predictor provider. The existing momentum.solve "
                "operation advances physical state with the full dt and "
                "cannot be repeated as an outer corrector.";
            node.id = "PressureSchedule."+node.id;
            return node;
        };
        fragment.name = "single-fluid pressure schedule without a provider";
        fragment.nodes.push_back(missing(
            PlanNodeKind::Update,"prepare","prepare pressure schedule"));
        PlanFragmentNode outer = loop(
            "PressureSchedule.outerCorrectors",
            "pressure-velocity outer correctors",request.outerCorrectors);
        outer.children.push_back(missing(
            PlanNodeKind::Assemble,"predictor.assemble",
            "assemble dedicated fixed-point predictor"));
        outer.children.push_back(missing(
            PlanNodeKind::Solve,"predictor.solve",
            "solve dedicated fixed-point predictor"));
        PlanFragmentNode corrections = loop(
            "PressureSchedule.pressureCorrectors","pressure correctors",
            request.pressureCorrectors);
        corrections.children.push_back(missing(
            PlanNodeKind::Update,"boundary.prepare",
            "prepare pressure boundary state"));
        PlanFragmentNode nonOrthogonal = loop(
            "PressureSchedule.nonOrthogonalCorrectors",
            "non-orthogonal passes",request.nonOrthogonalCorrectors+1);
        nonOrthogonal.children.push_back(missing(
            PlanNodeKind::Assemble,"assemble",
            "assemble pressure correction"));
        nonOrthogonal.children.push_back(missing(
            PlanNodeKind::Solve,"solve","solve pressure correction"));
        corrections.children.push_back(std::move(nonOrthogonal));
        corrections.children.push_back(missing(
            PlanNodeKind::Correct,"update.prepare",
            "prepare pressure update"));
        corrections.children.push_back(missing(
            PlanNodeKind::Correct,"velocity","correct velocity"));
        corrections.children.push_back(missing(
            PlanNodeKind::Correct,"flux","correct face flux"));
        corrections.children.push_back(missing(
            PlanNodeKind::Update,"correction.commit",
            "commit pressure correction"));
        outer.children.push_back(std::move(corrections));
        outer.children.push_back(missing(
            PlanNodeKind::ConvergenceCheck,"convergence.check",
            "check pressure-velocity convergence"));
        fragment.nodes.push_back(std::move(outer));
        fragment.nodes.push_back(missing(
            PlanNodeKind::Commit,"step.commit","commit pressure step"));
        return fragment;
    }

    fragment.name = "generic PISO time-step sequence";
    fragment.nodes.push_back(leaf(
        PlanNodeKind::Update,"PISO.prepare","prepare boundary and closure",
        OperationStage::Prepare));
    PlanFragmentNode outer = loop("PISO.outerCorrectors","PISO outer schedule",
                                  request.outerCorrectors);
    outer.children.push_back(leaf(
        PlanNodeKind::Assemble,"PISO.momentum.assemble",
        "assemble momentum predictor",OperationStage::MomentumAssemble,
        "E_MOMENTUM_PREDICTOR"));
    outer.children.push_back(leaf(
        PlanNodeKind::Solve,"PISO.momentum.solve","advance momentum predictor",
        OperationStage::MomentumSolve,"E_MOMENTUM_PREDICTOR"));
    PlanFragmentNode correctorLoop = loop(
        "PISO.pressureCorrectors","PISO pressure correctors",
        request.pressureCorrectors);
    PlanFragmentNode correction = sequence(
        "PISO.pressureCorrection","pressure correction sequence");
    correction.children.push_back(leaf(
        PlanNodeKind::Update,"PISO.pressure.prepare",
        "prepare pressure boundary state",
        OperationStage::PressureBoundaryPrepare));
    PlanFragmentNode nonOrthogonal = loop(
        "PISO.nonOrthogonalCorrectors","PISO non-orthogonal passes",
        request.nonOrthogonalCorrectors+1);
    nonOrthogonal.children.push_back(leaf(
        PlanNodeKind::Assemble,"PISO.pressure.assemble",
        "assemble pressure correction",OperationStage::PressureAssemble,
        "E_PRESSURE"));
    nonOrthogonal.children.push_back(leaf(
        PlanNodeKind::Solve,"PISO.pressure.solve","solve pressure correction",
        OperationStage::PressureSolve,"E_PRESSURE"));
    correction.children.push_back(std::move(nonOrthogonal));
    correction.children.push_back(leaf(
        PlanNodeKind::Correct,"PISO.velocity.correct","velocity correction",
        OperationStage::VelocityCorrect));
    correction.children.push_back(leaf(
        PlanNodeKind::Correct,"PISO.pressure.update","prepare pressure update",
        OperationStage::PressureUpdatePrepare));
    correction.children.push_back(leaf(
        PlanNodeKind::Correct,"PISO.flux.correct",
        "refresh derived face-flux state",OperationStage::FluxCorrect));
    correction.children.push_back(leaf(
        PlanNodeKind::Update,"PISO.pressure.commit",
        "commit pressure-corrected state",OperationStage::CorrectionCommit));
    correctorLoop.children.push_back(std::move(correction));
    outer.children.push_back(std::move(correctorLoop));
    fragment.nodes.push_back(std::move(outer));
    fragment.nodes.push_back(leaf(
        PlanNodeKind::Commit,"PISO.commit","commit corrected state",
        OperationStage::StepCommit));
    return fragment;
}

PlanFragment sharedPressurePlanFragment(
        const ExecutionPolicy& policy,
        const ExecutableEquationSystem& executable) {
    const auto leaf = [](const char* operation, const char* name) {
        PlanFragmentNode node;
        node.kind = PlanFragmentNode::Kind::Leaf;
        node.leafKind = PlanNodeKind::Update;
        node.operation = operation;
        node.id = std::string("EE.") + std::string(operation).substr(3);
        node.name = name;
        return node;
    };
    const auto loop = [](const char* id, const char* name, int count) {
        PlanFragmentNode node;
        node.kind = PlanFragmentNode::Kind::Loop;
        node.nodeKind = PlanNodeKind::Loop;
        node.id = id;
        node.name = name;
        node.repetitions = count;
        return node;
    };
    const auto declared = [&](const char* id) {
        return std::any_of(executable.operations.begin(),executable.operations.end(),
            [id](const ExecutableOperation& item) {
                return item.operation == id;
            });
    };
    PlanFragment fragment;
    fragment.id = "sharedPressure";
    fragment.rootId = "EE.step";
    fragment.consumedPolicies = {kSharedPressureScheduleId};
    fragment.includeTimePolicy = false;
    fragment.name = "Eulerian PIMPLE time step";
    fragment.nodes.push_back(leaf("ee.dt.compute","compute stable time step"));
    auto begin = leaf("ee.step.begin","reset diagnostics and prepare state");
    begin.id = "EE.begin";
    fragment.nodes.push_back(std::move(begin));
    auto outer = loop("EE.outer","Eulerian outer corrector",policy.repeatCount);
    outer.children.push_back(leaf("ee.interphase.compute","interphase coupling"));
    outer.children.push_back(leaf("ee.sources.assemble","phase sources"));
    if (declared("ee.turbulence.prepare")) {
        outer.children.push_back(leaf("ee.turbulence.prepare","turbulence closure preparation"));
    }
    outer.children.push_back(leaf("ee.sources.validate","validate source time step"));
    for (const auto& item : std::vector<std::pair<const char*,const char*>>{
            {"ee.momentum.diagonal","momentum diagonal"},
            {"ee.momentum.flux","interpolated momentum flux"},
            {"ee.faceFlux.canonical","canonical phase flux"},
            {"ee.continuity.assemble","phase continuity"},
            {"ee.boundary.prepare","boundary and halo"},
            {"ee.momentum.solve","momentum predictors"},
            {"ee.interphase.correct","semi-implicit interphase"},
            {"ee.boundary.afterMomentum","boundary and halo after momentum"},
            {"ee.diagonal.sync","momentum diagonal synchronization"},
            {"ee.momentum.flux.after","interpolated momentum flux"},
            {"ee.faceFlux.canonical.after","canonical phase flux"}}) {
        outer.children.push_back(leaf(item.first,item.second));
    }
    auto pressure = loop("EE.pressure","Eulerian pressure corrector",
                         policy.nestedRepeatCount);
    auto nonOrth = loop("EE.nonOrthogonal","Eulerian non-orthogonal correction",
                        policy.innerRepeatCount);
    for (const auto& item : std::vector<std::pair<const char*,const char*>>{
            {"ee.pressure.solve","pressure correction solve"},
            {"ee.pressure.publish","publish pressure correction"},
            {"ee.pressure.sync","pressure correction synchronization"},
            {"ee.phase.correct","phase correction"},
            {"ee.faceFlux.correct","canonical face-flux correction"},
            {"ee.boundary.afterPressure","boundary and halo after pressure"},
            {"ee.faceFlux.canonical.pressure","canonical phase flux"}}) {
        nonOrth.children.push_back(leaf(item.first,item.second));
    }
    pressure.children.push_back(std::move(nonOrth));
    outer.children.push_back(std::move(pressure));
    outer.children.push_back(leaf("ee.energy.solve","phase energy"));
    if (declared("ee.turbulence.solve")) {
        outer.children.push_back(leaf("ee.turbulence.solve","turbulence equations"));
    }
    outer.children.push_back(leaf("ee.boundary.final","boundary and halo"));
    outer.children.push_back(leaf("ee.outer.validate","outer-state validation"));
    fragment.nodes.push_back(std::move(outer));
    auto commit = leaf("ee.step.commit","commit time level and diagnostics");
    commit.id = "EE.commit";
    fragment.nodes.push_back(std::move(commit));
    fragment.nodes.push_back(leaf("ee.time.commit","commit dt, physical time and step"));
    return fragment;
}

namespace {

void collectMissingStages(const PlanFragmentNode& node,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing) {
    if (node.kind == PlanFragmentNode::Kind::Leaf) {
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

bool couplingPlanResolved(const PlanFragment& fragment,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing) {
    std::vector<std::string> found;
    for (const PlanFragmentNode& node : fragment.nodes) {
        collectMissingStages(node,executable,&found);
    }
    if (missing) *missing = found;
    return found.empty();
}

std::vector<std::string> couplingReportOperations(
        const CouplingPresetRequest& request) {
    std::vector<std::string> result;
    const PlanFragment fragment = couplingPlanFragment(request);
    const auto collect = [&](const auto& self,
                             const PlanFragmentNode& node) -> void {
        if (node.kind == PlanFragmentNode::Kind::Leaf) {
            const std::string id = node.missingOperation.empty()
                ? node.id : node.missingOperation;
            if (std::find(result.begin(),result.end(),id) == result.end()) {
                result.push_back(id);
            }
            return;
        }
        for (const PlanFragmentNode& child : node.children) {
            self(self,child);
        }
    };
    for (const PlanFragmentNode& node : fragment.nodes) collect(collect,node);
    return result;
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
        report.status = CouplingStatus::Inactive;
        report.reason =
            "momentum equation not present in the resolved equation system";
        return report;
    }
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw);
    if (!pressureConstraint) {
        report.status = CouplingStatus::Inactive;
        report.reason =
            "incompressibility / pressure-multiplier constraint not present";
        return report;
    }
    const std::string multiplier = pressureConstraint->multiplierUnknown;
    const bool pressureBound = multiplier.empty()
        ? rawHasUnknown(raw,"p")
        : rawHasUnknown(raw,multiplier);
    if (!pressureBound) {
        report.status = CouplingStatus::Inactive;
        report.reason =
            "the pressure constraint is not bound to a declared pressure "
            "unknown";
        return report;
    }
    // pressure-constraint formulation 需要生成算法压力方程；它由 transformer
    // 产生（E_PRESSURE / E_SHARED_PRESSURE），不是用户输入。
    const bool shared =
        pressureConstraint->id == "C_SHARED_PRESSURE";
    const bool hasPressureEquation =
        rawHasEquation(raw,"E_PRESSURE") || rawHasEquation(raw,"E_SHARED_PRESSURE");
    if (!hasPressureEquation) {
        report.derivedEquations = {
            shared ? "E_SHARED_PRESSURE" : "E_PRESSURE"};
    }
    if (request.outerCorrectors <= 0 || request.pressureCorrectors <= 0
        || request.nonOrthogonalCorrectors < 0) {
        report.status = CouplingStatus::Invalid;
        report.reason =
            "registered coupling preset has invalid corrector counts";
        return report;
    }
    // derived operations 由 formulation 声明；这里只报告片段引用的 stage，
    // 具体 OpId 在 fragment 解析阶段从 executable operation authority 得到。
    const PlanFragment fragment = couplingPlanFragment(request);
    for (const PlanFragmentNode& node : fragment.nodes) {
        const auto collect = [&](const auto& self,
                                 const PlanFragmentNode& item) -> void {
            if (item.kind == PlanFragmentNode::Kind::Leaf) {
                const std::string id = item.missingOperation.empty()
                    ? item.id : item.missingOperation;
                if (std::find(report.derivedOperations.begin(),
                              report.derivedOperations.end(),id)
                    == report.derivedOperations.end()) {
                    report.derivedOperations.push_back(id);
                }
                return;
            }
            for (const PlanFragmentNode& child : item.children) {
                self(self,child);
            }
        };
        collect(collect,node);
    }
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

    // 2. solve-plan fragment（外层/压力/非正交重复次数）
    ExecutionPolicy correction;
    correction.id = shared ? kSharedPressureScheduleId : kPressureScheduleId;
    correction.name = shared
        ? request.preset+" phase pressure coupling"
        : request.preset+" pressure-velocity coupling";
    const bool fixedPoint = request.presetKind != FDM::PressureCouplingPreset::PISO;
    correction.kind = fixedPoint
        ? ExecutionPolicyKind::PressureVelocityFixedPoint
        : ExecutionPolicyKind::SegregatedPressureCorrection;
    correction.strategyName = shared
        ? request.preset+" phase predictor/corrector"
        : request.preset;
    correction.strategyKind = fixedPoint
        ? FDM::SolveStrategyKind::PressureVelocityCoupling
        : FDM::SolveStrategyKind::PressureCorrection;
    if (shared) {
        // 共享压力 schedule 拥有全部相方程；缺一个都会让该方程没有 owning
        // solve block（validation 会拒绝）。
        correction.unknowns = {"p"};
        for (const EquationDescriptor& equation : raw.equations) {
            correction.equations.push_back(equation.id);
            for (const std::string& unknown : equation.solvedUnknowns) {
                if (unknown == "p") continue;
                if (std::find(correction.unknowns.begin(),
                              correction.unknowns.end(),unknown)
                    == correction.unknowns.end()) {
                    correction.unknowns.push_back(unknown);
                }
            }
        }
        // Algorithmic pressure equation 由 transformer 在 raw composition
        // 之后生成。
        correction.equations.push_back("E_SHARED_PRESSURE");
        correction.constraints = {"C_SHARED_PRESSURE","C_VOLUME_FRACTION"};
    } else {
        correction.equations = report.derivedEquations.empty()
            ? std::vector<std::string>{"E_PRESSURE"}
            : report.derivedEquations;
        correction.constraints = {"C_INCOMPRESSIBILITY"};
        correction.unknowns = {"pPrime","U"};
    }
    correction.priority = 20;
    correction.repeatCount = request.outerCorrectors;
    correction.nestedRepeatCount = request.pressureCorrectors;
    correction.innerRepeatCount = request.nonOrthogonalCorrectors+1;
    correction.origin = request.origin;
    system.addExecutionPolicy(std::move(correction));
    return report;
}

} // namespace SF::System
