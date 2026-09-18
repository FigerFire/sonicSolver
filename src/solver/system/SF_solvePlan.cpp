/// @file SF_solvePlan.cpp
/// @brief Structured solve-plan 编译及 transitional legacy block lowering。

#include "SF_solvePlan.h"
#include "SF_config.h"

#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>

namespace SF::System {
namespace {

SolvePlanNode leaf(
        PlanNodeKind kind,
        std::string id,
        std::string name,
        const ExecutionPolicy& policy,
        OpId operation = {}) {
    SolvePlanNode result;
    result.kind = kind;
    result.id = std::move(id);
    result.name = std::move(name);
    result.equations = policy.equations;
    result.constraints = policy.constraints;
    result.operation = operation;
    return result;
}

SolvePlanNode operationLeaf(
        PlanNodeKind kind,
        std::string id,
        std::string name,
        std::string equation,
        OpId operation) {
    SolvePlanNode result;
    result.kind = kind;
    result.id = std::move(id);
    result.name = std::move(name);
    if (!equation.empty()) result.equations.push_back(std::move(equation));
    result.operation = operation;
    return result;
}

ExecutionCapabilitySignature capabilitySignature(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies) {
    ExecutionCapabilitySignature result;
    result.pressureConstraint = hasConstraint(system,"C_INCOMPRESSIBILITY");
    result.constantDensity = std::any_of(
        system.unknowns.begin(),system.unknowns.end(),[](const UnknownDescriptor& unknown) {
            return unknown.id == "rho" && unknown.storageKey == "rhoConst";
        });
    result.conservativeState = std::any_of(
        system.unknowns.begin(),system.unknowns.end(),[](const UnknownDescriptor& unknown) {
            return unknown.storageBinding == StorageBinding::PackedDistributed
                && unknown.storageKey == "conservative";
        });
    result.momentumPredictor = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& equation) {
            return equation.operatorBinding == "momentum.predictor";
        });
    result.pressureCorrection = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& equation) {
            return equation.operatorBinding == "pressure.correction";
        });
    result.auxiliarySchedule = policies.size() != 2;
    const auto pressure = std::find_if(
        policies.begin(),policies.end(),[](const ExecutionPolicy& item) {
            return item.id == "S_PRESSURE";
        });
    const auto predictor = std::find_if(
        policies.begin(),policies.end(),[](const ExecutionPolicy& item) {
            return item.id == "S_PREDICTOR"
                && item.kind == ExecutionPolicyKind::ExplicitStages;
        });
    if (pressure != policies.end()) {
        result.pressureCorrectors = pressure->repeatCount;
        result.nonOrthogonalCorrectors = pressure->nestedRepeatCount-1;
    }
    if (predictor == policies.end()) result.momentumPredictor = false;
    return result;
}

bool supportsGenericPiso(
        const ExecutionCapabilitySignature& signature,
        const std::vector<ExecutionPolicy>& policies,
        std::string_view timeIntegrator) {
    const auto pressure = std::find_if(
        policies.begin(),policies.end(),[](const ExecutionPolicy& item) {
            return item.id == "S_PRESSURE";
        });
    return timeIntegrator == "Euler" && policies.size() == 2
        && pressure != policies.end()
        && pressure->strategyName == "PISO"
        && pressure->kind == ExecutionPolicyKind::SegregatedPressureCorrection
        && signature.pressureConstraint && signature.conservativeState
        && signature.momentumPredictor && signature.pressureCorrection
        && !signature.constantDensity && !signature.auxiliarySchedule
        && signature.pressureCorrectors > 0
        && signature.nonOrthogonalCorrectors == 0;
}

SolvePlanNode compileGenericPiso(
        const ExecutionPolicy& pressure) {
    SolvePlanNode root;
    root.kind = PlanNodeKind::Sequence;
    root.id = "PLAN_ROOT";
    root.name = "generic PISO time-step sequence";
    root.children.push_back(operationLeaf(
        PlanNodeKind::Update,"PISO.prepare","prepare boundary and closure",{},
        "pressure.prepare"));
    root.children.push_back(operationLeaf(
        PlanNodeKind::Assemble,"PISO.momentum.assemble",
        "assemble momentum predictor","E_MOMENTUM_PREDICTOR",
        "momentum.assemble"));
    root.children.push_back(operationLeaf(
        PlanNodeKind::Solve,"PISO.momentum.solve","advance momentum predictor",
        "E_MOMENTUM_PREDICTOR","momentum.solve"));

    SolvePlanNode loop;
    loop.kind = PlanNodeKind::Loop;
    loop.id = "PISO.pressureCorrectors";
    loop.name = "PISO pressure correctors";
    loop.repetitions = pressure.repeatCount;
    SolvePlanNode correction;
    correction.kind = PlanNodeKind::Sequence;
    correction.id = "PISO.pressureCorrection";
    correction.name = "pressure correction sequence";
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Update,"PISO.pressure.prepare",
        "prepare pressure boundary state",{},
        "pressure.boundary.prepare"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Assemble,"PISO.pressure.assemble",
        "assemble pressure correction","E_PRESSURE",
        "pressure.assemble"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Solve,"PISO.pressure.solve","solve pressure correction",
        "E_PRESSURE","pressure.solve"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Correct,"PISO.velocity.correct","velocity correction",{},
        "velocity.correct"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Correct,"PISO.pressure.update","prepare pressure update",{},
        "pressure.update.prepare"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Correct,"PISO.flux.correct",
        "refresh derived face-flux state",{},"flux.correct"));
    correction.children.push_back(operationLeaf(
        PlanNodeKind::Update,"PISO.pressure.commit",
        "commit pressure-corrected state",{},
        "pressure.correction.commit"));
    loop.children.push_back(std::move(correction));
    root.children.push_back(std::move(loop));
    root.children.push_back(operationLeaf(
        PlanNodeKind::Commit,"PISO.commit","commit corrected state",{},
        "pressure.step.commit"));
    return root;
}

// This is a direct lowering of EulerianEulerian::PressureStepper::stepImpl.
// Loop counts and operation order are compilation data; callbacks only invoke
// the existing numerical helpers.
SolvePlanNode compileEulerianPimple(
        const ExecutionPolicy& policy,
        bool prepareTurbulence,
        bool solveTurbulence) {
    const auto op = [](std::string id, std::string name, OpId operation) {
        return operationLeaf(PlanNodeKind::Update,std::move(id),std::move(name),{},std::move(operation));
    };
    SolvePlanNode root;
    root.kind = PlanNodeKind::Sequence;
    root.id = "EE.step";
    root.name = "Eulerian PIMPLE time step";
    root.children.push_back(op("EE.dt.compute","compute stable time step","ee.dt.compute"));
    root.children.push_back(op("EE.begin","reset diagnostics and prepare state","ee.step.begin"));
    SolvePlanNode outer;
    outer.kind = PlanNodeKind::Loop;
    outer.id = "EE.outer";
    outer.name = "Eulerian outer corrector";
    outer.repetitions = policy.repeatCount;
    for (const auto& item : std::vector<std::pair<std::string,std::string>>{
            {"interphase.compute","interphase coupling"}, {"sources.assemble","phase sources"},
            {"momentum.diagonal","momentum diagonal"},
            {"momentum.flux","interpolated momentum flux"}, {"faceFlux.canonical","canonical phase flux"},
            {"continuity.assemble","phase continuity"}, {"boundary.prepare","boundary and halo"},
            {"momentum.solve","momentum predictors"}, {"interphase.correct","semi-implicit interphase"},
            {"boundary.afterMomentum","boundary and halo after momentum"}, {"diagonal.sync","momentum diagonal synchronization"},
            {"momentum.flux.after","interpolated momentum flux"}, {"faceFlux.canonical.after","canonical phase flux"}}) {
        outer.children.push_back(op("EE."+item.first,item.second,"ee."+item.first));
    }
    if (prepareTurbulence) {
        outer.children.insert(
            outer.children.begin()+2,
            op("EE.turbulence.prepare","turbulence closure preparation",
               "ee.turbulence.prepare"));
    }
    outer.children.insert(
        outer.children.begin()+(prepareTurbulence ? 3 : 2),
        op("EE.sources.validate","validate source time step",
           "ee.sources.validate"));
    SolvePlanNode pressure;
    pressure.kind = PlanNodeKind::Loop; pressure.id = "EE.pressure";
    pressure.name = "Eulerian pressure corrector"; pressure.repetitions = policy.nestedRepeatCount;
    SolvePlanNode nonOrth;
    nonOrth.kind = PlanNodeKind::Loop; nonOrth.id = "EE.nonOrthogonal";
    nonOrth.name = "Eulerian non-orthogonal correction"; nonOrth.repetitions = policy.innerRepeatCount;
    for (const auto& item : std::vector<std::pair<std::string,std::string>>{
            {"pressure.solve","pressure correction solve"}, {"pressure.publish","publish pressure correction"},
            {"pressure.sync","pressure correction synchronization"}, {"phase.correct","phase correction"},
            {"faceFlux.correct","canonical face-flux correction"}, {"boundary.afterPressure","boundary and halo after pressure"},
            {"faceFlux.canonical.pressure","canonical phase flux"}}) {
        nonOrth.children.push_back(op("EE."+item.first,item.second,"ee."+item.first));
    }
    pressure.children.push_back(std::move(nonOrth)); outer.children.push_back(std::move(pressure));
    outer.children.push_back(op("EE.energy.solve","phase energy","ee.energy.solve"));
    if (solveTurbulence) {
        outer.children.push_back(op(
            "EE.turbulence.solve","turbulence equations","ee.turbulence.solve"));
    }
    for (const auto& item : std::vector<std::pair<std::string,std::string>>{
            {"boundary.final","boundary and halo"},
            {"outer.validate","outer-state validation"}}) {
        outer.children.push_back(op("EE."+item.first,item.second,"ee."+item.first));
    }
    root.children.push_back(std::move(outer));
    root.children.push_back(op("EE.commit","commit time level and diagnostics","ee.step.commit"));
    root.children.push_back(op("EE.time.commit","commit dt, physical time and step","ee.time.commit"));
    return root;
}

SolvePlanNode compilePolicy(
        const ExecutionPolicy& policy,
        std::string_view timeIntegrator) {
    if (policy.kind == ExecutionPolicyKind::SegregatedPressureCorrection
        || policy.kind == ExecutionPolicyKind::PressureVelocityFixedPoint) {
        SolvePlanNode loop;
        loop.kind = PlanNodeKind::Loop;
        loop.id = policy.id;
        loop.name = policy.name;
        loop.children.push_back(
            leaf(PlanNodeKind::Assemble,policy.id+".predictor","momentum predictor",policy));
        loop.children.push_back(
            leaf(PlanNodeKind::Solve,policy.id+".pressure","pressure correction solve",policy));
        loop.children.push_back(
            leaf(PlanNodeKind::Correct,policy.id+".pressureUpdate","pressure update",policy));
        loop.children.push_back(
            leaf(PlanNodeKind::Correct,policy.id+".velocity","velocity correction",policy));
        loop.children.push_back(
            leaf(PlanNodeKind::Correct,policy.id+".flux","flux correction",policy));
        loop.children.push_back(
            leaf(PlanNodeKind::ConvergenceCheck,policy.id+".convergence",
                 "convergence check",policy));
        return loop;
    }
    if (policy.kind == ExecutionPolicyKind::MonolithicKKT) {
        return leaf(
            PlanNodeKind::BlockSolve,policy.id,policy.name,policy,
            "ibm.kkt.solve");
    }
    if (policy.kind == ExecutionPolicyKind::ConstraintProjection) {
        // The existing built-in implementation performs multiplier evaluation
        // and state correction as one coherent local constraint projection.
        // One truthful leaf avoids inventing two fake callbacks around it.
        return leaf(
            PlanNodeKind::Correct,policy.id,policy.name,policy,
            "ibm.constraint.project");
    }
    if (policy.kind == ExecutionPolicyKind::BoundaryClosure) {
        return leaf(PlanNodeKind::Update,policy.id,policy.name,policy);
    }
    throw std::runtime_error("Unknown execution policy kind.");
}

SolvePlanNode compileExplicitStep(
        const std::vector<ExecutionPolicy>& policies,
        std::string_view timeIntegrator,
        std::set<std::string>& consumed) {
    const auto op = [](std::string id, std::string name, OpId operation) {
        return operationLeaf(PlanNodeKind::Update,std::move(id),
                             std::move(name),{},std::move(operation));
    };
    SolvePlanNode root;
    root.kind = PlanNodeKind::Sequence;
    root.id = "Explicit.step";
    root.name = "explicit time step";
    root.children.push_back(op(
        "Explicit.prepare","prepare physical step","flow.step.prepare"));
    root.children.push_back(op(
        "Explicit.dt","compute stable time step","flow.dt.compute"));
    root.children.push_back(op(
        "Explicit.begin","begin explicit step","flow.step.begin"));

    ExecutionPolicy merged;
    merged.id = "Explicit.stages";
    merged.name = std::string(timeIntegrator)+" fused explicit stages";
    bool foundExplicit = false;
    for (const auto& policy : policies) {
        if (policy.kind == ExecutionPolicyKind::ExplicitStages) {
            foundExplicit = true;
            consumed.insert(policy.id);
            merged.equations.insert(merged.equations.end(),
                                    policy.equations.begin(),policy.equations.end());
            merged.constraints.insert(merged.constraints.end(),
                                      policy.constraints.begin(),policy.constraints.end());
            merged.unknowns.insert(merged.unknowns.end(),
                                   policy.unknowns.begin(),policy.unknowns.end());
        } else if (policy.kind == ExecutionPolicyKind::BoundaryClosure) {
            // Boundary/interface closures are consumed by the fused explicit
            // RHS stage. DensityBasedRHS preserves the validated
            // boundary/halo/IBM ordering.
            consumed.insert(policy.id);
        }
    }
    if (!foundExplicit) {
        throw std::runtime_error("Explicit plan has no ExplicitStages policy.");
    }
    SolvePlanNode stages = leaf(
        PlanNodeKind::StageLoop,"Explicit.stages",merged.name,merged);
    stages.repetitions = FDM::explicitStageCount(
        FDM::parseTimeScheme(std::string(timeIntegrator)));
    stages.children.push_back(operationLeaf(
        PlanNodeKind::Update,"Explicit.stage.execute","execute explicit stage",
        {},"explicit.stage.execute"));
    root.children.push_back(std::move(stages));

    for (const auto& policy : policies) {
        if (policy.kind == ExecutionPolicyKind::ConstraintProjection
            || policy.kind == ExecutionPolicyKind::MonolithicKKT) {
            root.children.push_back(compilePolicy(policy,timeIntegrator));
            consumed.insert(policy.id);
        }
    }
    root.children.push_back(op(
        "Explicit.commit","commit physical state","flow.step.commit"));
    root.children.push_back(op(
        "Explicit.time.commit","commit time and step","time.commit"));
    return root;
}

bool hasEquationPrefix(
        const ExecutableEquationSystem& system,
        std::string_view prefix) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& item) {
            return item.id.compare(0,prefix.size(),prefix) == 0;
        });
}

bool hasTurbulenceClosure(const ExecutableEquationSystem& system) {
    return std::any_of(
        system.closures.begin(),system.closures.end(),
        [](const std::string& closure) {
            return closure.find("turbulence") != std::string::npos
                || closure.find("mu_t") != std::string::npos;
        });
}

void requireAllPoliciesConsumed(
        const std::vector<ExecutionPolicy>& policies,
        const std::set<std::string>& consumed) {
    for (const auto& policy : policies) {
        if (consumed.find(policy.id) == consumed.end()) {
            throw std::runtime_error(
                "SolvePlanner cannot lower active execution policy '"
                +policy.id+"' together with the selected structured plan.");
        }
    }
}

} // namespace

CompiledSolvePlan SolvePlanner::compile(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies,
        std::string_view timeIntegrator) {
    CompiledSolvePlan result;
    result.root.kind = PlanNodeKind::Sequence;
    result.root.id = "PLAN_ROOT";
    result.root.name = "time-step sequence";

    auto ordered = policies;
    std::stable_sort(ordered.begin(),ordered.end(),
        [](const ExecutionPolicy& left, const ExecutionPolicy& right) {
            return left.priority < right.priority;
        });
    for (const auto& policy : ordered) {
        if (policy.strategyKind == FDM::SolveStrategyKind::AlgebraicUpdate) {
            throw std::runtime_error(
                "Execution policy '"+policy.id+"' has no typed strategy.");
        }
        result.blocks.push_back({
            policy.id,policy.name,policy.strategyName,policy.equations,
            policy.constraints,policy.unknowns,policy.strategyKind});
    }

    const auto eulerian = std::find_if(ordered.begin(),ordered.end(),[](const ExecutionPolicy& policy) {
        return policy.id == "S_EE_PIMPLE";
    });
    if (eulerian != ordered.end()) {
        std::set<std::string> consumed{"S_EE_PIMPLE"};
        const auto turbulence = std::find_if(
            ordered.begin(),ordered.end(),[](const ExecutionPolicy& policy) {
                return policy.id == "S_TURBULENCE";
            });
        if (turbulence != ordered.end()) consumed.insert(turbulence->id);
        const bool solveTurbulence = hasEquationPrefix(system,"E_TURB_");
        const bool prepareTurbulence = solveTurbulence
            || hasTurbulenceClosure(system);
        result.root = compileEulerianPimple(
            *eulerian,prepareTurbulence,solveTurbulence);
        requireAllPoliciesConsumed(ordered,consumed);
        return result;
    }

    const ExecutionCapabilitySignature signature = capabilitySignature(system,ordered);
    if (supportsGenericPiso(signature,ordered,timeIntegrator)) {
        const auto pressure = std::find_if(
            ordered.begin(),ordered.end(),[](const ExecutionPolicy& item) {
                return item.id == "S_PRESSURE";
            });
        result.root = compileGenericPiso(*pressure);
        requireAllPoliciesConsumed(
            ordered,{"S_PREDICTOR","S_PRESSURE"});
        return result;
    }

    if (signature.constantDensity && signature.pressureConstraint) {
        return result;
    }

    const auto pressure = std::find_if(
        ordered.begin(),ordered.end(),
        [](const ExecutionPolicy& item) { return item.id == "S_PRESSURE"; });
    const auto predictor = std::find_if(
        ordered.begin(),ordered.end(),
        [](const ExecutionPolicy& item) { return item.id == "S_PREDICTOR"; });
    std::set<std::string> consumed;
    const bool hasExplicitStages = std::any_of(
        ordered.begin(),ordered.end(),
        [](const ExecutionPolicy& item) {
            return item.kind == ExecutionPolicyKind::ExplicitStages;
        });
    if (hasExplicitStages) {
        result.root = compileExplicitStep(ordered,timeIntegrator,consumed);
        requireAllPoliciesConsumed(ordered,consumed);
        return result;
    }
    for (const auto& policy : ordered) {
        if (pressure != ordered.end() && policy.id == "S_PREDICTOR") {
            consumed.insert(policy.id);
            continue;
        }
        if (policy.id == "S_PRESSURE" && predictor != ordered.end()) {
            ExecutionPolicy coupled = policy;
            coupled.equations.insert(
                coupled.equations.begin(),
                predictor->equations.begin(),predictor->equations.end());
            coupled.unknowns.insert(
                coupled.unknowns.begin(),
                predictor->unknowns.begin(),predictor->unknowns.end());
            result.root.children.push_back(
                compilePolicy(coupled,timeIntegrator));
        } else {
            result.root.children.push_back(
                compilePolicy(policy,timeIntegrator));
        }
        consumed.insert(policy.id);
    }

    requireAllPoliciesConsumed(ordered,consumed);

    return result;
}

ExecutionCapabilitySignature SolvePlanner::capabilities(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies) {
    return capabilitySignature(system,policies);
}

std::vector<OpId> SolvePlanner::requiredOperations(const CompiledSolvePlan& plan) {
    std::vector<OpId> result;
    const auto visit = [&](const auto& self, const SolvePlanNode& node) -> void {
        if (!node.operation.empty()
            && std::find(result.begin(),result.end(),node.operation) == result.end()) {
            result.push_back(node.operation);
        }
        for (const auto& child : node.children) self(self,child);
    };
    visit(visit,plan.root);
    return result;
}

} // namespace SF::System
