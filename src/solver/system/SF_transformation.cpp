/// @file SF_transformation.cpp
/// @brief Equation system transformer registry、match、validation 与 application。

#include "SF_transformation.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::System {
namespace {

class PressureConstraintTransformer final : public IEquationSystemTransformer {
public:
    PressureConstraintTransformer()
        : descriptor_{"pressureConstraint","pressure constraint",100,false,
                      {OriginKind::BuiltinPreset,"pressure coupling"}} {}

    const TransformationDescriptor& descriptor() const override {
        return descriptor_;
    }

    TransformationMatch match(const RawEquationSystem& raw) const override {
        if (!hasConstraint(raw,"C_INCOMPRESSIBILITY")) {
            return {TransformationState::RegisteredButNotApplicable,
                    "required incompressibility constraint not present"};
        }
        if (!hasEquation(raw,"E_MOMENTUM") || !hasUnknown(raw,"p")) {
            return {TransformationState::RegisteredButInvalid,
                    "requires Momentum equation and pressure multiplier p"};
        }
        return {TransformationState::RegisteredAndActive,
                "pressure constraint contract matched"};
    }

    void transform(
            const RawEquationSystem& raw,
            ExecutableEquationSystemBuilder& executable,
            std::vector<ExecutionPolicy>&,
            TransformationRecord& record) const override {
        UnknownDescriptor correction;
        correction.id = "pPrime";
        correction.name = "pressure correction";
        correction.role = UnknownRole::Algebraic;
        correction.storageBinding = StorageBinding::TransientWorkspace;
        correction.storageKey = "pressureCorrection";
        correction.runtimeStorageRequired = false;
        correction.boundaryRequired = false;
        correction.restartEligible = false;
        correction.outputEligible = false;
        correction.origin = {OriginKind::Generated,"pressureConstraint"};
        executable.addUnknown(std::move(correction));

        EquationDescriptor predictor{
            "E_MOMENTUM_PREDICTOR","momentum predictor","algorithmic",{"U"}};
        predictor.category = EquationCategory::AlgorithmicDerivedEquation;
        predictor.origin = {OriginKind::Generated,"pressureConstraint"};
        executable.addEquation(std::move(predictor),
            Equation::named("E_MOMENTUM_PREDICTOR",
                Equation::ddt({"U"}) + Equation::div({"momentumFlux"})
                    == Equation::Symbol{"zero"}));

        EquationDescriptor descriptor{
            "E_PRESSURE","pressure correction","constraint",{"pPrime","U"}};
        descriptor.category = EquationCategory::AlgorithmicDerivedEquation;
        descriptor.origin = {OriginKind::Generated,"pressureConstraint"};
        executable.addEquation(std::move(descriptor),
            Equation::named("E_PRESSURE",
                Equation::constraint({"pressureVelocityConsistency"})
                    == Equation::Symbol{"zero"}));
        executable.addOperator({
            "OP_PRESSURE_UPDATE","pressure update",{"pPrime"},{"p"},
            {OriginKind::Generated,"pressureConstraint"}});
        executable.addOperator({
            "OP_VELOCITY_CORRECTION","velocity correction",{"pPrime","U"},{"U"},
            {OriginKind::Generated,"pressureConstraint"}});
        executable.addOperator({
            "OP_FLUX_CORRECTION","flux correction",{"pPrime"},{"faceFlux"},
            {OriginKind::Generated,"pressureConstraint"}});
        // Pressure formulation 是"存在哪些 derived operation"的唯一 authority：
        // executable operation 的 id/stage/capability 只在这里声明一次。
        // Ordering 不在这里——plan fragment 只引用 stage。
        const Provenance source{OriginKind::Generated,"pressureConstraint"};
        const auto declare = [&](OperationStage stage, const char* id,
                                 const char* name,
                                 std::vector<OperationCapability> needs) {
            executable.addExecutableOperation(
                {id,name,stage,std::move(needs),source});
        };
        declare(OperationStage::Prepare,"pressure.prepare",
                "prepare pressure schedule",{OperationCapability::PressureSchedule});
        declare(OperationStage::MomentumAssemble,"momentum.assemble",
                "assemble momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::MomentumSolve,"momentum.solve",
                "advance momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::PressureBoundaryPrepare,
                "pressure.boundary.prepare","prepare pressure boundary state",
                {OperationCapability::PressureBoundary});
        declare(OperationStage::PressureAssemble,"pressure.assemble",
                "assemble pressure correction",{OperationCapability::PressureCorrection,
                                                  OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureSolve,"pressure.solve",
                "solve pressure correction",{OperationCapability::PressureCorrection,
                                               OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureUpdatePrepare,
                "pressure.update.prepare","prepare pressure update",
                {OperationCapability::PressureCorrection});
        declare(OperationStage::VelocityCorrect,"velocity.correct",
                "velocity correction",{OperationCapability::VelocityCorrection});
        declare(OperationStage::FluxCorrect,"flux.correct",
                "refresh derived face-flux state",{OperationCapability::FluxCorrection});
        declare(OperationStage::CorrectionCommit,
                "pressure.correction.commit",
                "commit pressure-corrected state",{OperationCapability::PressureCorrection});
        declare(OperationStage::StepCommit,"pressure.step.commit",
                "commit corrected state",{OperationCapability::PressureSchedule});
        executable.addCompiledEquation({
            "E_MOMENTUM_PREDICTOR",
            "momentum.predictor",
            {
                {"Q","conservative",0,5,ResourceAccessMode::ReadWrite,true,
                 SynchronizationRequirement::ReadHalo},
                {"R_momentum","residual",1,3,ResourceAccessMode::ReadWrite,false,
                 SynchronizationRequirement::WriteOwned}
            },
            false,true,{OriginKind::Generated,"pressureConstraint"}});
        executable.addCompiledEquation({
            "E_PRESSURE",
            "pressure.correction",
            {
                {"Q","conservative",0,5,ResourceAccessMode::ReadWrite,true,
                 SynchronizationRequirement::ReadHalo},
                {"pPrime","pressureCorrection",0,1,
                 ResourceAccessMode::ReadWrite,true,
                 SynchronizationRequirement::ReadHalo},
                {"A_p","pressureMatrix",0,1,ResourceAccessMode::Write,false,
                 SynchronizationRequirement::None},
                {"b_p","pressureRhs",0,1,ResourceAccessMode::Write,false,
                 SynchronizationRequirement::None}
            },
            true,true,{OriginKind::Generated,"pressureConstraint"}});
        record.generatedEquations.push_back("E_MOMENTUM_PREDICTOR");
        record.generatedEquations.push_back("E_PRESSURE");
        record.generatedOperators.insert(record.generatedOperators.end(),{
            "OP_PRESSURE_UPDATE","OP_VELOCITY_CORRECTION","OP_FLUX_CORRECTION"});
        for (const ExecutableOperation& operation : executable.operations()) {
            if (operation.origin.source == "pressureConstraint") {
                record.generatedOperations.push_back(operation.operation);
            }
        }
    }

private:
    TransformationDescriptor descriptor_;
};

class SharedPressureTransformer final : public IEquationSystemTransformer {
public:
    SharedPressureTransformer()
        : descriptor_{"sharedPressureConstraint","shared pressure constraint",
                      110,false,{OriginKind::BuiltinPreset,"eulerianEulerian"}} {}

    const TransformationDescriptor& descriptor() const override {
        return descriptor_;
    }

    TransformationMatch match(const RawEquationSystem& raw) const override {
        if (!hasConstraint(raw,"C_SHARED_PRESSURE")) {
            return {TransformationState::RegisteredButNotApplicable,
                    "shared-pressure constraint not present"};
        }
        const bool hasPhaseContinuity = std::any_of(
            raw.equations.begin(),raw.equations.end(),
            [](const EquationDescriptor& item) {
                return item.id.rfind("E_CONTINUITY.",0) == 0;
            });
        if (!hasUnknown(raw,"p") || !hasPhaseContinuity) {
            return {TransformationState::RegisteredButInvalid,
                    "requires shared pressure and per-phase continuity equations"};
        }
        return {TransformationState::RegisteredAndActive,
                "shared-pressure contract matched"};
    }

    void transform(
            const RawEquationSystem& raw,
            ExecutableEquationSystemBuilder& executable,
            std::vector<ExecutionPolicy>&,
            TransformationRecord& record) const override {
        EquationDescriptor descriptor{
            "E_SHARED_PRESSURE","shared pressure correction","constraint",{"p"}};
        descriptor.category = EquationCategory::AlgorithmicDerivedEquation;
        descriptor.origin = {OriginKind::Generated,"sharedPressureConstraint"};
        executable.addEquation(std::move(descriptor),
            Equation::named("E_SHARED_PRESSURE",
                Equation::constraint({"sharedPressure"})
                    == Equation::Symbol{"zero"}));
        executable.addOperator({
            "OP_PHASE_FLUX_CORRECTION","phase flux correction",{"p"},
            {"phaseFaceFlux"},{OriginKind::Generated,"sharedPressureConstraint"}});
        record.generatedEquations.push_back("E_SHARED_PRESSURE");
        record.generatedOperators.push_back("OP_PHASE_FLUX_CORRECTION");
        const Provenance source{OriginKind::Generated,"sharedPressureConstraint"};
        const auto declare = [&](const char* id, bool linear = false) {
            std::vector<OperationCapability> needs{
                OperationCapability::EulerianPhaseExecution};
            if (linear) needs.push_back(OperationCapability::PressureLinearSolve);
            executable.addExecutableOperation(
                {id,id,OperationStage::Prepare,std::move(needs),source});
            record.generatedOperations.push_back(id);
        };
        for (const char* id : {
                "ee.dt.compute", "ee.step.begin", "ee.interphase.compute",
                "ee.sources.assemble", "ee.sources.validate",
                "ee.momentum.diagonal", "ee.momentum.flux",
                "ee.faceFlux.canonical", "ee.continuity.assemble",
                "ee.boundary.prepare", "ee.momentum.solve",
                "ee.interphase.correct", "ee.boundary.afterMomentum",
                "ee.diagonal.sync", "ee.momentum.flux.after",
                "ee.faceFlux.canonical.after", "ee.pressure.publish",
                "ee.pressure.sync", "ee.phase.correct",
                "ee.faceFlux.correct", "ee.boundary.afterPressure",
                "ee.faceFlux.canonical.pressure", "ee.energy.solve",
                "ee.boundary.final", "ee.outer.validate",
                "ee.step.commit", "ee.time.commit"}) {
            declare(id);
        }
        declare("ee.pressure.solve",true);
        const bool solveTurbulence = std::any_of(
            raw.equations.begin(),raw.equations.end(),
            [](const EquationDescriptor& item) {
                return item.id.rfind("E_TURB_",0) == 0;
            });
        const bool prepareTurbulence = solveTurbulence || std::any_of(
            raw.closures.begin(),raw.closures.end(),
            [](const std::string& closure) {
                return closure.find("turbulence") != std::string::npos
                    || closure.find("mu_t") != std::string::npos;
            });
        if (prepareTurbulence) declare("ee.turbulence.prepare");
        if (solveTurbulence) declare("ee.turbulence.solve");
    }

private:
    TransformationDescriptor descriptor_;
};

class ImmersedConstraintTransformer final : public IEquationSystemTransformer {
public:
    ImmersedConstraintTransformer()
        : descriptor_{"immersedConstraint","immersed constraint",200,false,
                      {OriginKind::BuiltinPreset,"immersed boundary"}} {}

    const TransformationDescriptor& descriptor() const override {
        return descriptor_;
    }

    TransformationMatch match(const RawEquationSystem& raw) const override {
        const bool hasImmersedConstraint = std::any_of(
            raw.constraints.begin(),raw.constraints.end(),
            [](const ConstraintDescriptor& item) {
                return item.id.rfind("C_IBM",0) == 0;
            });
        if (!hasImmersedConstraint) {
            return {TransformationState::RegisteredButNotApplicable,
                    "variational IBM constraint not present"};
        }
        return {TransformationState::RegisteredAndActive,
                "immersed constraint contract matched"};
    }

    void transform(
            const RawEquationSystem&,
            ExecutableEquationSystemBuilder&,
            std::vector<ExecutionPolicy>&,
            TransformationRecord&) const override {
        // Current descriptor already supplies constraint equations. The registry
        // owns the transformation identity; numerical lowering remains in the
        // explicitly reported legacy IBM adapter.
    }

private:
    TransformationDescriptor descriptor_;
};

bool hasEquationId(
        const std::vector<EquationDescriptor>& equations,
        std::string_view id) {
    return std::any_of(equations.begin(),equations.end(),
        [&](const EquationDescriptor& item) { return item.id == id; });
}

bool hasUnknownId(
        const std::vector<UnknownDescriptor>& unknowns,
        std::string_view id) {
    return std::any_of(unknowns.begin(),unknowns.end(),
        [&](const UnknownDescriptor& item) { return item.id == id; });
}

} // namespace

ExecutableEquationSystemBuilder::ExecutableEquationSystemBuilder(
        const RawEquationSystem& raw) {
    system_.unknowns = raw.unknowns;
    system_.equations = raw.equations;
    system_.equationDefinitions = raw.equationDefinitions;
    system_.constraints = raw.constraints;
    system_.closures = raw.closures;
    system_.boundaries = raw.boundaries;
    system_.dependencies = raw.dependencies;
}

void ExecutableEquationSystemBuilder::addUnknown(UnknownDescriptor unknown) {
    if (hasUnknownId(system_.unknowns,unknown.id)) {
        throw std::runtime_error("Transformation generated duplicate unknown '"
                                 +unknown.id+"'.");
    }
    system_.unknowns.push_back(std::move(unknown));
}

void ExecutableEquationSystemBuilder::addEquation(
        EquationDescriptor descriptor, Equation::Definition definition) {
    if (descriptor.id != definition.name || hasEquationId(system_.equations,descriptor.id)) {
        throw std::runtime_error(
            "Transformation generated invalid/duplicate equation '"
            +descriptor.id+"'.");
    }
    system_.equations.push_back(std::move(descriptor));
    system_.equationDefinitions.add(std::move(definition));
}

void ExecutableEquationSystemBuilder::addOperator(
        GeneratedOperatorDescriptor operation) {
    if (operation.id.empty()) {
        throw std::runtime_error("Generated correction operator requires an id.");
    }
    const auto duplicate = std::find_if(
        system_.correctionOperators.begin(),system_.correctionOperators.end(),
        [&](const GeneratedOperatorDescriptor& item) {
            return item.id == operation.id;
        });
    if (duplicate != system_.correctionOperators.end()) {
        throw std::runtime_error(
            "Transformation generated duplicate operator '"+operation.id+"'.");
    }
    system_.correctionOperators.push_back(std::move(operation));
}

void ExecutableEquationSystemBuilder::addCompiledEquation(
        CompiledEquation equation) {
    if (equation.equationId.empty()
        || !hasEquationId(system_.equations,equation.equationId)) {
        throw std::runtime_error(
            "Compiled equation requires an executable equation definition.");
    }
    const auto duplicate = std::find_if(
        system_.compiledEquations.begin(),system_.compiledEquations.end(),
        [&](const CompiledEquation& item) {
            return item.equationId == equation.equationId;
        });
    if (duplicate != system_.compiledEquations.end()) {
        throw std::runtime_error(
            "Duplicate compiled equation '"+equation.equationId+"'.");
    }
    system_.compiledEquations.push_back(std::move(equation));
}

void ExecutableEquationSystemBuilder::addExecutableOperation(
        ExecutableOperation operation) {
    if (operation.operation.empty()) {
        throw std::runtime_error(
            "Executable operation requires a runtime operation id.");
    }
    const auto duplicate = std::find_if(
        system_.operations.begin(),system_.operations.end(),
        [&](const ExecutableOperation& item) {
            return item.operation == operation.operation;
        });
    if (duplicate != system_.operations.end()) {
        throw std::runtime_error(
            "Transformation declared duplicate executable operation '"
            +operation.operation+"'.");
    }
    system_.operations.push_back(std::move(operation));
}

void TransformerRegistry::registerTransformer(
        std::unique_ptr<IEquationSystemTransformer> transformer) {
    if (!transformer) throw std::runtime_error("Null equation-system transformer.");
    if (find(transformer->descriptor().id)) {
        throw std::runtime_error(
            "Duplicate equation-system transformer '"
            +transformer->descriptor().id+"'.");
    }
    transformers_.push_back(std::move(transformer));
}

const IEquationSystemTransformer* TransformerRegistry::find(
        std::string_view id) const {
    const auto found = std::find_if(
        transformers_.begin(),transformers_.end(),
        [&](const auto& item) { return item->descriptor().id == id; });
    return found == transformers_.end() ? nullptr : found->get();
}

ExecutableEquationSystem TransformationPipeline::apply(
        const RawEquationSystem& raw,
        const std::vector<TransformationDescriptor>& requests,
        const TransformerRegistry& registry,
        std::vector<ExecutionPolicy>& policies,
        std::vector<TransformationRecord>& records) {
    ExecutableEquationSystemBuilder executable(raw);
    auto ordered = requests;
    std::stable_sort(ordered.begin(),ordered.end(),
        [](const auto& left, const auto& right) {
            return left.priority < right.priority;
        });
    for (const auto& request : ordered) {
        TransformationRecord record;
        record.descriptor = request;
        const auto* transformer = registry.find(request.id);
        if (!transformer) {
            record.state = TransformationState::NotRegistered;
            record.reason = "no transformer implementation is registered";
            records.push_back(record);
            throw std::runtime_error(
                "Requested equation-system transformer '"+request.id
                +"' is not registered.");
        }
        const TransformationMatch match = transformer->match(raw);
        record.state = match.state;
        record.reason = match.reason;
        if (match.state == TransformationState::RegisteredAndActive) {
            transformer->transform(raw,executable,policies,record);
            record.state = TransformationState::Applied;
        } else if (request.explicitlyRequested
                   || match.state == TransformationState::RegisteredButInvalid) {
            records.push_back(record);
            throw std::runtime_error(
                "Equation-system transformer '"+request.id+"' cannot apply: "
                +match.reason+".");
        }
        records.push_back(std::move(record));
    }
    return executable.finish();
}

std::unique_ptr<IEquationSystemTransformer> makePressureConstraintTransformer() {
    return std::make_unique<PressureConstraintTransformer>();
}

std::unique_ptr<IEquationSystemTransformer> makeSharedPressureTransformer() {
    return std::make_unique<SharedPressureTransformer>();
}

std::unique_ptr<IEquationSystemTransformer> makeImmersedConstraintTransformer() {
    return std::make_unique<ImmersedConstraintTransformer>();
}

} // namespace SF::System
