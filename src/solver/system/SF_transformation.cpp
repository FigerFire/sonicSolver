#include "core/system/SF_operationIds.h"
/// @file SF_transformation.cpp
/// @brief Equation system transformer registry、match、validation 与 application。

#include "SF_transformation.h"
#include "SF_eulerianRelations.h"
#include "SF_pressureCoupling.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::System {
namespace {

class PressureConstraintTransformer final : public IEquationSystemTransformer {
public:
    explicit PressureConstraintTransformer(PressureCouplingInputs inputs,bool fixedTime)
        : inputs_(std::move(inputs)), fixedTime_(fixedTime), descriptor_{"pressureConstraint","pressure constraint",100,false,
                      {OriginKind::BuiltinPreset,"pressure coupling"}} {}

    const TransformationDescriptor& descriptor() const override {
        return descriptor_;
    }

    TransformationMatch match(const RawEquationSystem& raw) const override {
        if (!hasConstraint(raw,inputs_.constraint)) {
            return {TransformationState::RegisteredButNotApplicable,
                    "required incompressibility constraint not present"};
        }
        if (!hasEquation(raw,inputs_.predictor.equation) || !hasUnknown(raw,"p")) {
            return {TransformationState::RegisteredButInvalid,
                    "requires Momentum equation and pressure multiplier p"};
        }
        return {TransformationState::RegisteredAndActive,
                "pressure constraint contract matched"};
    }

    void transform(
            const RawEquationSystem& raw,
            ExecutableEquationSystemBuilder& executable,
            TransformationRecord& record) const override {
        std::vector<std::string> priorEquations;
        for (const auto& equation:executable.system().registry.entries()) priorEquations.push_back(equation.id);
        const auto velocity = std::find_if(
            raw.state.symbols().begin(),raw.state.symbols().end(),
            [](const StateSymbol& unknown) { return unknown.id == "U"; });
        if (velocity == raw.state.symbols().end()) {
            throw std::runtime_error(
                "Pressure constraint requires a declared mathematical U unknown.");
        }
        const bool primitiveTarget=inputs_.representation==PressurePredictorKind::Velocity;
        if (inputs_.predictorTarget.symbol!=(primitiveTarget ? "U" : "rhoU"))
            throw std::runtime_error("Pressure transformation requires an explicit momentum solution target.");
        const auto relation=[&](const EquationCall& call,Equation equation) {
            if (call.equation.empty()) throw std::runtime_error("Pressure transformation requires an explicit relation reference.");
            equation.id=call.equation;executable.addEquation(std::move(equation));
        };
        if (!primitiveTarget) {
            const auto equations=conservativePressureRelations();
            const std::vector<EquationCall> references{inputs_.correction,inputs_.velocityUpdate,inputs_.pressureUpdate,
                inputs_.fluxUpdate,inputs_.publication,inputs_.relaxation,inputs_.fluxConsistency,inputs_.convergence};
            for (std::size_t i=0;i<equations.size();++i) {
                auto equation=equations[i];
                if (!fixedTime_ && (equation.id=="relaxIterate" || equation.id=="restoreFlux" || equation.id=="checkConvergence")) continue;
                relation(references[i],std::move(equation));
            }
        }
        if (primitiveTarget) {
            using Expr=FormulaExpr;
            const Provenance origin{OriginKind::Generated,"pressureConstraint"};
            // These formulas describe the operations implemented by the
            // current pressure provider. rAU, HbyA, and faceResponse are
            // derived numerical workspace, not physical state authorities.
            relation(inputs_.correction,{"pSimple",
                Expr::negate(Expr::op(Math::Op::Div,{
                    Expr::multiply(PressureMathematics::symbol("faceResponse"),
                        Expr::op(Math::Op::Grad,{PressureMathematics::symbol("pPrime")}))},
                    "pressureLaplacian")),
                Expr::negate(PressureMathematics::op("faceVolumeDivergence",{PressureMathematics::symbol("correctedFlux")},
                                     "continuityDefect")),origin});
            relation(inputs_.pressureUpdate,{"correctP",PressureMathematics::symbol("p"),
                Expr::add(PressureMathematics::symbol("p"),PressureMathematics::symbol("pPrime")),origin});
            relation(inputs_.velocityUpdate,{"correctU",PressureMathematics::symbol("U"),
                Expr::subtract(PressureMathematics::symbol("U"),
                    Expr::multiply(PressureMathematics::symbol("rAU"),
                        Expr::op(Math::Op::Grad,{PressureMathematics::symbol("pPrime")},
                                 "correction.pressureGradient"))),origin});
            relation(inputs_.fluxUpdate,{"correctFluxp",PressureMathematics::symbol("phi"),
                Expr::subtract(PressureMathematics::symbol("phi"),
                    PressureMathematics::op("pressureFlux",{
                        PressureMathematics::symbol("faceResponse"),PressureMathematics::symbol("pPrime")},
                        "correction.faceFlux")),origin});
        }
        if (primitiveTarget && fixedTime_) {
            using Expr=FormulaExpr;
            const Provenance origin{OriginKind::Generated,"pressure coupling relations"};
            relation(inputs_.relaxation,{"relaxIterate",PressureMathematics::symbol("iterate"),
                PressureMathematics::op("relax",{PressureMathematics::symbol("iterate"),PressureMathematics::symbol("laggedIterate")}),origin});
            relation(inputs_.fluxConsistency,{"restoreFlux",PressureMathematics::symbol("phi"),
                PressureMathematics::op("consistentFlux",{PressureMathematics::symbol("U"),PressureMathematics::symbol("p")}),origin});
            relation(inputs_.convergence,{"checkConvergence",PressureMathematics::symbol("converged"),
                PressureMathematics::op("residualConvergence",{PressureMathematics::symbol("iterate"),PressureMathematics::symbol("laggedIterate")}),origin});
        }
        executable.addOperator({
            "OP_PRESSURE_UPDATE","pressure update",{"pPrime"},{"p"},
            {OriginKind::Generated,"pressureConstraint"}});
        executable.addOperator({
            "OP_VELOCITY_CORRECTION","velocity correction",
            {"pPrime",primitiveTarget ? "U" : "rhoU"},{primitiveTarget ? "U" : "rhoU"},
            {OriginKind::Generated,"pressureConstraint"}});
        executable.addOperator({
            "OP_FLUX_CORRECTION","flux correction",{"pPrime"},
            {primitiveTarget ? "faceFlux" : "fluxValidity"},
            {OriginKind::Generated,"pressureConstraint"}});
        for (const auto& equation:executable.system().registry.entries())
            if (std::find(priorEquations.begin(),priorEquations.end(),equation.id)==priorEquations.end())
                record.generatedEquations.push_back(equation.id);
        record.generatedOperators.insert(record.generatedOperators.end(),{
            "OP_PRESSURE_UPDATE","OP_VELOCITY_CORRECTION","OP_FLUX_CORRECTION"});
    }

private:
    PressureCouplingInputs inputs_;
    bool fixedTime_;
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
            raw.registry.entries().begin(),raw.registry.entries().end(),
            [](const Equation& item) {
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
            TransformationRecord& record) const override {
        for (auto equation:eulerianPressureRelations()) {
            equation.origin={OriginKind::Generated,"sharedPressureConstraint"};
            executable.addEquation(equation);
            record.generatedEquations.push_back(equation.id);
        }

    }

private:
    TransformationDescriptor descriptor_;
};

bool hasUnknownId(
        const std::vector<StateSymbol>& unknowns,
        std::string_view id) {
    return std::any_of(unknowns.begin(),unknowns.end(),
        [&](const StateSymbol& item) { return item.id == id; });
}

} // namespace

ExecutableEquationSystemBuilder::ExecutableEquationSystemBuilder(
        const RawEquationSystem& raw) {
    system_.immersed = raw.immersed;
    system_.state = raw.state;
    system_.registry = raw.registry;
    system_.constraints = raw.constraints;
    system_.closures = raw.closures;
    system_.boundaries = raw.boundaries;
    system_.boundaryClosures = raw.boundaryClosures;
    system_.dependencies = raw.dependencies;
}

ExecutableEquationSystemBuilder::ExecutableEquationSystemBuilder(const ExecutableEquationSystem& composed)
    : system_(composed) {}

void ExecutableEquationSystemBuilder::addState(StateSymbol unknown) {
    if (hasUnknownId(system_.state.symbols(),unknown.id)) {
        throw std::runtime_error("Transformation generated duplicate unknown '"
                                 +unknown.id+"'.");
    }
    system_.state.add(std::move(unknown));
}

void ExecutableEquationSystemBuilder::replaceEquation(Equation formula) {
    system_.state.validateSymbols(formula.lhs);system_.state.validateSymbols(formula.rhs);
    system_.registry.replace(std::move(formula));
}

void ExecutableEquationSystemBuilder::addEquation(Equation formula) {
    system_.state.validateSymbols(formula.lhs);system_.state.validateSymbols(formula.rhs);
    system_.registry.add(std::move(formula));
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
        || !system_.registry.contains(equation.equationId)) {
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
            transformer->transform(raw,executable,record);
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

std::unique_ptr<IEquationSystemTransformer> makePressureConstraintTransformer(PressureCouplingInputs inputs,bool fixedTime) {
    return std::make_unique<PressureConstraintTransformer>(std::move(inputs),fixedTime);
}

std::unique_ptr<IEquationSystemTransformer> makeSharedPressureTransformer() {
    return std::make_unique<SharedPressureTransformer>();
}



} // namespace SF::System
