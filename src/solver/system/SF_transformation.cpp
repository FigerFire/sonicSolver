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
    explicit PressureConstraintTransformer(std::string momentumTarget)
        : momentumTarget_(std::move(momentumTarget)), descriptor_{"pressureConstraint","pressure constraint",100,false,
                      {OriginKind::BuiltinPreset,"pressure coupling"}} {}

    const TransformationDescriptor& descriptor() const override {
        return descriptor_;
    }

    TransformationMatch match(const RawEquationSystem& raw) const override {
        if (!hasConstraint(raw,"C_INCOMPRESSIBILITY")) {
            return {TransformationState::RegisteredButNotApplicable,
                    "required incompressibility constraint not present"};
        }
        if (!hasEquation(raw,"momentum") || !hasUnknown(raw,"p")) {
            return {TransformationState::RegisteredButInvalid,
                    "requires Momentum equation and pressure multiplier p"};
        }
        return {TransformationState::RegisteredAndActive,
                "pressure constraint contract matched"};
    }

    void transform(
            const RawEquationSystem& raw,
            ExecutableEquationSystemBuilder& executable,
            std::vector<LegacyExecutionPolicy>&,
            TransformationRecord& record) const override {
        const auto velocity = std::find_if(
            raw.state.symbols().begin(),raw.state.symbols().end(),
            [](const StateSymbol& unknown) { return unknown.id == "U"; });
        if (velocity == raw.state.symbols().end()) {
            throw std::runtime_error(
                "Pressure constraint requires a declared mathematical U unknown.");
        }
        const bool primitiveTarget=momentumTarget_=="U";
        if (!primitiveTarget && momentumTarget_!="rhoU")
            throw std::runtime_error("Pressure transformation requires an explicit momentum solution target.");
        if (!primitiveTarget) {
            for (auto equation:conservativePressureRelations()) {
                record.generatedEquations.push_back(equation.id);
                executable.addEquation(std::move(equation));
            }
        }
        if (primitiveTarget) {
            using Expr=FormulaExpr;
            const Provenance origin{OriginKind::Generated,"pressureConstraint"};
            // These formulas describe the operations implemented by the
            // current pressure provider. rAU, HbyA, and faceResponse are
            // derived numerical workspace, not physical state authorities.
            executable.addEquation({"pSimple",
                Expr::negate(Expr::op("div",{
                    Expr::multiply(Expr::symbol("faceResponse"),
                        Expr::op("grad",{Expr::symbol("pPrime")}))},
                    "pressureLaplacian")),
                Expr::negate(Expr::op("div",{Expr::symbol("correctedFlux")},
                                     "continuityDefect")),origin});
            executable.addEquation({"correctP",Expr::symbol("p"),
                Expr::add(Expr::symbol("p"),Expr::symbol("pPrime")),origin});
            executable.addEquation({"correctU",Expr::symbol("U"),
                Expr::subtract(Expr::symbol("U"),
                    Expr::multiply(Expr::symbol("rAU"),
                        Expr::op("grad",{Expr::symbol("pPrime")},
                                 "correction.pressureGradient"))),origin});
            executable.addEquation({"correctFluxp",Expr::symbol("phi"),
                Expr::subtract(Expr::symbol("phi"),
                    Expr::op("pressureFlux",{
                        Expr::symbol("faceResponse"),Expr::symbol("pPrime")},
                        "correction.faceFlux")),origin});
        }
        if (primitiveTarget) {
            using Expr=FormulaExpr;
            const Provenance origin{OriginKind::Generated,"pressure coupling relations"};
            executable.addEquation({"relaxIterate",Expr::symbol("iterate"),
                Expr::op("relax",{Expr::symbol("iterate"),Expr::symbol("laggedIterate")}),origin});
            executable.addEquation({"restoreFlux",Expr::symbol("phi"),
                Expr::op("consistentFlux",{Expr::symbol("U"),Expr::symbol("p")}),origin});
            executable.addEquation({"checkConvergence",Expr::symbol("converged"),
                Expr::op("residualConvergence",{Expr::symbol("iterate"),Expr::symbol("laggedIterate")}),origin});
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
        if (primitiveTarget) record.generatedEquations.push_back("pSimple");
        record.generatedOperators.insert(record.generatedOperators.end(),{
            "OP_PRESSURE_UPDATE","OP_VELOCITY_CORRECTION","OP_FLUX_CORRECTION"});
    }

private:
    std::string momentumTarget_;
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
            raw.legacyEquations.begin(),raw.legacyEquations.end(),
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
            std::vector<LegacyExecutionPolicy>&,
            TransformationRecord& record) const override {
        for (auto equation:eulerianPressureRelations()) {
            equation.origin={OriginKind::Generated,"sharedPressureConstraint"};
            equation.authored=true;
            if (equation.id=="E_SHARED_PRESSURE") {
                EquationDescriptor descriptor{equation.id,"shared volume-pressure relation","constraint",{"p"}};
                descriptor.category=EquationCategory::AlgorithmicDerivedEquation;
                executable.addEquation(std::move(descriptor),eulerianBackendDefinition(equation));
                executable.replaceEquation(equation);
            } else executable.addEquation(equation);
            record.generatedEquations.push_back(equation.id);
        }

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
            std::vector<LegacyExecutionPolicy>&,
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
        const std::vector<StateSymbol>& unknowns,
        std::string_view id) {
    return std::any_of(unknowns.begin(),unknowns.end(),
        [&](const StateSymbol& item) { return item.id == id; });
}

} // namespace

ExecutableEquationSystemBuilder::ExecutableEquationSystemBuilder(
        const RawEquationSystem& raw) {
    system_.state = raw.state;
    system_.legacyEquations = raw.legacyEquations;
    system_.legacyDefinitions = raw.legacyDefinitions;
    system_.registry = raw.registry;
    system_.constraints = raw.constraints;
    system_.closures = raw.closures;
    system_.boundaries = raw.boundaries;
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

void ExecutableEquationSystemBuilder::addEquation(
        EquationDescriptor descriptor, SF::Equation::Definition definition) {
    if (descriptor.id != definition.name || hasEquationId(system_.legacyEquations,descriptor.id)) {
        throw std::runtime_error(
            "Transformation generated invalid/duplicate equation '"
            +descriptor.id+"'.");
    }
    system_.registry.add(formulaFromEquation(definition,descriptor.origin));
    system_.legacyEquations.push_back(std::move(descriptor));
    system_.legacyDefinitions.add(std::move(definition));
}

void ExecutableEquationSystemBuilder::replaceEquation(Equation formula) {
    system_.registry.replace(std::move(formula));
}

void ExecutableEquationSystemBuilder::addEquation(Equation formula) {
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
        std::vector<LegacyExecutionPolicy>& policies,
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

std::unique_ptr<IEquationSystemTransformer> makePressureConstraintTransformer(std::string momentumTarget) {
    return std::make_unique<PressureConstraintTransformer>(std::move(momentumTarget));
}

std::unique_ptr<IEquationSystemTransformer> makeSharedPressureTransformer() {
    return std::make_unique<SharedPressureTransformer>();
}

std::unique_ptr<IEquationSystemTransformer> makeImmersedConstraintTransformer() {
    return std::make_unique<ImmersedConstraintTransformer>();
}

} // namespace SF::System
