#include "SF_eulerianCoupling.h"
#include "SF_eulerianRelations.h"
#include "SF_eulerianAssembly.h"
#include "SF_eulerianTurbulence.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
namespace SF::System {
namespace {
using E=FormulaExpr;
SolvePlanNode operation(const char* id,const std::vector<CompiledMathRef>& calls={}) {
    SolvePlanNode node;node.kind=PlanNodeKind::Update;node.id=std::string("EE.")+std::string(id).substr(3);
    node.name=node.id;node.operation=id;node.provider="flow.eulerian-pressure";node.equationCalls=calls;return node;
}
// These extensions feed the existing source ledger; arbitrary AST changes cannot inherit its kernel.
bool sourceExtensions(const E& expression,const E& base,const std::string& method,
        std::vector<EulerianSourceExtension>& extensions) {
    Equation a{"check",expression,E::constantValue(0),{}},b{"check",base,E::constantValue(0),{}};
    if (canonicalFormula(a)==canonicalFormula(b)) return true;
    if (expression.kind()!=E::Kind::Add || expression.arguments().size()!=2) return false;
    const auto& term=expression.arguments()[1];
    if (term.kind()!=E::Kind::Operator || term.operation()!=Math::Op::Source || term.arguments().size()!=1) return false;
    const auto& symbol=term.arguments()[0];
    if (symbol.kind()!=E::Kind::Symbol || !(method=="EulerianPhaseMomentum" && (symbol.name()=="gravity" || symbol.name()=="MRF")
        || method=="EulerianPhaseEnthalpy" && symbol.name()=="wallHeat")) return false;
    const auto extension=symbol.name()=="gravity" ? EulerianSourceExtension::Gravity
        : symbol.name()=="MRF" ? EulerianSourceExtension::MRF : EulerianSourceExtension::WallHeat;
    if (std::find(extensions.begin(),extensions.end(),extension)!=extensions.end()) return false;
    extensions.push_back(extension);
    return sourceExtensions(expression.arguments()[0],base,method,extensions);
}
class EulerianMethod final:public IProvider {
public:
    explicit EulerianMethod(std::string id):id_(std::move(id)) {}
    std::string_view id() const override { return id_; }
    std::string_view runtimeProvider() const override { return "flow.eulerian-pressure"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& step,const NumericalBinding& binding) const override {
        if (!hasConstraint(system,"C_SHARED_PRESSURE") || !hasConstraint(system,"C_VOLUME_FRACTION"))
            throw std::runtime_error("Eulerian provider requires shared pressure and volume closure.");
        const bool continuity=id_=="EulerianPhaseContinuity",momentum=id_=="EulerianPhaseMomentum",energy=id_=="EulerianPhaseEnthalpy";
        const bool pressure=id_=="EulerianSharedPressureCorrection",correction=id_=="EulerianPhaseCorrection",flux=id_=="EulerianPhaseFluxCorrection";
        CompiledEquationCall result;result.source=step;result.target.symbol=step.target.symbol;result.target.kind=step.target.kind;
        result.target.viewOwner=StateViewOwner::NumericalProvider;result.equationMethod=id_;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.calls={{step.equation,targetText(step.target)}};
        std::vector<const StateSymbol*> phases;
        for (const auto& symbol:system.state.symbols())
            if (symbol.id.rfind("phaseMass.",0)==0) phases.push_back(&symbol);
        const auto slot=[](const StateSymbol* phase) {
            const auto end=phase->storageKey.find('.');
            if (phase->storageKey.rfind("phase",0)!=0 || end==std::string::npos)
                throw std::runtime_error("Invalid Eulerian phase storage contract: "+phase->id);
            return std::stoi(phase->storageKey.substr(5,end-5));
        };
        std::sort(phases.begin(),phases.end(),[&](const auto* a,const auto* b) { return slot(a)<slot(b); });
        for (std::size_t i=0;i<phases.size();++i) {
            const auto prefix="phase"+std::to_string(i)+".";
            const auto suffix=phases[i]->id.substr(std::string("phaseMass").size());
            for (const auto& layout:std::vector<std::pair<std::string,std::string>>{
                    {phases[i]->id,"mass"},{"momentum"+suffix,"momentum"},{"enthalpy"+suffix,"enthalpy"}}) {
                const auto& state=system.state.at(layout.first);
                if (slot(phases[i])!=static_cast<int>(i) || state.storageKey!=prefix+layout.second
                    || state.storageBinding!=StorageBinding::ProviderDistributed || state.componentOffset!=0
                    || state.components!=(layout.second=="momentum" ? 3 : 1) || state.role!=StateRole::Transported)
                    throw std::runtime_error("Eulerian phase backing differs from the declared PhaseSystem order.");
            }
        }
        if (phases.size()<2) throw std::runtime_error("Eulerian provider requires at least two phase states.");
        std::vector<CompiledMathRef> required;
        for (const auto* phase:phases) {
            const auto suffix=phase->id.substr(std::string("phaseMass").size());
            const auto symbol=continuity ? phase->id : momentum || correction ? "momentum"+suffix
                : energy ? "enthalpy"+suffix : "volumeFaceFlux"+suffix;
            if (!pressure) required.push_back({continuity ? "E_CONTINUITY"+suffix : momentum ? "momentum"+suffix
                : energy ? "E_ENTHALPY"+suffix : correction ? "correctPhase"+suffix : "correctPhaseFlux"+suffix,symbol});
        }
        if (correction) required.push_back({"correctSharedP","p"});
        if (!pressure) {
            if (binding.inputs.size()!=required.size()) throw std::runtime_error("Eulerian fusion requires every phase target.");
            for (std::size_t i=0;i<required.size();++i)
                if (binding.inputs[i]!=required[i].target) throw std::runtime_error("Eulerian fusion phase order/target mismatch.");
            if (std::none_of(required.begin(),required.end(),[&](const auto& ref) { return ref.equation==step.equation && ref.target==step.target.symbol; }))
                throw std::runtime_error("Eulerian equation target is outside its declared group.");
            result.fusionKey=id_;result.fusionMembers=required;
        }
        const auto& shared=system.state.at("p");
        if (shared.storageKey!="pressure" || shared.components!=1 || shared.storageBinding!=StorageBinding::ProviderDistributed
            || shared.role!=StateRole::Algebraic || shared.readOnly())
            throw std::runtime_error("Eulerian shared pressure must alias the unique PhaseSystem pressure.");
        const auto& actual=system.registry.at(step.equation);
        CompiledEulerianAssemblyContract contract;
        contract.equation=step.equation;contract.target=step.target;
        contract.method=id_;contract.provider=std::string(runtimeProvider());
        Equation expected;
        if (pressure || step.equation=="correctSharedP") {
            expected=eulerianPressureRelations()[pressure ? 0 : 1];
        } else {
            const auto member=std::find_if(required.begin(),required.end(),[&](const auto& item) {
                return item.equation==step.equation;
            });
            const auto index=static_cast<std::size_t>(member-required.begin());
            contract.phaseSlot=index;contract.phaseName=phases.at(index)->id.substr(std::string("phaseMass.").size());
            const auto& suffix=contract.phaseName;
            const auto relations=eulerianPhaseRelations(suffix,dependsOn(actual.rhs,"referencePhaseRemainder"));
            expected=relations[continuity ? 0 : momentum ? 1 : energy ? 2 : correction ? 3 : 4];
        }
        Equation lhsA{actual.id,actual.lhs,E::constantValue(0),{}},lhsB{expected.id,expected.lhs,E::constantValue(0),{}};
        if (canonicalFormula(lhsA)!=canonicalFormula(lhsB) || !sourceExtensions(actual.rhs,expected.rhs,id_,contract.extensions))
            throw std::runtime_error("Unsupported Eulerian WHAT/kernel mathematical contract: "+step.equation);
        std::reverse(contract.extensions.begin(),contract.extensions.end());
        using R=EulerianAssemblyRelation;
        contract.relation=continuity ? (dependsOn(actual.rhs,"referencePhaseRemainder") ? R::ReferencePhaseClosure : R::PhaseMassBalance)
            : momentum ? R::PhaseMomentumBalance : energy ? R::PhaseEnthalpyBalance
            : pressure ? R::SharedPressureCorrection : flux ? R::PairedFaceFluxCorrection
            : step.equation=="correctSharedP" ? R::SharedPressureUpdate : R::PhaseMomentumCorrection;
        if (continuity) {
            const auto references=std::count_if(phases.begin(),phases.end(),[&](const auto* phase) {
                return dependsOn(system.registry.at("E_CONTINUITY."+phase->id.substr(std::string("phaseMass.").size())).rhs,"referencePhaseRemainder");
            });
            if (references!=1) throw std::runtime_error("Eulerian volume closure requires exactly one reference-phase contract.");
        }
        if (pressure) {
            if (step.equation!="E_SHARED_PRESSURE" || step.target.symbol!="deltaP" || step.target.kind!=TargetKind::Correction)
                throw std::runtime_error("Shared pressure requires registered deltaP workspace.");
            result.target.workspace="pressureCorrection";
        } else if (flux) {
            if (step.target.kind!=TargetKind::Workspace) throw std::runtime_error("Phase face flux is numerical workspace.");
            const auto suffix=step.target.symbol.substr(std::string("volumeFaceFlux").size());
            result.target.workspace=system.state.at("phaseMass"+suffix).storageKey;
            result.target.workspace.replace(result.target.workspace.size()-4,4,"volumeFaceFlux");
            result.target.resources={{step.target.symbol,result.target.workspace,0,3,ResourceAccessMode::Write}};
            result.target.workspaceLocation=VariableLocation::EulerianFace;
            result.target.workspaceOwnership=OwnershipKind::CanonicalFace;
        } else {
            const auto& state=system.state.at(step.target.symbol);
            if (step.target.kind!=TargetKind::Physical || state.storageBinding!=StorageBinding::ProviderDistributed
                || state.role==StateRole::Derived || state.components!=(momentum || (correction && step.target.symbol!="p") ? 3 : 1))
                throw std::runtime_error("Eulerian target must alias the existing writable provider backing.");
        }
        std::vector<const char*> ids;
        if (continuity) {
            result.stepInputs={{OpIds::EeDtCompute,std::string(runtimeProvider())},{OpIds::EeStepBegin,std::string(runtimeProvider())}};
            result.iterationInputs={{OpIds::EeInterphaseCompute,std::string(runtimeProvider())},{OpIds::EeSourcesAssemble,std::string(runtimeProvider())}};
            result.physicalPublication={OpIds::EeStepCommit,std::string(runtimeProvider())};
            ids.insert(ids.end(),{OpIds::EeSourcesValidate,OpIds::EeMomentumDiagonal,OpIds::EeMomentumFlux,
                OpIds::EeFaceFluxCanonical,OpIds::EeContinuityAssemble,OpIds::EeBoundaryPrepare});
        } else if (momentum) ids={OpIds::EeMomentumSolve,OpIds::EeInterphaseCorrect,OpIds::EeBoundaryAfterMomentum,
            OpIds::EeDiagonalSync,OpIds::EeMomentumFluxAfter,OpIds::EeFaceFluxCanonicalAfter};
        else if (pressure) ids={OpIds::EePressureSolve,OpIds::EePressurePublish,OpIds::EePressureSync};
        else if (correction) ids={OpIds::EePhaseCorrect};
        else if (flux) ids={OpIds::EeFaceFluxCorrect,OpIds::EeBoundaryAfterPressure,OpIds::EeFaceFluxCanonicalPressure};
        else {
            ids={OpIds::EeEnergySolve};
            result.iterationOutputs={{OpIds::EeBoundaryFinal,std::string(runtimeProvider())},{OpIds::EeOuterValidate,std::string(runtimeProvider())}};
        }
        result.fragment.id=id_;result.fragment.kind=PlanNodeKind::Sequence;
        for (const auto* op:ids) {
            auto leaf=operation(op,result.calls);
            result.fragment.children.push_back(std::move(leaf));
            result.operations.push_back({op,op,OperationStage::Prepare,
                std::string(op)==OpIds::EePressureSolve ? std::vector<OperationCapability>{OperationCapability::EulerianPhaseExecution,OperationCapability::PressureLinearSolve}
                    : std::vector<OperationCapability>{OperationCapability::EulerianPhaseExecution},actual.origin});
        }
        if (continuity) for (const auto* op:{OpIds::EeInterphaseCompute,OpIds::EeSourcesAssemble,OpIds::EeDtCompute,OpIds::EeStepBegin,OpIds::EeStepCommit,OpIds::EeBoundaryFinal,OpIds::EeOuterValidate})
            result.operations.push_back({op,op,OperationStage::Prepare,{OperationCapability::EulerianPhaseExecution},actual.origin});
        result.reads={"p"};result.writes={step.target.symbol};
        for (const auto* phase:phases) {
            const auto suffix=phase->id.substr(std::string("phaseMass").size());
            for (const auto* variable:{"phaseMass","momentum","enthalpy","alpha","rho","U","h","T"})
                result.reads.push_back(std::string(variable)+suffix);
            if (!pressure && !flux) {
                // Recovery republishes cached primitives and reference-phase closure into the existing arrays.
                for (const auto* variable:{"phaseMass","momentum","enthalpy","alpha","rho","U","h","T"})
                    result.writes.push_back(std::string(variable)+suffix);
            }
        }
        if (correction) result.writes.push_back("p");
        result.temporalMethod="physical-step segregated Eulerian update";
        result.requirements={"explicit shared-pressure PIMPLE HOW","phase Upwind; existing HYPRE backend","existing owner/COPY/SUM/canonical-face semantics"};
        result.providerContract=std::make_shared<const CompiledEulerianAssemblyContract>(std::move(contract));
        return result;
    }
private:std::string id_;
};
}
void addEulerianMethods(ProviderRegistry& providers) {
    addEulerianTurbulenceMethods(providers);
    static const EulerianMethod continuity("EulerianPhaseContinuity"),momentum("EulerianPhaseMomentum"),
        pressure("EulerianSharedPressureCorrection"),correction("EulerianPhaseCorrection"),
        flux("EulerianPhaseFluxCorrection"),energy("EulerianPhaseEnthalpy");
    for (const auto* method:{&continuity,&momentum,&pressure,&correction,&flux,&energy}) providers.add(*method);
}
}
