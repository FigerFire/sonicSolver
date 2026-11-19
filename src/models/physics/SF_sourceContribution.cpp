/// @file SF_sourceContribution.cpp
/// @brief Gravity, MRF, and wall-heat mathematical source contributions.

#include "SF_sourceContribution.h"

#include "core/system/SF_systemContribution.h"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace SF::Physics::SourceContribution {
namespace {

struct SourceTarget {
    FDM::SourceKind kind;
    const char* symbol;
    bool energy;
};

constexpr std::array<SourceTarget,3> sourceKinds{{
    {FDM::SourceKind::Gravity,"gravity",false},
    {FDM::SourceKind::MRF,"MRF",false},
    {FDM::SourceKind::WallHeat,"wallHeat",true}
}};

} // namespace

void contribute(
        System::SystemContribution& system,
        const FDM::SourceConfig& config,const std::vector<Target>& destinations) {
    for (FDM::SourceKind kind : config.enabled) {
        const auto item = std::find_if(
            sourceKinds.begin(),sourceKinds.end(),
            [kind](const SourceTarget& value) { return value.kind == kind; });
        if (item == sourceKinds.end()) {
            throw std::runtime_error(
                "No equation contribution is registered for SourceKind.");
        }
        system.recordContribution(
            std::string("model.")+item->symbol,
            std::string(item->symbol)+" equation source contribution");
        bool applied=false;
        for(const auto& target:destinations) {
            const bool energy=target.form==Form::EnergyDensity;
            if(energy!=item->energy)continue;
            const auto units=energy?Math::signature(0,{1,-1,-3,0,0,0,0})
                :target.form==Form::Acceleration?Math::signature(1,{0,1,-2,0,0,0,0}):Math::signature(1,{1,-2,-2,0,0,0,0});
            auto source=System::FormulaExpr::op(Math::Op::Source,{System::FormulaExpr::symbol(item->symbol,units)},item->symbol);
            System::MathematicalExtension extension;extension.families={target.equation};
            extension.rhs=std::move(source);extension.exactEquation=true;
            system.mathematicalExtensions.push_back(std::move(extension));applied=true;
        }
        if(!applied)throw std::runtime_error(std::string(item->symbol)+": no explicit compatible source target equation");
    }
}
std::vector<Target> targets(const System::EquationRegistry& equations) {
    std::vector<Target> result;
    for(const auto& equation:equations.entries()) {
        // These are the existing preset's declared source roles, not field identities.
        bool role=false;for(const auto* family:{"momentum","energy","enthalpy","E_ENTHALPY"})
            role=role || equation.id==family || equation.id.rfind(std::string(family)+".",0)==0;
        if(!role)continue;
        std::vector<const System::FormulaExpr*> transients;
        const auto visit=[&](const auto& self,const System::FormulaExpr& expr)->void {
            if(expr.kind()==System::FormulaExpr::Kind::Operator && expr.operation()==Math::Op::Ddt)transients.push_back(&expr);
            for(const auto& child:expr.arguments())self(self,child);
        };visit(visit,equation.lhs);
        if(transients.size()!=1 || transients.front()->arguments().size()!=1 || transients.front()->arguments()[0].kind()!=System::FormulaExpr::Kind::Symbol)throw std::runtime_error("Source target requires an explicit typed form: "+equation.id);
        const auto signature=transients.front()->signature();
        if(Math::compatible(signature,Math::signature(1,{0,1,-2,0,0,0,0})))result.push_back({equation.id,Form::Acceleration});
        else if(Math::compatible(signature,Math::signature(1,{1,-2,-2,0,0,0,0})))result.push_back({equation.id,Form::ForceDensity});
        else if(Math::compatible(signature,Math::signature(0,{1,-1,-3,0,0,0,0})))result.push_back({equation.id,Form::EnergyDensity});
        else throw std::runtime_error("Unsupported source target units; supply an explicit typed form: "+equation.id);
    }
    return result;
}
} // namespace SF::Physics::SourceContribution
