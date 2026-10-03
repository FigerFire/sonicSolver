#include "SF_formulaStorage.h"

#include <stdexcept>
#include <utility>

namespace SF::System {

FormulaValues bindFormulaCellValues(
        const StateRealization& state, const CompiledTarget& output,
        std::function<double(const std::string&,int,int,int,int)> boundaryValue) {
    if (output.symbol.empty() || output.kind==TargetKind::Workspace)
        throw std::runtime_error("Cell Equation binding requires one declared Output symbol.");
    const auto& target=state.at(output.symbol);
    const bool physical=output.kind==TargetKind::Physical;
    const auto outputFields=physical ? target.fields
        : state.view(output.symbol,output.kind==TargetKind::Working
            ? StateViewKind::Working : StateViewKind::Correction).fields;
    const int outputOffset=physical ? target.descriptor->componentOffset : 0;
    if (!target.descriptor || target.descriptor->location!=VariableLocation::EulerianCell
        || outputFields.size()!=1 || !outputFields.front()
        || !outputFields.front()->geometry)
        throw std::runtime_error("Generic Equation cell binding requires one realized cell patch.");
    const Field& geometry=*outputFields.front()->geometry;
    const int nx=geometry.MX()-2*geometry.NG();
    const int ny=geometry.MY()-2*geometry.NG();
    const int nz=geometry.MZ()-2*geometry.NG();
    if (nx<=0 || ny<=0 || nz<=0 || !boundaryValue)
        throw std::runtime_error("Generic Equation cell binding requires valid geometry and boundary closure.");
    const auto fullCell=[field=&geometry,nx,ny](int cell) {
        const int i=cell%nx+field->NG();
        const int j=(cell/nx)%ny+field->NG();
        const int k=cell/(nx*ny)+field->NG();
        return field->getIdx(i,j,k);
    };
    FormulaValues values;
    values.ownedCells=nx*ny*nz;
    values.components=target.descriptor->components;
    values.read=[&state,fullCell,count=values.ownedCells,outputFields,outputOffset,targetId=output.symbol](
            const std::string& symbol,int cell,int component) {
        if (cell<0 || cell>=count) throw std::out_of_range("Equation cell read is out of range.");
        if (symbol==targetId) {
            const int c=outputFields.front()->components==1 ? 0 : component;
            return outputFields.front()->read(fullCell(cell),outputOffset+c);
        }
        const auto& resolved=state.at(symbol);
        if (resolved.descriptor->constantValue)
            return *resolved.descriptor->constantValue;
        if (resolved.fields.size()!=1 || !resolved.fields.front())
            throw std::runtime_error("Equation symbol '"+symbol+"' lacks single-patch storage.");
        const int c=resolved.descriptor->components==1 ? 0 : component;
        if (c<0 || c>=resolved.descriptor->components)
            throw std::runtime_error("Equation symbol component is out of range.");
        return resolved.fields.front()->read(
            fullCell(cell),resolved.descriptor->componentOffset+c);
    };
    values.write=[&state,fullCell,count=values.ownedCells,
                  targetId=output.symbol,outputFields,outputOffset,physical](
            const std::string& symbol,int cell,int component,double value) {
        if (symbol!=targetId || cell<0 || cell>=count
            || component<0 || component>=state.at(symbol).descriptor->components)
            throw std::runtime_error("Equation write violates declared Output layout.");
        outputFields.front()->write(fullCell(cell),outputOffset+component,value);
        // A generic Equation may update the conservative Field. Derived EOS
        // values must not survive a write through this non-owning view.
        if (physical) outputFields.front()->geometry->invalidateThermodynamicCache();
    };
    values.dof=[targetId=output.symbol,components=values.components,
                count=values.ownedCells](
            const std::string& symbol,int cell,int component) {
        if (symbol!=targetId || cell<0 || cell>=count
            || component<0 || component>=components)
            throw std::runtime_error("Equation DOF violates declared Output.");
        return LinearAlgebra::GlobalDofId::make(
            LinearAlgebra::GlobalDofSpace::ScalarTransport,cell,component);
    };
    values.boundaryValue=std::move(boundaryValue);
    return values;
}

} // namespace SF::System
