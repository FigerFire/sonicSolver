/// @file SF_kkt.cpp
/// @brief 压力、速度、表面乘子与可选刚体广义速度的单体 KKT 装配。

#include "solver/algorithm/pressureBased/SF_kkt.h"

#include "SF_fluidStateModel.h"
#include "SF_immersedSystem.h"
#include "SF_thermodynamicClosure.h"
#include "constraint/variational/SF_constraintBlock.h"
#include "solver/algorithm/pressureBased/SF_kktDofs.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace SF::PressureBased {
namespace {
struct Coefficient { int index = -1; double value = 0.0; };
bool unknown(const Field& field, int i, int j, int k) {
    const int ng=field.NG();
    return i>=ng && i<ng+field.NX()
        && j>=ng && j<ng+field.NY()
        && k>=ng && k<ng+field.NZ()
        && field.CellFlag(i,j,k)==FLUID_CELL
        && !field.isSolverBoundaryPoint(i,j,k);
}
std::array<int,3> momentum(const Field& field) {
    if (field.hasStateModel()) return {
        field.stateModel()->momentumIndex(0),
        field.stateModel()->momentumIndex(1),
        field.stateModel()->momentumIndex(2)};
    return {RU,RV,RW};
}
double density(const Field& field, int i, int j, int k) {
    const int count=field.hasStateModel()
        ?field.stateModel()->densityVariableCount():1;
    double result=0.0;
    for (int variable=0; variable<count; ++variable)
        result+=field(i,j,k,variable);
    if (!std::isfinite(result)||result<=0.0)
        throw std::runtime_error("surface KKT found invalid density.");
    return result;
}
double volume(const Field& field, int i, int j, int k) {
    const double inverse=field.Jac(i,j,k);
    if (!std::isfinite(inverse)||inverse==0.0)
        throw std::runtime_error("surface KKT found invalid Jacobian.");
    return 1.0/std::abs(inverse);
}
std::array<double,3> metric(
        const Field& field,int axis,int i,int j,int k) {
    if(axis==0)return {field.XiX(i,j,k),field.XiY(i,j,k),field.XiZ(i,j,k)};
    if(axis==1)return {field.EtX(i,j,k),field.EtY(i,j,k),field.EtZ(i,j,k)};
    return {field.ZeX(i,j,k),field.ZeY(i,j,k),field.ZeZ(i,j,k)};
}
void offset(int axis,int sign,int& i,int& j,int& k) {
    i=j=k=0;
    if(axis==0)i=sign; else if(axis==1)j=sign; else k=sign;
}
Vector3 velocity(const Field& field,const std::array<int,3>& m,
                 int i,int j,int k) {
    const double rho=density(field,i,j,k);
    return {field(i,j,k,m[0])/rho,field(i,j,k,m[1])/rho,
            field(i,j,k,m[2])/rho};
}
double divergenceAt(const Field& field,const std::array<int,3>& m,
                    int i,int j,int k) {
    double result=0.0;
    for(int axis=0;axis<3;++axis){
        if(!Math::isDirectionActiveIndex(axis))continue;
        int ip=0,jp=0,kp=0,im=0,jm=0,km=0;
        offset(axis,1,ip,jp,kp);offset(axis,-1,im,jm,km);
        const Vector3 upper=velocity(field,m,i+ip,j+jp,k+kp);
        const Vector3 lower=velocity(field,m,i+im,j+jm,k+km);
        const auto met=metric(field,axis,i,j,k);
        result+=0.5*(met[0]*(upper.x-lower.x)
                    +met[1]*(upper.y-lower.y)
                    +met[2]*(upper.z-lower.z));
    }
    return result;
}
double component(const Vector3& value,int index) {
    return index==0?value.x:(index==1?value.y:value.z);
}
void setComponent(Vector3& value,int index,double componentValue) {
    if(index==0)value.x=componentValue;
    else if(index==1)value.y=componentValue;
    else value.z=componentValue;
}
void add(LinearAlgebra::GlobalDofRow& row,
         LinearAlgebra::GlobalDofId column,double value) {
    if(value==0.0)return;
    row.add(column,value);
}
} // namespace

MonolithicKKT::MonolithicKKT(
        FDM::PressureCorrectionConfig config,double idealGasGamma)
    :config_(std::move(config)),idealGasGamma_(idealGasGamma) {}

void MonolithicKKT::prepareLinearSystem(
        const std::vector<LinearAlgebra::GlobalDofId>& dofs) {
    if(dofs.empty())
        throw std::runtime_error("surface KKT has an empty linear system.");
    if(linearSystem_&&numberedDofs_==dofs)return;
    std::vector<std::pair<LinearAlgebra::GlobalDofId,std::int64_t>> entries;
    entries.reserve(dofs.size());
    for(std::size_t row=0;row<dofs.size();++row)
        entries.emplace_back(dofs[row],static_cast<std::int64_t>(row));
    numbering_=std::make_unique<LinearAlgebra::StaticDistributedNumbering>(
        0,static_cast<std::int64_t>(dofs.size())-1,
        static_cast<std::int64_t>(dofs.size()),entries);
    linearSystem_=std::make_unique<LinearAlgebra::DistributedLinearSystem>(
        config_.linear.pressure,*numbering_);
    numberedDofs_=dofs;
}

KKTCorrectionSummary MonolithicKKT::correct(
        Field& field,double targetTime,double dt,
        FDM::IImmersedConstraint& constraint) {
    if(!std::isfinite(dt)||dt<=0.0)
        throw std::runtime_error("surface KKT requires positive finite dt.");
    if(config_.linear.pressure.method!=FDM::KrylovMethod::FlexGMRES)
        throw std::runtime_error(
            "surface KKT is nonsymmetric and requires pressure solver "
            "flexGMRES in system/solverProperties.");
    if(config_.coupling.momentumRelaxation!=1.0
        ||config_.coupling.pressureRelaxation!=1.0)
        throw std::runtime_error(
            "surface KKT requires momentumRelaxation=1 and "
            "pressureRelaxation=1 so the solved constraints are not altered.");

    const FDM::IImmersedSystem* immersed=constraint.systemProvider();
    if(!immersed) {
        throw std::runtime_error(
            "surface KKT requires an immersed-system descriptor; the "
            "constraint adapter did not provide one.");
    }
    const auto& descriptor=immersed->algorithmDescriptor();
    if(!descriptor.monolithic
       ||!descriptor.variational.fluidKineticIncrement
       ||!descriptor.variational.incompressibilityConstraint
       ||!descriptor.variational.immersedNoSlipConstraint) {
        throw std::runtime_error(
            "surface KKT requires a monolithic constrained variational "
            "descriptor with kinetic, pressure and IBM no-slip terms.");
    }

    const auto& surface=constraint.prepareMonolithicSystem(field,targetTime,dt);
    const auto m=momentum(field);
    std::vector<int> cells;
    std::vector<int> local((size_t)field.TotalSize(),-1);
    const int ng=field.NG();
    for(int k=ng;k<ng+field.NZ();++k)
        for(int j=ng;j<ng+field.NY();++j)
            for(int i=ng;i<ng+field.NX();++i)
                if(unknown(field,i,j,k)){
                    local[(size_t)field.getIdx(i,j,k)]=(int)cells.size();
                    cells.push_back(field.getIdx(i,j,k));
                }
    if(cells.empty()||surface.points.empty())
        throw std::runtime_error("surface KKT has no fluid or surface DOFs.");

    const auto localCell = [&](int fieldIndex) {
        if (fieldIndex < 0
            || fieldIndex >= static_cast<int>(local.size())
            || local[static_cast<size_t>(fieldIndex)] < 0) {
            throw std::runtime_error(
                "surface KKT received a surface interpolation cell that is "
                "not a local fluid unknown.");
        }
        return local[static_cast<size_t>(fieldIndex)];
    };

    const int n=(int)cells.size();
    if(config_.reference.referenceCell<0
       ||config_.reference.referenceCell>=n)
        throw std::runtime_error(
            "surface KKT referenceCell is outside the pressure block.");
    const int markers=(int)surface.points.size();
    const KKTSubspace subspace=activeKKTSubspace(surface.solidEquation);
    const auto& velocityComponents=subspace.velocityComponents;
    const auto& solidComponents=subspace.solidComponents;
    const int vectorCount=(int)velocityComponents.size();
    const int qCount=(int)solidComponents.size();
    std::vector<std::int64_t> constraintEntities;
    constraintEntities.reserve(surface.points.size());
    for(const auto& point:surface.points) {
        if (!point.globalConstraintId.valid()) {
            throw std::runtime_error(
                "surface KKT received an invalid GlobalConstraintDofId.");
        }
        constraintEntities.push_back(point.globalConstraintId.value());
    }
    const KKTLayout layout(n,std::move(constraintEntities),
                           velocityComponents,solidComponents);
    std::vector<double> rho((size_t)n),cellVolume((size_t)n);
    std::vector<Vector3> predicted((size_t)n);
    std::vector<double> divergenceStar((size_t)n,0.0);
    std::vector<std::vector<Coefficient>> divergenceRows((size_t)n);
    std::vector<std::vector<Coefficient>> pressureTranspose((size_t)(vectorCount*n));
    KKTCorrectionSummary summary;

    for(int row=0;row<n;++row){
        int i=0,j=0,k=0;field.getIJK(cells[(size_t)row],i,j,k);
        rho[(size_t)row]=density(field,i,j,k);
        cellVolume[(size_t)row]=volume(field,i,j,k);
        predicted[(size_t)row]=velocity(field,m,i,j,k);
        const double div=divergenceAt(field,m,i,j,k);
        for(int axis=0;axis<3;++axis){
            if(!Math::isDirectionActiveIndex(axis))continue;
            const auto met=metric(field,axis,i,j,k);
            for(int componentSlot=0;componentSlot<vectorCount;++componentSlot){
                const int componentIndex=velocityComponents[(size_t)componentSlot];
                for(int sign:{-1,1}){
                    int di=0,dj=0,dk=0;offset(axis,sign,di,dj,dk);
                    const int neighbour=field.getIdx(i+di,j+dj,k+dk);
                    const int neighbourLocal=local[(size_t)neighbour];
                    if(neighbourLocal<0)continue;
                    const double value=0.5*sign
                        *met[(size_t)componentIndex];
                    const int dof=vectorCount*neighbourLocal+componentSlot;
                    divergenceRows[(size_t)row].push_back({dof,value});
                    pressureTranspose[(size_t)dof].push_back({row,
                        -cellVolume[(size_t)row]*value});
                }
            }
        }
        summary.maxDivergenceBefore=std::max(
            summary.maxDivergenceBefore,std::abs(div));
        divergenceStar[(size_t)row]=div;
    }

    for (const auto& point:surface.points) {
        if(!std::isfinite(point.prescribedVelocity.x)
           ||!std::isfinite(point.prescribedVelocity.y)
           ||!std::isfinite(point.prescribedVelocity.z))
            throw std::runtime_error(
                "surface KKT received an invalid prescribed velocity.");
    }
    // J/J^T 的 marker 拓扑属于 IBM variational layer；pressure KKT 只将
    // 这个 block 映射到自己的 fluid/pressure rows。
    const auto constraintBlock=IBM::Variational::SurfaceConstraintBlock::build(
        surface,velocityComponents,[&](int cell) { return localCell(cell); });
    std::vector<std::vector<Coefficient>> spread(
        (size_t)(vectorCount*n));
    for (const auto& row:constraintBlock.rows()) {
        for (const auto& term:row) {
            spread[(size_t)(vectorCount*term.localCell+term.component)].push_back(
                {vectorCount*term.marker+term.component,term.adjoint});
        }
    }

    // 增广 Lagrangian 只改变原 KKT 的 fluid Hessian 和已知约束误差右端：
    // (A + gamma J^T M_L J) du - J^T M_L Lambda
    //     = gamma J^T M_L (U_s - J u*).
    // Lambda 仍是未知量，因此这不是 Brinkman/penalty 替代。
    std::vector<std::vector<Coefficient>> augmentation(
        (size_t)(vectorCount*n));
    std::vector<double> augmentationRhs(
        (size_t)(vectorCount*n),0.0);
    const double gamma=surface.augmentationCoefficient;
    if (!std::isfinite(gamma) || gamma<0.0) {
        throw std::runtime_error(
            "surface KKT received an invalid augmentation coefficient.");
    }
    if (gamma>0.0 && surface.solidEquation.solveGeneralizedVelocity) {
        throw std::runtime_error(
            "coupled-rigid augmented KKT requires u-q augmentation cross "
            "blocks, which are not implemented in this serial build.");
    }
    if (gamma>0.0) {
        for (int marker=0;marker<markers;++marker) {
            const auto& point=surface.points[static_cast<std::size_t>(marker)];
            Vector3 interpolated;
            for (int c=0;c<vectorCount;++c) {
                const int componentIndex=velocityComponents[(size_t)c];
                double value=0.0;
                for (const auto& term:constraintBlock.row(
                        marker,c)) {
                    value += component(predicted[(size_t)term.localCell],
                                       componentIndex)*term.interpolation;
                }
                setComponent(interpolated,componentIndex,value);
            }
            const Vector3 error=point.prescribedVelocity-interpolated;
            for (int c=0;c<vectorCount;++c) {
                    const int componentIndex=velocityComponents[(size_t)c];
                    for (const auto& rowTerm:constraintBlock.row(marker,c)) {
                        const int rowDof=vectorCount*rowTerm.localCell+c;
                        augmentationRhs[(size_t)rowDof] += gamma*point.measure
                            *rowTerm.interpolation*component(error,componentIndex);
                        for (const auto& columnTerm:
                             constraintBlock.row(marker,c)) {
                        augmentation[(size_t)rowDof].push_back({
                            vectorCount*columnTerm.localCell+c,
                            gamma*point.measure*rowTerm.interpolation
                                *columnTerm.interpolation});
                        }
                    }
            }
        }
    }

    prepareLinearSystem(layout.orderedDofs());
    LinearAlgebra::GlobalDofSystem matrix;
    matrix.rows.reserve(layout.orderedDofs().size());
    for(int row=0;row<n;++row){
        LinearAlgebra::GlobalDofRow equation(
            layout.pressure(row));
        int i=0,j=0,k=0;field.getIJK(cells[(size_t)row],i,j,k);
        if(row==config_.reference.referenceCell){
            add(equation,layout.pressure(row),1.0);
            equation.setRightHandSide(0.0);
        }else{
            for(const auto& item:divergenceRows[(size_t)row]){
                add(equation,layout.velocity(
                    item.index/vectorCount,
                    velocityComponents[(size_t)(item.index%vectorCount)]),
                    cellVolume[(size_t)row]*item.value);
            }
            const double soundSquared=field.hasStateModel()
                ?std::pow(field.thermodynamicState(i,j,k).soundSpeed,2)
                :idealGasGamma_*Boundary::pressureAt(field,i,j,k)
                    /rho[(size_t)row];
            if(!std::isfinite(soundSquared)||soundSquared<=0.0)
                throw std::runtime_error("surface KKT found invalid sound speed.");
            add(equation,layout.pressure(row),
                cellVolume[(size_t)row]/(rho[(size_t)row]*soundSquared*dt));
            equation.setRightHandSide(
                -cellVolume[(size_t)row]*divergenceStar[(size_t)row]);
        }
        matrix.rows.push_back(std::move(equation));
    }
    for(int dof=0;dof<vectorCount*n;++dof){
        const int cell=dof/vectorCount;
        const int componentIndex=
            velocityComponents[(size_t)(dof%vectorCount)];
        LinearAlgebra::GlobalDofRow equation(
            layout.velocity(cell,componentIndex));
        add(equation,layout.velocity(cell,componentIndex),
            rho[(size_t)cell]*cellVolume[(size_t)cell]/dt);
        for(const auto& item:augmentation[(size_t)dof])
            add(equation,layout.velocity(
                item.index/vectorCount,
                velocityComponents[(size_t)(item.index%vectorCount)]),item.value);
        for(const auto& item:pressureTranspose[(size_t)dof])
            add(equation,layout.pressure(item.index),item.value);
        for(const auto& item:spread[(size_t)dof])
            add(equation,layout.constraint(
                item.index/vectorCount,
                velocityComponents[(size_t)(item.index%vectorCount)]),item.value);
        equation.setRightHandSide(augmentationRhs[(size_t)dof]);
        matrix.rows.push_back(std::move(equation));
    }
    for(int marker=0;marker<markers;++marker){
        const auto& point=surface.points[(size_t)marker];
        Vector3 interpolated;
        for(int c=0;c<vectorCount;++c){
            const int componentIndex=velocityComponents[(size_t)c];
            double value=0.0;
            for(const auto& term:constraintBlock.row(marker,c))
                value += component(
                    predicted[(size_t)term.localCell],componentIndex)
                    *term.interpolation;
            setComponent(interpolated,componentIndex,value);
        }
        const Vector3 target=point.prescribedVelocity;
        for(int c=0;c<vectorCount;++c){
            const int componentIndex=velocityComponents[(size_t)c];
            LinearAlgebra::GlobalDofRow equation(
                layout.constraint(marker,componentIndex));
            for(const auto& term:constraintBlock.row(marker,c)){
                add(equation,layout.velocity(
                    term.localCell,componentIndex),term.interpolation);
            }
            if(surface.solidEquation.solveGeneralizedVelocity){
                for(int q=0;q<qCount;++q){
                    const int solidComponent=solidComponents[(size_t)q];
                    add(equation,layout.solid(solidComponent),
                        -component(point.solidVelocityBasis[
                            (size_t)solidComponent],componentIndex));
                }
            }
            equation.setRightHandSide(
                component(target-interpolated,componentIndex));
            matrix.rows.push_back(std::move(equation));
        }
    }
    if(surface.solidEquation.solveGeneralizedVelocity){
        for(int solidRow=0;solidRow<qCount;++solidRow){
            const int rowComponent=solidComponents[(size_t)solidRow];
            LinearAlgebra::GlobalDofRow equation(
                layout.solid(rowComponent));
            for(int q=0;q<qCount;++q){
                const int columnComponent=solidComponents[(size_t)q];
                add(equation,layout.solid(columnComponent),
                    surface.solidEquation.lhs[
                        (size_t)(6*rowComponent+columnComponent)]);
            }
            for(int marker=0;marker<markers;++marker){
                const auto& point=surface.points[(size_t)marker];
                for(int lambdaSlot=0;lambdaSlot<vectorCount;++lambdaSlot){
                    const int lambdaComponent=
                        velocityComponents[(size_t)lambdaSlot];
                    const double value=point.measure*component(
                        point.solidVelocityBasis[(size_t)rowComponent],
                        lambdaComponent);
                    add(equation,layout.constraint(
                        marker,lambdaComponent),value);
                }
            }
            equation.setRightHandSide(
                surface.solidEquation.rhs[(size_t)rowComponent]);
            equation.setInitialGuess(
                surface.solidEquation.initialGuess[(size_t)rowComponent]);
            matrix.rows.push_back(std::move(equation));
        }
    }

    const auto solution=linearSystem_->solve(matrix);
    summary.iterations=solution.iterations;
    summary.relativeResidual=solution.relativeResidual;
    summary.structureRebuilds=linearSystem_->statistics().structureRebuilds;
    summary.linearSolves=linearSystem_->statistics().solves;
    double kineticIncrementFunctional=0.0;
    for(int cell=0;cell<n;++cell){
        int i=0,j=0,k=0;field.getIJK(cells[(size_t)cell],i,j,k);
        Vector3 corrected=predicted[(size_t)cell];
        Vector3 increment;
        for(int c:velocityComponents)setComponent(corrected,c,
            component(corrected,c)+solution.solution[
                layout.velocitySlot(cell,c)]);
        for(int c:velocityComponents)setComponent(increment,c,
            solution.solution[layout.velocitySlot(cell,c)]);
        kineticIncrementFunctional += 0.5*rho[(size_t)cell]
            *cellVolume[(size_t)cell]*dot(increment,increment)/dt;
        for(int c:velocityComponents)
            field(i,j,k,m[(size_t)c])=rho[(size_t)cell]*component(corrected,c);
        const double pressure=Boundary::pressureAt(field,i,j,k)
            +solution.solution[layout.pressureSlot(cell)];
        if(!std::isfinite(pressure)||pressure<=0.0)
            throw std::runtime_error("surface KKT produced non-positive pressure.");
        if(field.hasStateModel()){
            std::vector<double> state((size_t)field.NVar());
            for(int variable=0;variable<field.NVar();++variable)
                state[(size_t)variable]=field(i,j,k,variable);
            field(i,j,k,field.stateModel()->energyIndex())=
                field.stateModel()->totalEnergyFromPressure(
                    state.data(),field.NVar(),pressure);
        }else{
            const double kinetic=0.5*rho[(size_t)cell]*dot(corrected,corrected);
            field(i,j,k,E)=pressure/(idealGasGamma_-1.0)+kinetic;
        }
    }
    field.invalidateThermodynamicCache();
    FDM::ImmersedKKTState state;
    state.kineticIncrementFunctional=kineticIncrementFunctional;
    for(int dof=0;dof<vectorCount*n;++dof){
        const std::size_t matrixRow=layout.velocitySlot(
            dof/vectorCount,
            velocityComponents[(size_t)(dof%vectorCount)]);
        const auto& row=matrix.rows[matrixRow];
        double residual=-row.rightHandSide();
        for(const auto& entry:row.coefficients()){
            residual += entry.value*solution.solution[
                static_cast<std::size_t>(numbering_->row(entry.column))];
        }
        state.maximumStationarityResidual=std::max(
            state.maximumStationarityResidual,std::abs(residual));
    }
    state.surfaceMultiplier.resize((size_t)markers);
    for(int marker=0;marker<markers;++marker)
        for(int c:velocityComponents)setComponent(
            state.surfaceMultiplier[(size_t)marker],c,
            solution.solution[layout.constraintSlot(marker,c)]);
    if(surface.solidEquation.solveGeneralizedVelocity){
        state.hasGeneralizedVelocity=true;
        state.generalizedVelocity=surface.solidEquation.initialGuess;
        for(int q:solidComponents)state.generalizedVelocity[(size_t)q]=
            solution.solution[layout.solidSlot(q)];
    }
    summary.immersed=constraint.acceptMonolithicSolution(
        field,targetTime,dt,state);
    for(int row=0;row<n;++row){
        int i=0,j=0,k=0;field.getIJK(cells[(size_t)row],i,j,k);
        const double div=divergenceAt(field,m,i,j,k);
        summary.maxDivergenceAfter=std::max(
            summary.maxDivergenceAfter,std::abs(div));
    }
    return summary;
}

} // namespace SF::PressureBased
