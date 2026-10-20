#pragma once
/// Public host interprets no model contract and owns no numerical workspace.
#include "core/interfaces/SF_solverStepper.h"
#include "core/config/SF_config.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"
namespace SF::SolverAlgorithm {
class SingleFluidStepper : public FDM::INavierStokesStepper {
public:
    SingleFluidStepper(FDM::SolverConfig,const System::ExecutableEquationSystem&,
        const System::CompiledNumericalSystem&,const System::CompiledSolvePlan&,
        const System::RuntimeRequirements&);
    // Callback captures require a stable host address for the entire binding lifetime.
    SingleFluidStepper(const SingleFluidStepper&)=delete;
    SingleFluidStepper& operator=(const SingleFluidStepper&)=delete;
    SingleFluidStepper(SingleFluidStepper&&)=delete;
    SingleFluidStepper& operator=(SingleFluidStepper&&)=delete;
    void bindServices(FDM::SolverServices) override;
    void bindSolvePlan(const System::CompiledSolvePlan&) override;
    void prepare(FDM::SolverState&) override;
    FDM::StepResult advance(FDM::SolverState&) override;
    double physicalTime() const;
    double timeStep() const;
private:
    FDM::SolverConfig config_;
    const System::ExecutableEquationSystem& executable_;
    const System::CompiledNumericalSystem& numerics_;
    const System::CompiledSolvePlan& solve_;
    const System::RuntimeRequirements& runtime_;
    FDM::SolverServices services_;
    State::StateBundle* state_=nullptr;
    // Identity-only invariant receipt: never used to enumerate or execute physical state.
    std::vector<Field*> boundPatchIdentities_;
    bool servicesBound_=false,solvePlanBound_=false,prepared_=false;
    double maximumTimeStep_=0;
    System::StateRealization realizedState_;
    Run::OpRegistry operations_;
    void emitTimeStep();
};
}
