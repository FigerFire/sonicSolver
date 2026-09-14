/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.11-----------*/

/// @file SF_densityBasedTime.cpp
/// @brief Density-based patches and registered scalars share one explicit tableau.
///
/// Qn remains owned by StateBundle.  This file only constructs temporary RK
/// stage state and combines existing RHS snapshots.  The supplied RHS callback
/// completes every patch's spatial work and its canonical/MPI barriers before
/// a stage update; it is deliberately the only place that knows those details.

#include "solver/algorithm/SF_densityBasedTime.h"
#include "solver/algorithm/SF_highOrderTrace.h"

#include "methods/numerics/structured/SF_structured.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::SolverAlgorithm::DensityBasedTime {
namespace {

struct ScalarRKStorage {
    std::vector<double> q0, k1, k2, k3, k4;
};

struct RKStorage {
    std::vector<double> q0, k1, k2, k3, k4;
    std::vector<ScalarRKStorage> registered;
};

bool usesExplicitTableau(State::UpdatePolicy policy) {
    return policy == State::UpdatePolicy::Explicit
        || policy == State::UpdatePolicy::BoundedExplicit
        || policy == State::UpdatePolicy::HamiltonJacobi;
}

void traceUpdatedConservative(
        const char* checkpoint,
        const std::vector<Field*>& fields,
        const State::StateBundle& state) {
    if (HighOrderTrace::activeFor(state.step)) {
        HighOrderTrace::conservative(checkpoint, fields);
    }
}

void applyScalarStage(ScalarField& scalar,
                      const std::vector<double>& q0,
                      const std::vector<double>& rhs,
                      double coeff,
                      double dt) {
    if (q0.size() != scalar.values().size() || rhs.size() != q0.size()) {
        throw std::runtime_error("integrated scalar stage size mismatch.");
    }
    for (size_t n = 0; n < q0.size(); ++n) {
        scalar.values()[n] = q0[n] + coeff * dt * rhs[n];
    }
}

void applyScalarFinal(ScalarField& scalar,
                      const ScalarRKStorage& storage,
                      double dt) {
    const size_t size = storage.q0.size();
    if (size != scalar.values().size()
        || storage.k1.size() != size || storage.k2.size() != size
        || storage.k3.size() != size || storage.k4.size() != size) {
        throw std::runtime_error("registered scalar RK4 size mismatch.");
    }
    for (size_t n = 0; n < size; ++n) {
        scalar.values()[n] = storage.q0[n]
            + dt/6.0*(storage.k1[n]+2.0*storage.k2[n]
                      +2.0*storage.k3[n]+storage.k4[n]);
    }
}

void validateRegistered(const State::IntegratedVariable& variable,
                        const char* stage) {
    const auto& bounds = variable.descriptor.bounds;
    for (double value : variable.value->values()) {
        if (!std::isfinite(value)) {
            throw std::runtime_error(
                "Transported variable '" + variable.descriptor.name
                + "' is non-finite after " + stage + ".");
        }
        if (bounds.enabled && (value < bounds.lower || value > bounds.upper)) {
            throw std::runtime_error(
                "Transported variable '" + variable.descriptor.name
                + "' violates bounds after " + stage
                + "; multi-block integration does not clamp it.");
        }
    }
}

void snapshotRegistered(State::VariableRegistry& registry,
                        const Field& field,
                        std::vector<ScalarRKStorage>& storage) {
    registry.validateFor(field);
    storage.resize(registry.size());
    for (size_t n = 0; n < registry.size(); ++n) {
        const auto policy =
            registry.variables()[n].descriptor.updatePolicy;
        if (usesExplicitTableau(policy)) {
            storage[n].q0 = registry.variables()[n].value->values();
        }
    }
}

void captureRegisteredRHS(State::VariableRegistry& registry,
                          std::vector<ScalarRKStorage>& storage,
                          int stage) {
    if (storage.size() != registry.size()) {
        throw std::runtime_error("transported-variable registry changed during RK4.");
    }
    for (size_t n = 0; n < registry.size(); ++n) {
        const auto policy =
            registry.variables()[n].descriptor.updatePolicy;
        if (!usesExplicitTableau(policy)) {
            continue;
        }
        std::vector<double>* destination = nullptr;
        if (stage == 1) destination = &storage[n].k1;
        else if (stage == 2) destination = &storage[n].k2;
        else if (stage == 3) destination = &storage[n].k3;
        else if (stage == 4) destination = &storage[n].k4;
        else throw std::runtime_error("invalid registered RK stage.");
        *destination = *registry.variables()[n].rhs;
    }
}

void applyRegisteredStage(State::VariableRegistry& registry,
                          const std::vector<ScalarRKStorage>& storage,
                          int stage, double coeff, double dt,
                          const char* label) {
    for (size_t n = 0; n < registry.size(); ++n) {
        const auto policy =
            registry.variables()[n].descriptor.updatePolicy;
        if (!usesExplicitTableau(policy)) {
            continue;
        }
        const std::vector<double>* rhs = nullptr;
        if (stage == 1) rhs = &storage[n].k1;
        else if (stage == 2) rhs = &storage[n].k2;
        else if (stage == 3) rhs = &storage[n].k3;
        else throw std::runtime_error("invalid registered RK stage update.");
        applyScalarStage(*registry.variables()[n].value,
                         storage[n].q0, *rhs, coeff, dt);
        validateRegistered(registry.variables()[n], label);
    }
}

void applyRegisteredFinal(State::VariableRegistry& registry,
                          const std::vector<ScalarRKStorage>& storage,
                          double dt) {
    for (size_t n = 0; n < registry.size(); ++n) {
        const auto policy =
            registry.variables()[n].descriptor.updatePolicy;
        if (!usesExplicitTableau(policy)) {
            continue;
        }
        applyScalarFinal(*registry.variables()[n].value, storage[n], dt);
        validateRegistered(registry.variables()[n], "RK4 final");
    }
}

size_t storageIndex(const Field& field, int i, int j, int k, int v) {
    return (((size_t)k * field.MY() + (size_t)j) * field.MX() + (size_t)i)
        * (size_t)field.NVar() + (size_t)v;
}

void snapshotField(const Field& field, std::vector<double>& dst) {
    dst.assign((size_t)field.MX() * field.MY() * field.MZ()
               * (size_t)field.NVar(), 0.0);
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            dst[storageIndex(field, i, j, k, v)] = field(i, j, k, v);
        }
    });
}

void snapshotRHS(const Field& field, const Residual& residual,
                 std::vector<double>& dst) {
    dst.assign((size_t)field.MX() * field.MY() * field.MZ()
               * (size_t)field.NVar(), 0.0);
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            dst[storageIndex(field, i, j, k, v)] =
                Math::spatialResidual(field, residual, i, j, k, v);
        }
    });
}

void applyStage(Field& field,
                const std::vector<double>& q0,
                const std::vector<double>& kn,
                double coeff,
                double dt) {
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            const size_t idx = storageIndex(field, i, j, k, v);
            field(i, j, k, v) = q0[idx] - coeff * dt * kn[idx];
        }
    });
}

void applyFinalRK4(Field& field,
                   const RKStorage& storage,
                   double dt) {
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            const size_t idx = storageIndex(field, i, j, k, v);
            field(i, j, k, v) = storage.q0[idx]
                - dt / 6.0 * (storage.k1[idx]
                              + 2.0 * storage.k2[idx]
                              + 2.0 * storage.k3[idx]
                              + storage.k4[idx]);
        }
    });
}

void applySSPStage(Field& field,
                   const std::vector<double>& q0,
                   const std::vector<double>& rhs,
                   double baseWeight,
                   double eulerWeight,
                   double dt) {
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            const size_t idx = storageIndex(field, i, j, k, v);
            const double euler = field(i, j, k, v) - dt * rhs[idx];
            field(i, j, k, v) = baseWeight * q0[idx]
                + eulerWeight * euler;
        }
    });
}

void applyRegisteredSSPStage(
        State::VariableRegistry& registry,
        const std::vector<ScalarRKStorage>& storage,
        double baseWeight,
        double eulerWeight,
        double dt,
        const char* label) {
    for (size_t n = 0; n < registry.size(); ++n) {
        const auto& variable = registry.variables()[n];
        if (!usesExplicitTableau(variable.descriptor.updatePolicy)) continue;
        auto& values = variable.value->values();
        if (!variable.rhs || storage[n].q0.size() != values.size()
            || variable.rhs->size() != values.size()) {
            throw std::runtime_error(
                "multi-block SSP-RK3 storage mismatch for '"
                + variable.descriptor.name + "'.");
        }
        for (size_t id = 0; id < values.size(); ++id) {
            const double euler = values[id] + dt * (*variable.rhs)[id];
            values[id] = baseWeight * storage[n].q0[id]
                + eulerWeight * euler;
        }
        validateRegistered(variable, label);
    }
}

State::VariableRegistry* registryFor(
        Field& field,
        const std::vector<Field*>& fields,
        State::StateBundle& state,
        FDM::IEquationSystemCoupling* equationSystem) {
    if (equationSystem) {
        if (auto* registry = equationSystem->variables(field)) return registry;
    }
    if (fields.size() == 1 && !state.transported.empty()) {
        return &state.transported;
    }
    return nullptr;
}

} // namespace

void stepEuler(const std::vector<Field*>& fields,
               std::vector<PatchWorkspace>& workspaces,
               State::StateBundle& state,
               FDM::IEquationSystemCoupling* equationSystem,
               const AssembleRHS& assembleRHS,
               const Publish& publish,
               const Validate& validate) {
    assembleRHS(fields, workspaces, state.time);
    for (size_t index = 0; index < fields.size(); ++index) {
        Field* field = fields[index];
        if (!field) continue;
        Math::applyDivergence(*field, workspaces[index].residual, state.dt);
        field->invalidateThermodynamicCache();
        if (State::VariableRegistry* registry =
                registryFor(*field, fields, state, equationSystem)) {
                registry->validateFor(*field);
                for (const auto& variable : registry->variables()) {
                    const auto policy = variable.descriptor.updatePolicy;
                    if (!usesExplicitTableau(policy)) {
                        continue;
                    }
                    const std::vector<double> old = variable.value->values();
                    applyScalarStage(*variable.value, old, *variable.rhs,
                                     1.0, state.dt);
                    validateRegistered(variable, "Euler final");
                }
        }
    }
    traceUpdatedConservative("Q after Euler update", fields, state);
    publish(fields);
    validate(fields, "Euler final");
}

void stepSSPRK3(const std::vector<Field*>& fields,
                std::vector<PatchWorkspace>& workspaces,
                State::StateBundle& state,
                FDM::IEquationSystemCoupling* equationSystem,
                const AssembleRHS& assembleRHS,
                const Publish& publish,
                const Validate& validate) {
    std::vector<RKStorage> storage(fields.size());
    for (size_t n = 0; n < fields.size(); ++n) {
        if (!fields[n]) continue;
        snapshotField(*fields[n], storage[n].q0);
        if (auto* registry = registryFor(
                *fields[n], fields, state, equationSystem)) {
            snapshotRegistered(*registry, *fields[n], storage[n].registered);
        }
    }

    const auto advanceStage = [&](double stageFraction,
                                  double baseWeight,
                                  double eulerWeight,
                                  const char* label) {
        assembleRHS(fields, workspaces, state.time + stageFraction * state.dt);
        for (size_t n = 0; n < fields.size(); ++n) {
            Field* field = fields[n];
            if (!field) continue;
            snapshotRHS(*field, workspaces[n].residual, storage[n].k1);
            applySSPStage(*field, storage[n].q0, storage[n].k1,
                          baseWeight, eulerWeight, state.dt);
            field->invalidateThermodynamicCache();
            if (auto* registry = registryFor(*field, fields, state, equationSystem)) {
                applyRegisteredSSPStage(
                    *registry, storage[n].registered,
                    baseWeight, eulerWeight, state.dt, label);
            }
        }
        traceUpdatedConservative(label, fields, state);
        publish(fields);
        validate(fields, label);
    };

    advanceStage(0.0, 0.0, 1.0, "SSP-RK3 stage 1");
    advanceStage(1.0, 0.75, 0.25, "SSP-RK3 stage 2");
    advanceStage(0.5, 1.0 / 3.0, 2.0 / 3.0, "SSP-RK3 final");
}

void stepRK4(const std::vector<Field*>& fields,
             std::vector<PatchWorkspace>& workspaces,
             State::StateBundle& state,
             FDM::IEquationSystemCoupling* equationSystem,
             const AssembleRHS& assembleRHS,
             const Publish& publish,
             const Validate& validate) {
    std::vector<RKStorage> storage(fields.size());
    for (size_t n = 0; n < fields.size(); ++n) {
        if (!fields[n]) continue;
        snapshotField(*fields[n], storage[n].q0);
        if (auto* registry = registryFor(
                *fields[n], fields, state, equationSystem)) {
            snapshotRegistered(*registry, *fields[n], storage[n].registered);
        }
    }

    assembleRHS(fields, workspaces, state.time);
    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) snapshotRHS(*fields[n], workspaces[n].residual, storage[n].k1);
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                captureRegisteredRHS(*registry, storage[n].registered, 1);
        }
    }
    validate(fields, "RK4 stage 1");

    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) applyStage(*fields[n], storage[n].q0, storage[n].k1, 0.5, state.dt);
        if (fields[n]) fields[n]->invalidateThermodynamicCache();
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                applyRegisteredStage(*registry, storage[n].registered,
                                     1, 0.5, state.dt, "RK4 stage 2");
        }
    }
    traceUpdatedConservative("Q after RK4 stage 2 update", fields, state);
    publish(fields);
    validate(fields, "RK4 stage 2");
    assembleRHS(fields, workspaces, state.time + 0.5 * state.dt);
    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) snapshotRHS(*fields[n], workspaces[n].residual, storage[n].k2);
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                captureRegisteredRHS(*registry, storage[n].registered, 2);
        }
    }

    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) applyStage(*fields[n], storage[n].q0, storage[n].k2, 0.5, state.dt);
        if (fields[n]) fields[n]->invalidateThermodynamicCache();
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                applyRegisteredStage(*registry, storage[n].registered,
                                     2, 0.5, state.dt, "RK4 stage 3");
        }
    }
    traceUpdatedConservative("Q after RK4 stage 3 update", fields, state);
    publish(fields);
    validate(fields, "RK4 stage 3");
    assembleRHS(fields, workspaces, state.time + 0.5 * state.dt);
    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) snapshotRHS(*fields[n], workspaces[n].residual, storage[n].k3);
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                captureRegisteredRHS(*registry, storage[n].registered, 3);
        }
    }

    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) applyStage(*fields[n], storage[n].q0, storage[n].k3, 1.0, state.dt);
        if (fields[n]) fields[n]->invalidateThermodynamicCache();
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                applyRegisteredStage(*registry, storage[n].registered,
                                     3, 1.0, state.dt, "RK4 stage 4");
        }
    }
    traceUpdatedConservative("Q after RK4 stage 4 update", fields, state);
    publish(fields);
    validate(fields, "RK4 stage 4");
    assembleRHS(fields, workspaces, state.time + state.dt);
    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) snapshotRHS(*fields[n], workspaces[n].residual, storage[n].k4);
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                captureRegisteredRHS(*registry, storage[n].registered, 4);
        }
    }

    for (size_t n = 0; n < fields.size(); ++n) {
        if (fields[n]) applyFinalRK4(*fields[n], storage[n], state.dt);
        if (fields[n]) fields[n]->invalidateThermodynamicCache();
        if (fields[n]) {
            if (auto* registry = registryFor(
                    *fields[n], fields, state, equationSystem))
                applyRegisteredFinal(*registry, storage[n].registered, state.dt);
        }
    }
    traceUpdatedConservative("Q after RK4 final update", fields, state);
    publish(fields);
    validate(fields, "RK4 final");
}

void advance(
    const std::vector<Field*>& fields,
    std::vector<PatchWorkspace>& workspaces,
        State::StateBundle& state,
        FDM::TimeScheme scheme,
        FDM::IEquationSystemCoupling* equationSystem,
        const AssembleRHS& assembleRHS,
        const Publish& publish,
        const Validate& validate) {
    if (scheme == FDM::TimeScheme::Euler) {
        stepEuler(fields, workspaces, state, equationSystem, assembleRHS, publish, validate);
    } else if (scheme == FDM::TimeScheme::SSPRK3) {
        stepSSPRK3(fields, workspaces, state, equationSystem, assembleRHS, publish, validate);
    } else {
    stepRK4(fields, workspaces, state, equationSystem, assembleRHS, publish, validate);
    }
}

} // namespace SF::SolverAlgorithm::DensityBasedTime
