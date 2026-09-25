# Phase 1B numerics ownership plan

This classification was made before modifying numerical sources.  A file is
classified by the semantics it implements, rather than by its current folder.
`B` means CFD discretization ownership; it does not imply that every `B` file
can move in Phase 1B without creating a newly direct `models -> solver` or
`infrastructure -> solver` edge.

| Current file(s) | Classification | Target in this phase | Reason |
| --- | --- | --- | --- |
| `convection/SF_WENO.h` | B | `solver/discretization/convection/` | Acquires Field stencils, selects directions, handles IBM/ILW closure, and stores face flux. |
| `convection/SF_WENO3.{h,cpp}` | A | retain in `methods/numerics/convection/` | Fixed scalar stencil reconstruction and nonlinear weights only. |
| `convection/SF_WENO5.{h,cpp}` | A | retain in `methods/numerics/convection/` | Fixed scalar stencil reconstruction and nonlinear weights only. |
| `convection/SF_WENO7.{h,cpp}` | A | retain in `methods/numerics/convection/` | Fixed scalar stencil reconstruction and nonlinear weights only. |
| `convection/SF_TENO5.{h,cpp}` | A | retain in `methods/numerics/convection/` | Fixed scalar stencil reconstruction and nonlinear weights only. |
| `convection/SF_riemann.{h,cpp}` | A | retain in `methods/numerics/convection/` | Array-only Euler eigensystem and matrix kernels; no Field traversal. |
| `Flux/SF_face.h` | B | deferred | Field face geometry, IBM-cell interpretation, and flux storage.  Moving it would expose the existing `infrastructure/mpi/SF_haloExchange.cpp -> numerical flux` dependency as a new prohibited `infrastructure -> solver` include. |
| `Flux/SF_roe.{h,cpp}` | B | deferred with `SF_face.h` | CFD conservative face flux and Field writes. |
| `Flux/SF_lw.{h,cpp}` | B | deferred with `SF_face.h` | CFD face flux and IBM-cell checks. |
| `Flux/SF_lxF.{h,cpp}` | B | deferred with `SF_face.h` | CFD face flux; used by halo exchange. |
| `Flux/SF_sw.{h,cpp}` | B | deferred with `SF_face.h` | CFD face flux. |
| `Flux/SF_lowOrder.{h,cpp}` | B | deferred with `SF_face.h` | CFD face flux, face orientation, and IBM-cell checks. |
| `Flux/SF_rusanovEOS.{h,cpp}` | B | deferred with `SF_face.h` | Equation-set face flux and canonical-face geometry. |
| `structured/SF_iteration.h` | B | deferred | Field traversal plus active solver dimension; moving it makes numerous model and infrastructure callers direct solver consumers. |
| `structured/SF_stencil.h` | B | deferred with iteration | Field stencil acquisition and ghost-depth contract. |
| `structured/SF_canonicalFace.h` | B | deferred with iteration | Canonical face geometry and owner/neighbour sign semantics. |
| `structured/SF_faceGeometry.h` | B | deferred with canonical face | Field metric access. |
| `structured/SF_faceAssembly.h` | B | deferred with canonical face | Canonical-owner face-flux storage and residual assembly. |
| `structured/SF_residualAssembly.h` | B | deferred | Field residual assembly; its `applyDivergence` performs a time update. |
| `structured/SF_structured.h` | B | deferred umbrella | Existing solver `SF_structured.h` has the same basename; a move needs a later consumer migration, not a forwarding header. |
| `viscous/SF_newtonian.h` | A | `methods/math/SF_newtonian.h` | Stateless tensor stress and traction algebra. |
| `viscous/SF_viscous.h` | B | deferred | Central gradient/viscous-flux assembly.  Turbulence models currently call `velocityGradientAt`; moving it now would create new direct model-to-solver edges. |
| `scalar/SF_scalarTransport.h` | B + C | deferred, documented split | Boundary closure and RHS are B; `advanceForwardEuler` is C.  Current multiphase models use both through one public header, so splitting now would either add model-to-solver dependencies or create a compatibility facade. |
| `time/SF_Euler.{h,cpp}` | C | retain | Explicit time update. |
| `time/SF_RK4.{h,cpp}` | C | retain | RK stage handling and state commit. |
| `time/SF_time.h` | C | retain | Euler/RK/SSPRK lifecycle and stage utilities. |
| `state/SF_eulerState.h` | A | retain | Array-only Euler state/frame transforms. |
| `state/SF_equationState.h` | B | deferred | Field/equation-set state validation. |
| `state/SF_physicalState.h` | B | deferred | Field traversal and physical-state diagnostics. |
| `immersed/SF_constraintAlgebra.h` | A | retain | Row/weight algebra without Field or MPI. |
| `immersed/SF_interpolation.h` | A | retain | Generic row interpolation. |
| `immersed/SF_regularizedKernel.h` | A | retain | Scalar compact-support kernel. |
| `immersed/SF_spreading.h` | A | retain | Generic row contribution generation. |
| `immersed/SF_immersed.h` | A | retain | Pure immersed-kernel umbrella. |
| `immersed/SF_immersedOperators.h` | A, legacy include | deferred | Existing compatibility forwarding header; IBM source architecture is outside Phase 1B. |
| `SF_utility.h` | B + A | deferred | Scalar EOS helpers are A, but Field/EOS overloads make this a mixed Field-dependent header used by models and solver code. |
| `SF_numerics.h` | mixed umbrella | update only as required | It currently exposes both B and C APIs, so it cannot become a pure-math entrypoint in this phase. |

## Planned safe ownership moves

1. Move the Field-aware WENO assembly header into the existing solver
   convection directory and update its direct solver consumer.  The four
   WENO/TENO scalar kernels and the Riemann array kernels remain in `methods`.
2. Move only the Field-independent Newtonian stress/traction algebra to the
   existing `methods/math` module.
3. Do not hide the remaining B files behind forwarding headers.  Their
   non-solver callers are the concrete dependency blockers listed above and
   will remain explicit deferred work.

No time integration, state source, solver service, lifecycle, IBM algorithm,
or MPI ownership semantics will be changed.
