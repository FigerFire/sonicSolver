# MixerVessel2DLevelSet

This case is a structured-grid sonicSolver adaptation of:

`/Volumes/OpenFOAM-v1912/tutorials/multiphase/multiphaseEulerFoam/mixerVessel2D`

## Mapping

- OpenFOAM `MRF1.cellZone rotor` -> `rotorZone` set in the self-contained `constant/mesh.sfm`.
- OpenFOAM `omega constant 10.472` -> `constant/MRFProperties`.
- OpenFOAM tank radius `R=1` with `scale=0.1` -> physical radius `0.10 m`.
- OpenFOAM MRF radius `ri=0.5*(0.05+0.07)` -> `rotorZone` radius `0.06 m`.
- OpenFOAM laminar setup -> `constant/turbulenceProperties` disables turbulence.
- OpenFOAM gravity `(0 0 0)` -> `constant/g` is zero.

## FDM Adaptation

The original OpenFOAM case is a multi-block circular O-grid with Eulerian
volume fractions for water, oil, mercury, and air. The current sonicSolver
Level Set module is a two-phase sign-field module, so this case uses:

- `water` as the default liquid phase.
- `airCore` as the gas phase set.
- The same center-block plus four outer arc blocks style as the project blockMeshDict examples.
- `frontBack` as the empty direction for a 2D slab; it is generated as a helper
  set from the `front` and `back` block faces.

`system/controlDict` advances one physical time step. The case therefore checks
the multi-patch EquationSystem lifecycle, Level-Set halo exchange, MRF source
coupling, and unified time driver instead of only checking mesh/config parsing.

Regenerate the mesh with the sonicSolver mesh generator, then add the
Level-Set/MRF helper sets:

```bash
./build/sonicSolver test/mixerVessel2DLevelSet_createMesh
python3 test/mixerVessel2DLevelSet/generateMixerVesselCase.py
```
