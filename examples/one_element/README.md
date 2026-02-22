# One Element Validation Example

## Overview

This example demonstrates that the stress response obtained using the standalone emulator is numerically identical to the response obtained through Abaqus FEA.

The purpose of this test is to verify deterministic equivalence between:

- The Abaqus UMAT interface (`umat_abaqus.for`)
- The standalone constitutive emulator (`umat_standalone.for` and `emulator.f90`)

> **Important**: Due to the highly restrictive boundary conditions imposed on the lateral faces in this test, the resulting stress levels are higher than those observed in a classical uniaxial tensile test with free lateral contraction.

---

## Abaqus Model

The Abaqus simulation can be executed using the file `HCP-UMAT.inp` located in this directory.

### Model Configuration

- Geometry: 1 x 1 x 1 cube  
- Element type: single C3D8 element  
- Crystal structure: Hexagonal Close-Packed (HCP)  
- Material parameters representative of Ti-6Al-4V  

Material parameters were extracted from:

Arcidiacono MF, Rahimi S.  
*A surface integrity-informed crystal-plasticity based modelling of fatigue crack initiation in aerospace-grade Ti-6Al-4V.*  
Materials Science and Technology. 2025;0(0).  
doi:10.1177/02670836251365326

---

## Constitutive Settings

In this example:

- Small strain theory is used (`PROPS(146) = 0`)  
- `NLGEOM` is deactivated in Abaqus  
- The element is constrained on all faces except the loaded face  
- A uniaxial displacement is applied  
- Maximum applied strain: 0.005  
- The response remains purely elastic  

---

## Boundary Conditions and Emulator Consistency

The emulator reproduces the exact strain path imposed in the Abaqus model.

If alternative boundary conditions are implemented in the Abaqus input file, the strain evolution logic in `emulator.f90` must be modified so that the incremental strain history is consistent with the imposed constraints.

Maintaining consistent strain evolution is essential for deterministic equivalence between the standalone and Abaqus executions.

---

## Results

A comparison between stress component **S11** and strain component **E11** obtained from Abaqus FEA and the standalone emulator is shown in:

`stress_S11_strain_E11_comparison.png`

The results are numerically identical, confirming correct implementation and interface consistency between the two execution paths.

This example serves as a minimal verification benchmark for validating future modifications to the constitutive model or interface layers.