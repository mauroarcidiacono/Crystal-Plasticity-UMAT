# Crystal-Plasticity-UMAT
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18732831-blue)](https://doi.org/10.5281/zenodo.18732831) [![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Crystal plasticity UMAT for finite element analysis in Abaqus FEA.

The repository contains an Abaqus FEA UMAT and an emulator file written in Fortran90. It is possible to use the emulator to test the UMAT outside of Abaqus FEA. 

The original UMAT was published by Dr. Yonggang Huang and the original files can be found in [`huang_umat_archive`](./huang_umat_archive).

The UMAT has been modularised and further comments were added for clarity. Kinematic hardening was integrated, allowing to include up to three backstresses, and it was extended to model the hexagonal close-packed (HCP) crystal structure. Additionally, unit tests were added to ensure the proper functioning of core numerical routines.

A single-element benchmark is provided in `examples/one_element`, demonstrating numerical equivalence between the Abaqus FEA UMAT implementation and the standalone emulator under matched strain paths and constitutive settings.

> **Note**: the NONCUBICCRYSTAL subroutine is called when PROPS(9) is not zero. As of now, the subroutine contains a change of basis for a Ti-6Al-4V HCP crystal structure (the aspect ratio is defined in the subroutine). To use a different crystal structure, the change of basis matrix, normal vectors, direction vectors and number of slip systems per set must be modified. Otherwise, the flow control can be changed where NONCUBICCRYSTAL is called to use a different subroutine to define the crystal structure. 

## Runtime Interfaces

Two UMAT interfaces are provided:
- `umat_abaqus.for` – Abaqus-compatible version
- `umat_standalone.for` – Version used with the standalone emulator

Inside the Abaqus FEA environment, the file that should be selected as entry point is the `umat_abaqus.for`. The standalone version is to be used in the emulator.

## Documentation

- `UMAT-NOTES.md` – Detailed walkthrough of `core.for` algorithms and the code structure.
- `UMAT-VARIABLES.md` – Complete variable and parameter dictionary.
- `Mathematical_Notes_UMAT.pdf` – Full derivation of the constitutive equations, Newton–Raphson scheme, consistent tangent operator, and kinematic hardening formulation.

## How to Compile and Run the UMAT with the Emulator

Before compiling, update the input files:
- `cp_params.csv` contains the crystal plasticity material parameters. This file must be modified to match the material and model configuration you intend to use.
- `emulator.f90` defines the applied strain, number of state variables, number of properties, and tensor dimensions. These values must be consistent with your UMAT configuration and `cp_params.csv`.

The emulator does not perform automatic consistency checks. The user is responsible for ensuring that the parameters and array sizes are coherent.

Follow these steps to compile and run the emulator:

1. **Create the Object File for the Subroutine and the XIT subroutine**  
    Compile the subroutine into an object file:
    ```
    gfortran -c umat_standalone.for
    ```

    And compile the XIT subroutine:
    ```
    gfortran -c xit_emulator.for
    ```

2. **Create the Object File for the Emulator**  
    Compile the emulator file (e.g., `emulator.f90`) into an object file:
    ```
    gfortran -c emulator.f90
    ```

3. **Link the Object Files to Create the Executable**  
    Link the subroutine object file (`umat.o`) and the emulator object file (`emulator.o`) to create the executable:
    ```
    gfortran -o emulator umat_standalone.o xit_emulator.o emulator.o
    ```

4. **Run the Emulator**
    Execute the compiled emulator program:
    ```
    emulator
    ```

## Material Parameters

The example material parameter file (`cp_params.csv`) and the single-element validation case are based on crystal plasticity parameters reported in:

> Arcidiacono MF, Rahimi S.  
> *A surface integrity-informed crystal-plasticity based modelling of fatigue crack initiation in aerospace-grade Ti-6Al-4V.*  
> Materials Science and Technology. 2025;0(0).  
> https://doi.org/10.1177/02670836251365326

## Testing

Unit tests for key numerical subroutines are provided. To run the tests and validate the subroutines, refer to the `README.md` inside the `/tests` directory for step-by-step instructions.

## Cite this repository

If you use this code in your work, please cite it using the following Zenodo record:

**DOI:** [10.5281/zenodo.18732831](https://doi.org/10.5281/zenodo.18732831)

### BibTeX

```bibtex
@software{arcidiacono2026crystal,
  author    = {Arcidiacono, Mauro F. and Rahimi, Salaheddin},
  title     = {Crystal Plasticity UMAT for Abaqus FEA},
  version   = {v1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18732831},
  url       = {https://doi.org/10.5281/zenodo.18732831}
}
