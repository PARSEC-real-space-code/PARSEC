# Pseudopotential Algorithm of Real-Sapce Electronic Structure Calculation (PARSEC)

## Current Release

The latest release is PARSEC 1.4. It is the core engine of PARSEC that solve for the solutions of Kohn–Sham equations.

## About PARSEC

PARSEC is a computer code that solves the Kohn-Sham equations by expressing electron wave-functions directly in real space, without the use of explicit basis sets. It uses norm-conserving pseudopotentials (Troullier-Martins and other varieties). It is designed for ab initio quantum-mechanical calculations of the electronic structure of matter, within density-functional theory.

PARSEC is optimized for massively parallel computing environment, but it is also compatible with serial machines. A finite-difference approach is used for the calculation of spatial derivatives. Owing to the sparsity of the Hamiltonian matrix, the Kohn-Sham equations are solved by direct diagonalization, with the use of extremely efficient sparse-matrix eigensolvers. 

Some of its features are:
- Choice of boundary conditions: periodic (on all three directions), or confined.
- Structural relaxation.
- Simulated annealing.
- Langevin molecular dynamics.
- Polarizability calculations (confined-system boundary conditions only).
- Spin-orbit coupling.
- Non-collinear magnetism.

## How to Compile
One can compile PARSEC under the folder src/ by 

```bash
make MACH=ccm_intel -j8
```

where *ccm_intel* is the name of the machine-dependent configuration file in src/config. One may need to change some parameters to suit one's programming environment. 

One may want to start with some examples (e.g. benchmarks/Si28H36 for 0D systems, tests/bcc for 3D systems, etc). For instance, one can go to benchmarks/Si28H36, and run the simulation by

```bash
mpiexec -n 8 ../../src/parsec-mpif90-ifort-17.0.4.mpi
```

## Terms of Usage

If you publish work using our code, please cite some of the following articles:

J. R. Chelikowsky, N. Troullier, and Y. Saad, Finite-difference pseudo potential method: Electronic structure calculations without a basis, Phys. Rev. Lett. 72, 1240 (1994).

J. R. Chelikowsky, The pseudopotential-density functional method applied to nanostructures, J. of Phys. D 33, R33 (2000).

L. Kronik, A. Makmal, M.L. Tiago, M.M.G. Alemany, M. Jain, X. Huang, Y. Saad, and J.R. Chelikowsky, PARSEC-the pseudopotential algorithm for real-space electronic structure calculations: Recent Advances and novel applications to nano-structures, physica status solidi (b) 243, 1063 (2006).

Y. Saad, J.R. Chelikowsky and S.M. Shontz, Numerical methods for electronic structure calculations of materials, SIAM Rev. 52, 3 (2010). 
