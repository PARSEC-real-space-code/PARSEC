# Pseudopotential algorithm for real-space electronic structure calculations (PARSEC)

## Current release

This version is PARSEC 1.4. It is the core engine of PARSEC that solves the Kohn–Sham equations.

## About PARSEC

PARSEC is a computer code that solves the Kohn-Sham equations by expressing electron wave functions directly in real space without the use of explicit basis sets. It uses norm-conserving pseudopotentials, such as Troullier-Martins and other varieties. It is designed for *ab initio* quantum mechanical calculations of the electronic structure of matter within density functional theory.

PARSEC is optimized for massively parallel computing environments but is also compatible with serial machines. A finite-difference approach is used for the calculation of spatial derivatives. Owing to the sparsity of the Hamiltonian matrix, the Kohn-Sham equations are solved by direct diagonalization using extremely efficient sparse-matrix eigensolvers. 

Some of its features are:
- Choice of boundary conditions: periodic (in all three directions) or confined
- Structural relaxation
- Simulated annealing
- Langevin molecular dynamics
- Polarizability calculations (confined-system boundary conditions only)
- Spin-orbit coupling
- Non-collinear magnetism

## How to compile
One can compile PARSEC under the folder [src](src) by:
```bash
make MACH=ubuntu_intel -j8
```
where `ubuntu_intel` is the name of the machine-dependent configuration file in [src/config](src/config). You may need to adjust some parameters to suit your programming environment. 

You can start with some examples (e.g., [examples/benchmarks/0d_Si28H36](examples/benchmarks/0d_Si28H36) for 0D systems, [examples/tests/bcc](examples/tests/bcc) for 3D systems, etc.). For instance, navigate to [examples/benchmarks/0d_Si28H36](examples/benchmarks/0d_Si28H36) and run the simulation using
```bash
mpirun -np 8 ../../../src/parsec-ubuntu_intel-ifx-2024.2.0.mpi
```

## Terms of usage

If you publish work using our code, please cite some of the following papers:

1.  James R. Chelikowsky, Norm Troullier, and Yousef Saad, *Finite-difference-pseudopotential method: electronic structure calculations without a basis*, [Physical Review Letters **72**, 1240](https://doi.org/10.1103/PhysRevLett.72.1240) (1994).
2.  James R. Chelikowsky, *The pseudopotential-density functional method applied to nanostructures*, [Journal of Physics D: Applied Physics **33**, R33](https://doi.org/10.1088/0022-3727/33/8/201) (2000).
3.  Leeor Kronik, Adi Makmal, Murilo L. Tiago, Manuel M. G. Alemany, Manish Jain, Xiangyang Huang, Yousef Saad, and James R. Chelikowsky, *PARSEC – the pseudopotential algorithm for real-space electronic structure calculations: recent advances and novel applications to nano-structures*, [Physica Status Solidi (b) **243**, 1063](https://doi.org/10.1002/pssb.200541463) (2006).
4.  Yousef Saad, James R. Chelikowsky, and Suzanne M. Shontz, *Numerical methods for electronic structure calculations of materials*, [SIAM Review **52**, 3](https://doi.org/10.1137/060651653) (2010). 
