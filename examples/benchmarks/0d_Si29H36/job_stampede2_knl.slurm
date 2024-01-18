#!/bin/bash
#SBATCH -J Si29H36
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -t 00:30:00
#SBATCH -A Quantum-Theory-of-Ma

module load gcc
module list
pwd
date

ibrun ~/tt4/parsec/src/parsec-mpif90-gfortran-7.1.0.mpi
#ibrun ~/tt4/parsec/src/parsec-stampede2_knl_intel-ifort-18.0.2.mpi

date
