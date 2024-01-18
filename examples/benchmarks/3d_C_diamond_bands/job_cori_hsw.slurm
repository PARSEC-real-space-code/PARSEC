#!/bin/bash
#SBATCH -J C_diamond
#SBATCH --error=myjob.e%j
#SBATCH --output=myjob.o%j
##SBATCH --qos=debug
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --constraint=haswell
#SBATCH --reservation=vESW_2020_Thurs

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=1

module li
date

srun ~/PARSEC/parsec/src/parsec-cori_hsw_intel-ftn-ifort-19.0.3.199.mpi

date
