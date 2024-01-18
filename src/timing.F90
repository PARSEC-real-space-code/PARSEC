!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! ** legacy text:
! **This is a machine-dependent subroutine that should be modified
! **when porting to new systems.
!
! For MPI we use MPI_WTIME and for serial we use cpu_time which 
! is probably system dependent
!---------------------------------------------------------------
subroutine mysecond(t)

  use constants
#ifdef MPI
  use mpi
#endif
  implicit none
  real(dp), intent(out) :: t 
#ifdef MPI
!$OMP MASTER
  t = MPI_WTIME()
!$OMP END MASTER
#else
!$OMP MASTER
  call cpu_time(t)
!$OMP END MASTER
#endif

end subroutine mysecond
!===============================================================
