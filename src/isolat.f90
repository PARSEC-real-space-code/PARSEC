!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine checks whether there are insulated atoms -
! i.e., atoms whose position prevents any bonding with any other
! atoms and may therefore lead to non-physical results.
!
! It must be placed INSIDE the atom movement loop, so that 
! the atom positions  will be rechecked after each molecular 
! dynamics or energy minimization step.
!
! WARNING!
! At present the only allowed grid setup involves one sphere 
! encompassing the entire cluster. Therefore the only check
! needed is that no atom is outside this boundary sphere. Should
! other grid constructions (e.g., an individual sphere around
! each atom) be implemented, the test(s) should be updated
! accordingly.
!
!---------------------------------------------------------------
subroutine isolat(clust,grid,nper,ierr)

  use constants
  use cluster_module 
  use grid_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(in) :: clust
  ! the grid
  type (grid_data), intent(in) :: grid
  ! Number of periodic directions, the first nper directions in the
  ! set (x,y,z) are assumed to be periodic.
  integer, intent(in) :: nper
  ! error flag, 290 < ierr < 301
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! index which is false if all is well, true if any isolated atom found
  logical isolated_found, iamisolated
  ! distance of coordinate from origin
  real(dp) :: dist 
  ! counter
  integer i,j,rmax
  ! temp r-vactor 
  real(dp) :: r_vec (3)

  ! cluster calculation outside-domain check
  logical, external :: outside_domain

  !---------------------------------------------------------------
  !
  ! PERIODIC SYSTEMS:
  ! For each atom, find its distance from the origin, if this
  ! distance is larger than rmax, flag a problem for the specific
  ! atom and for the whole cluster.
  !
  !---------------------------------------------------------------
  !
  ! CLUSTER SYSTEMS:
  ! Use the function outside_domain (also used in grid_partition)
  ! to see if any atoms are outside the computational domain.
  !
  !---------------------------------------------------------------
  !
  isolated_found = .false.
  rmax=grid%rmax
   do i=1, clust%atom_num
     r_vec = (/clust%xatm(i),clust%yatm(i),clust%zatm(i)/)

     ! use sphere-type criteria for periodic boundary conditions
     if(nper > 0) then
        r_vec(1:nper) = zero
        dist = sqrt (dot_product (r_vec,r_vec))
        iamisolated = (dist > rmax)
     else
        ! cluster calculation, check the point against the domain
        iamisolated = outside_domain(grid, r_vec)
     endif

     if(iamisolated) then
        if (.not. isolated_found) then
           isolated_found = .true.
           write(7,'(/,a,/)') 'ERROR: isolated atoms found:'
        endif
        write(7,20) i
     endif
   enddo

  if (isolated_found) then
     ! If problem found - report it and terminate.
     write(7,*) 'STOP in isolat.'
     ierr = 291
     return
  endif

10 format(1x,'Atom #',i4,2x,'is outside the boundary sphere!')
20 format(1x,'Atom #',i4,2x,'is outside the cluster domain!')

end subroutine isolat
!===============================================================
