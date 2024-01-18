!===============================================================
!
! Copyright (C) 2015 Finite Difference Research Group
! This file is part of parsec, http://parsec.ices.utexas.edu/
!
! This subroutine writes processes the requested grid data and
! prints it to file. Currently uses cube format (see below).
! write_cube subroutine adapted in part from cube.f90 distributed
! in the Quantum ESPRESSO package, 2014.
!
! OFER: Probably does NOT work with symmetry on!
!
! author: O. Sinai (Weizmann), 2015
!
!---------------------------------------------------------------
subroutine export_grid_data(export_griddata_flag,parallel,pbc,grid,elec_st,pot)
  use constants
  use cluster_module
  use grid_module
  use electronic_struct_module
  use potential_module
  use pseudo_potential_module
  use parallel_data_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  !  flag that determines which data to export
  integer, intent(in) :: export_griddata_flag(MAX_EXPORT_OPTS)
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  potential related data
  type (potential), intent(inout) :: pot
  !
  ! Work variables:
  !
  ! data on grid
  real(dp) :: gridded_data(parallel%mydim)
  ! data collected to master processor, 1d
  real(dp) :: collected_data(grid%nwedge)
  ! data collected to master processor, 3d
  real(dp) :: grid3d_data(grid%n1,grid%n2,grid%n3)
  ! Coordinates of minimal x,y,z (origin)
  real(dp) :: origin_xyz(3)
  ! counters
  integer iflag,curflag,i,j,k,igrid,i3d,j3d,k3d
  ! output file name
  character (len=10) :: outfile_name
  ! Debug block work vars
  real(dp) :: dzabs, dzz
  integer jgrid

  !---------------------------------------------------------------

  ! Loop over data to print
  iflag = 0
  curflag = -1
  do while (curflag /= 0)
     iflag = iflag+1
     curflag = export_griddata_flag(iflag)
     if (curflag == 0) exit
     ! Identify current data to print
     gridded_data = zero
     select case (curflag)
     case (CHGDENS)
        call dcopy(parallel%mydim,elec_st%rho(:,1),1,gridded_data,1)
        outfile_name = 'charge_den'
     case (VEXT)
        call dcopy(parallel%mydim,pot%vion,1,gridded_data,1)
        outfile_name = 'v_external'
     case (VHAR)
        call dcopy(parallel%mydim,pot%vhart,1,gridded_data,1)
        outfile_name = 'v__hartree'
     case (VEXTHAR)
        gridded_data = pot%vion + pot%vhart
        outfile_name = 'v_elecstat'
     end select
   
!     ! OS_DEBUG block: create well-defined data.
!     if (export_griddata_flag == CHGDENS) then
!        do jgrid = 1,parallel%mydim
!           igrid = jgrid + parallel%irows(parallel%group_iam) - 1
!   
!           dzz = (grid%shift(3) + grid%kz(igrid))*grid%step(3) &
!                - 10                             ! z-position of valley minimum
!           dzabs = abs(dzz)
!   
!           ! Simply draw abs(z) position from minimum
!           gridded_data(jgrid) = dzabs
!        enddo
!     endif
!     ! End debug block
    
     ! Collect grid data over all processors.
     if (parallel%iammaster) collected_data = zero
     call collect_function(parallel,gridded_data)
     if (parallel%iammaster) call dcopy(grid%nwedge,parallel%ftmp,1,collected_data,1)
    
     if (parallel%iammaster) then
        ! Collect all data in a 3d grid
        grid3d_data = zero
        origin_xyz = zero
        i3d = 0
        do i = grid%nxmin+grid%norder, grid%nxmax-grid%norder
           i3d = i3d+1
           j3d = 0
           do j = grid%nymin+grid%norder, grid%nymax-grid%norder
              j3d = j3d+1
              k3d = 0
              do k = grid%nzmin+grid%norder, grid%nzmax-grid%norder
                 k3d = k3d+1
                 ! Index on general 1d grid
                 igrid = grid%indexw(i,j,k)
    
                 ! Save origin of data
                 if ((i3d == 1) .and. (j3d == 1) .and. (k3d == 1)) then
                    origin_xyz(1) = (grid%shift(1) + i)*grid%step(1)
                    origin_xyz(2) = (grid%shift(2) + j)*grid%step(2)
                    origin_xyz(3) = (grid%shift(3) + k)*grid%step(3)
                 endif
    
                 ! Take data only from points inside grid boundaries
                 if ((igrid == 0) .or. (igrid == grid%ndim+1)) cycle
    
!                 grid3d_data(i3d,j3d,k3d) = (grid%shift(3) + grid%kz(igrid))*grid%step(3)
                 grid3d_data(i3d,j3d,k3d) = collected_data(igrid)
              enddo
           enddo
        enddo

!        ! OS_DEBUG block
!        ! Build 3d grid
!        grid3d_data = zero
!        do i = 1,grid%n1
!           do j = 1,grid%n2
!              do k = 1,grid%n3
!                 if (export_griddata_flag == CHGDENS) then
!!                    grid3d_data(i,j,k) = dble(k)**2 - 30.0d0*dble(k)
!                    dzz = (grid%shift(3) + k)*grid%step(3) &
!                          - 10                             ! z-position of valley minimum
!                    grid3d_data(i,j,k) = abs(dzz)
!                 endif
!              enddo
!           enddo
!        enddo
!        ! End debug block
     
        ! Call printing subroutine
        write(7,*) 'Exporting ', outfile_name, ' grid data to file'
        call write_cube(pbc,grid,grid3d_data,origin_xyz,outfile_name)
        write(7,*) 'Finished exporting ', outfile_name

     endif

  enddo

end subroutine export_grid_data

!===============================================================
!
! This subroutine writes grid data to file. Currently exports in
! cube format with dummy data for atoms: a single H at origin.
! Adapted in part from cube.f90 distributed in the Quantum
! ESPRESSO package, 2014.
!
! author: O. Sinai (Weizmann), 2015
!
!---------------------------------------------------------------
subroutine write_cube(pbc,grid,grid3d_data,origin_xyz,outfile_name)
  use constants
  use cluster_module
  use grid_module
  use potential_module
  use pseudo_potential_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  data on 3d grid
  real(dp), intent(in) :: grid3d_data(grid%n1,grid%n2,grid%n3)
  ! Coordinates of minimal x,y,z (origin)
  real(dp), intent(in) :: origin_xyz(3)
  !  output file name
  character (len=*), intent(in) :: outfile_name
  !
  ! Work variables:
  !
  ! counters, file pointers, dummy variables
  integer i,j,k,outfile_unit,num_atoms,atomic_number
  ! Dummy variable
  real(dp) :: atomic_charge
  ! output file name with suffix
  character (len=len(outfile_name)+5) :: outfile_name_suffix

  !---------------------------------------------------------------
  !C     WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
  !C     TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
  !C     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
  !C
  !C     LINE   FORMAT      CONTENTS
  !C     ===============================================================
  !C      1     A           TITLE
  !C      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
  !C      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
  !C      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
  !C      #ATOMS LINES OF ATOM COORDINATES:
  !C      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
  !C      REST: 6E13.5      CUBE DATA
  !C
  !C     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS

  ! Dummy number
  num_atoms = 1

  ! Open file with given name
  outfile_name_suffix = outfile_name//'.cube'
  open(unit=outfile_unit,file=outfile_name_suffix,form='formatted')

  write(outfile_unit,*) 'Cubefile created from PARSEC calculation'
  write(outfile_unit,*) 'Volumetric '//outfile_name//' data'
  write(outfile_unit,'(I5,3F12.6)') num_atoms, (origin_xyz(i),i=1,3)
  if (pbc%per > 0) then
     write(outfile_unit,'(I5,3F12.6)') grid%n1, (pbc%latt_vec(i,1)/dble(grid%n1),i=1,3)
  else
     write(outfile_unit,'(I5,3F12.6)') grid%n1, grid%step(1), zero, zero
  endif
  if (pbc%per > 1) then
     write(outfile_unit,'(I5,3F12.6)') grid%n2, (pbc%latt_vec(i,2)/dble(grid%n2),i=1,3)
  else
     write(outfile_unit,'(I5,3F12.6)') grid%n2, zero, grid%step(2), zero
  endif
  if (pbc%per > 2) then
     write(outfile_unit,'(I5,3F12.6)') grid%n3, (pbc%latt_vec(i,3)/dble(grid%n3),i=1,3)
  else
     write(outfile_unit,'(I5,3F12.6)') grid%n3, zero, zero, grid%step(3)
  endif
 
  ! Write out (dummy) atom data
  do i = 1, num_atoms
     atomic_number = 1
     atomic_charge = dble(atomic_number)
     ! Atomic_charge could be alternatively set to valence charge.
     ! Positions are in cartesian coordinates and a.u.
     !
     ! So do this properly, should:
     ! * Find a way to work out the correct atom numbers.
     ! * Find a way to wrap atom coords back into the cell.
     write(outfile_unit,'(I5,5F12.6)') atomic_number, atomic_charge, zero, zero, zero 
  enddo

  ! Write out gridded data
  do i = 1,grid%n1
     do j = 1,grid%n2
        write(outfile_unit,'(6E13.5)') (grid3d_data(i,j,k),k = 1,grid%n3)
     enddo
  enddo

  ! Close file
  close(outfile_unit)

end subroutine write_cube
!===============================================================

!  !------ Original (only slightly modified) code from Amir's subroutine ----
!  ! position vectors and matrices
!  real(dp) :: u(3),r(3)
!  real(dp) :: tmpmat(3,3)
!
!  ! Load normalized lattice vectors
!  tmpmat = pbc%avec_norm
!
!  ! Go over given data and write to file
!  ! Go over all grid points on this processor
!  do jgrid = 1,parallel%mydim
!     ! Determine the point index in the general scheme of things as
!     ! the local index + the first general index stored on processor - 1
!     igrid = jgrid + parallel%irows(parallel%group_iam) - 1
!     u(1) = (grid%shift(1) + grid%kx(igrid))*grid%step(1)
!     u(2) = (grid%shift(2) + grid%ky(igrid))*grid%step(2)
!     u(3) = (grid%shift(3) + grid%kz(igrid))*grid%step(3)
!
!     call matvec3('N',tmpmat,u,r)
!
!     write(outfile_unit,25) r(1),r(2),r(3),gridded_data(jgrid)
!  enddo
!
!25 format(3(f15.9,','),e12.6)
