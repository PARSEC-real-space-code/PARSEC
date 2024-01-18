!===================================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Define the distribution of grid points among processors so
!  that each one receives the same number of points (or the
!  difference in number of points held by all processors is at
!  most one).
!
!-------------------------------------------------------------------
subroutine grid_partition(grid,pbc,rsymm,parallel,ierr)

  use constants
  use grid_module
  use pbc_module
  use symmetry_module
  use parallel_data_module

  implicit none
  !
  !  Input/Output variables:
  !
  !  grid related data
  type (grid_data), intent(inout) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  !  symmetry operations in reduced (Abelian) group
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  error flag, 270 < ierr < 281
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  integer, dimension(1):: ierrvec
  !  cutoff radius used in the determination of grid points
  real(dp) :: rcut
  !  number of neighbors used for derivation (on one side)
  integer nord2
  !  actual size of hamiltonian, dimension of grid arrays
  integer ndim,ldn
  !  bounds of grid along each direction
  integer nxmin,nxmax,nymin,nymax,nzmin,nzmax
  !  variable for distance between points, maximum distance
  real(dp) :: dist,dmax
  !  temporary 1-d index of each point
  integer nef, neibs_num, nef2
  !  wrapping indexes
  integer iwrap, iwrap2, iwrap3
  !  counters
  integer ii,jj,i,j,k,j00,k00,ny1,ny2,nz1,nz2,iold
  !  variables for the irreducible wedge
  integer itrans, ir(3), itest, nwedge
  !  Cartesian coordinates of a grid point and its equivalent one
  real(dp) :: rr(3), rw(3), rstep
  !  arrays for irreducible wedge
  integer, allocatable :: ninv(:),rinv(:,:)

  !  number of procs available
  integer nnodes
  !  processor rank
  integer inode
  !  Indexing of distributed blocks
  integer irows(0:parallel%group_size)
  !  maximum number of FFT grid points (to be stored in pbc structure)
  integer maxdfft
  !  allocation check
  integer alcstat
  !
  ! grid limits - defines lower and upper limits of the grid for
  ! calculation of wrap around values
  integer :: grid_limits(3,6)
  !
  ! gwrap holds values of wrap around offset for the grid
  integer :: gwrap(3), gwrap2(3)
  ! idx1 is a temporary array that holds calculated grid_limits
  integer :: idx1(3,2), idx2(3,2)
  !
  !  External function:
  !
  integer, external :: Part
  logical, external :: outside_domain

  !-------------------------------------------------------------------

  ierr = 0
  nnodes = parallel%group_size
  nord2 = grid%norder
  !
  !  PART1: Setup the GRID
  !
  !  determine minimal number of mesh points (per each axis, per
  !  each side of origin) that would encompass a sphere of size rmax
  !  (confined system), or a box of side rmax (periodic system)
  !  inside which is a grid with a spacing of h. Note that for a
  !  periodic system rmax is the WHOLE side of the box, not half of
  !  it, so it's really comparable to the diameter of the sphere!

  nxmin = -grid%n1/2
  nxmax = grid%n1 + nxmin -1
  nymin = -grid%n2/2
  nymax = grid%n2 + nymin -1
  nzmin = -grid%n3/2
  nzmax = grid%n3 + nzmin -1

  nxmax = nxmax + nord2
  nxmin = nxmin - nord2
  nymax = nymax + nord2
  nymin = nymin - nord2
  nzmax = nzmax + nord2
  nzmin = nzmin - nord2

  call grid_set_index (nxmin,nxmax,nymin,nymax,nzmin,nzmax,grid)

  if (pbc%is_on .and. parallel%iammaster) then
     maxdfft = (pbc%n1 + grid%norder)*(pbc%n2 + grid%norder)* &
          (pbc%n3 + grid%norder)
     call pbc_set_maxdfft (maxdfft, pbc)
  endif

  nxmin = -grid%n1/2
  nxmax = grid%n1 + nxmin -1
  nymin = -grid%n2/2
  nymax = grid%n2 + nymin -1
  nzmin = -grid%n3/2
  nzmax = grid%n3 + nzmin -1
  !
  !  Assign spatial coordinates (x,y,z) to each array point (i,j,k).
  !  To each point, assign a 1d index (the counter is
  !  nef). Build the 3d->1d conversion arrays, kx, ky, kz. The
  !  actual position (atomic units, cartesian coordinates, relative
  !  to the origin) of each grid point is: xx=h*(kx+shift),
  !  yy=h*(ky+shift), zz=h*(kz+shift).
  !
  !  start by initializing the 3d -> 1d conversion index
  grid%indexg = 0
  !
  !  FOR CONFINED SYSTEM ONLY: must drop grid points outside a
  !  sphere of radius grid%rmax. For periodic system, use a cut-off
  !  radius (rcut) that certainly encompasses the periodic cell.
  !
  rcut = grid%rmax * grid%rmax
  nef = 0
  do i = nxmax, nxmin, -1
     rr(1) = (grid%shift(1) + i)*grid%step(1)
     do j = nymax, nymin, -1
        rr(2) = (grid%shift(2) + j)*grid%step(2)
        do k = nzmax, nzmin, -1
           rr(3) = (grid%shift(3) + k)*grid%step(3)

           ! Use sphere type criteria for periodic boundary conditions
           if(pbc%is_on) then
               rr(1:pbc%per) = zero
               dist = dot_product(rr,rr)
               if (dist > rcut) cycle
           else
               ! cluster calculation, check the point against the domain
               if(outside_domain(grid, rr)) cycle
           endif
           nef = nef + 1
           grid%indexg(i,j,k) = nef
        enddo
     enddo
  enddo
  !
  !  The total number of points inside the sphere is the size of the
  !  Hamiltonian of the system, ndim. Allocate arrays accordingly.
  !
  ndim = nef
  call grid_set_ndim (ndim, grid)

  !
  !
  ! Finding neighbor points and tagging them.
  if(grid%hartree_neibs_flag) then
   neibs_num=0
   grid%neibs_index=0
   do i = grid%nxmax, grid%nxmin, -1
      rr(1) = (grid%shift(1) + i)*grid%step(1)
      do j = grid%nymax, grid%nymin, -1
         rr(2) = (grid%shift(2) + j)*grid%step(2)
         do k = grid%nzmax, grid%nzmin, -1
            rr(3) = (grid%shift(3) + k)*grid%step(3)
            if(grid%indexg(i,j,k)==0) then ! point is not inside the domain 
             do jj=nord2, 1, -1 ! this could be replaced by just nord2 - AMIR
                               ! the only case were a neighbor is a 
                               ! neighbor with jj<nord2 and not a neighbor
                               ! with nord2 is problematic by itself.
              rstep=jj*grid%step(3)
              rw(1:2)=rr(1:2)
           
              rw(3)=rr(3)-rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if
           
              rw(3)=rr(3)+rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if

              rw(3)=rr(3)
              rw(1)=rr(1)-rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if

              rw(1)=rr(1)+rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if

              rw(1)=rr(1)
              rw(2)=rr(2)-rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if

              rw(2)=rr(2)+rstep
              if(.not. outside_domain(grid,rw)) then ! has neighbor inside
                 neibs_num=neibs_num+1
                 grid%neibs_index(i,j,k)=neibs_num
                 exit
              end if

             
             end do ! nord2 loop
            end if ! end of if point is outside domain
         end do ! end k loop
      end do ! end j loop
   end do  ! end i loop 

   grid%neibs_num=neibs_num

   call grid_set_ext_neibs(grid, parallel%procs_num)
  
   write(9,*) 'calculating grid neighbors'
   do i = grid%nxmax, grid%nxmin, -1
      do j = grid%nymax, grid%nymin, -1
         do k = grid%nzmax, grid%nzmin, -1
          if(grid%neibs_index(i,j,k) .ne. 0) then
            jj=grid%neibs_index(i,j,k)
            grid%neibs_fx(jj)=i
            grid%neibs_fy(jj)=j
            grid%neibs_fz(jj)=k
            write(9,*) jj,i,j,k
          end if
         end do
      end do
   end do

   write(7,*) 'finished building of external neighbors array'
   write(7,*) 'ndim is: ', ndim,' neibs_num is: ', neibs_num 
  endif ! if(grid%hartree_neibs_flag)

  !
  !  Define the set of grid points in the irreducible wedge, IW.
  !  The IW contains all grid points that can reconstruct the
  !  original grid by the application of symmetry operations, and it
  !  is such that no two points in the IW are equivalent to each
  !  other by the application of any symmetry operation.
  !
  !  Information about the wedge is contained in arrays indexw,
  !  rindex, rtrans:
  !  the n-th point in the IW has coordinates i,j,k so that
  !  indexw(i,j,k)=n. If indexw(i,j,k)=0, this point is not in the IW
  !
  !  if rindex(n)=m and rtrans(n)=i, the n-th point in original grid
  !  is equivalent to the m-th point in IW upon symmetry operation i.
  !
  !  dmax stores the maximum distance between a grid point and its
  !  equivalent in the wedge after a symmetry operation. Ideally,
  !  these two point should coincide, but numerical roundoff or the
  !  fractional translation symm%tnp (see subroutine symop) may interfere.
  dmax = zero
  nwedge = 0
  do i = nxmax, nxmin, -1
     rr(1) = (grid%shift(1) + i)*grid%step(1)
     do j = nymax, nymin, -1
        rr(2) = (grid%shift(2) + j)*grid%step(2)
        do k = nzmax, nzmin, -1
           rr(3) = (grid%shift(3) + k)*grid%step(3)
           if (grid%indexg(i,j,k) == 0) cycle
           do itrans = 1, rsymm%ntrans
              call symop(rsymm,itrans,rr,rw)
              do ii = 1, 3
                 ir(ii)=nint(rw(ii)/grid%step(ii)-grid%shift(ii))
                 rw(ii)=rw(ii) - (ir(ii) + grid%shift(ii))*grid%step(ii)
              enddo
              dist = sqrt(dot_product(rw,rw))
              dmax = max(dist,dmax)
              if (dist > 1.d-4) then
                 write(9,*) ' ERROR: grid point displaced from ' &
                      ,'image. Irreducible wedge may be wrong!'
                 write(9,*) ' Must ignore symmetry operations!'
                 write(9,*) rr
                 write(9,*) rw
                 write(9,*) ir
                 write(9,*) itrans
                 ierr = 271
                 return
              endif
              itest = grid%indexw(ir(1),ir(2),ir(3))
              if ( itest /= 0) exit
           enddo
           if (itest == 0) then
              nwedge = nwedge + 1
              grid%indexw(i,j,k) = nwedge
              itest = nwedge
              itrans = 1
           endif
           grid%rindex(grid%indexg(i,j,k)) = itest
           grid%rtrans(grid%indexg(i,j,k)) = itrans
        enddo
     enddo
  enddo
  grid%nwedge = nwedge
  ierrvec = ierr
  call pisum(ierrvec,1,parallel%groups_num,parallel%gmaster_comm)
  ierr = ierrvec(1)
  if (ierr > 0) then
     ierr = 271
     return
  endif

  call grid_set_wedge (grid,parallel%procs_num)

  ldn = nwedge/nnodes + nnodes
  parallel%ldn = ldn
  parallel%ndim = ndim
  parallel%nwedge = nwedge

  write(9,*)
  write(9,*) 'Setup messages:'
  write(9,*) '---------------'
  write(9,*) 'The Hamiltonian matrix size is', ndim
  write(9,*) ' reduced size is ',nwedge
  write(9,*) ' maximum distance between grid points and their images is'
  write(9,'(2x,g11.4,a)') dmax,'  [bohr]'

  write(9,*)
  write(9,15) 6*nord2
  if (parallel%iammaster) then
  write(7,*)
  write(7,*) 'Setup messages:'
  write(7,*) '---------------'
  write(7,*) 'The Hamiltonian matrix size is', ndim
  write(7,*) ' reduced size is ',nwedge
  write(7,*) ' maximum distance between grid points' &
       ,' and their images is'
  write(7,'(2x,g11.4,a)') dmax,'  [bohr]'

  write(7,*)
  write(7,15) 6*nord2
  endif
15 format(1x,'There are',1x,i2,1x, &
        'laplacian-related non-diagonal elements per row') 
  !
  !  PART2: Partitioning The Domain
  !
  !  At this point we are ready to do any kind of partitioning
  !  geometrical or based on the nwedge-array.
  !  Part(i,j,k,..) = processor number that the (i,j,k) resides on
  !
  irows(:) = 0
  !
  !  Find how many rows each processor has and how many equivalent
  !  points each grid point in IW has.
  !  Note: since grid%indexw already has information about which
  !  points to keep, we can run over (-nx,nx),(-ny,ny),(-nz,nz) for
  !  periodic and non-periodic boundary conditions. For equivalent
  !  points, search over all grid points where indexg in non-zero,
  !  since those ones define the full grid.
  !
  call grid_set_ist(grid,nnodes)

  allocate (ninv(grid%nwedge),stat=alcstat)
  call alccheck('ninv',grid%nwedge,alcstat) 
  ninv = 0
  allocate (rinv(rsymm%ntrans,grid%nwedge),stat=alcstat)
  call alccheck('rinv',rsymm%ntrans*grid%nwedge,alcstat)
  do i = nxmax, nxmin, -1
     do j = nymax, nymin, -1
        do k = nzmax, nzmin, -1
           if (grid%indexg(i,j,k) /= 0) then
              ii = grid%rindex(grid%indexg(i,j,k))
              ninv(ii) = ninv(ii) + 1
              rinv(ninv(ii),ii) = grid%indexg(i,j,k)
              grid%fx(grid%indexg(i,j,k)) = i
              grid%fy(grid%indexg(i,j,k)) = j
              grid%fz(grid%indexg(i,j,k)) = k
           endif
           if (grid%indexw(i,j,k) /= 0) then
              irows(Part(grid,nnodes,i,j,k,.true.)) &
                   =irows(Part(grid,nnodes,i,j,k,.true.)) + 1
           endif
        enddo
     enddo
  enddo
  if (maxval(ninv) > rsymm%ntrans .or. minval(ninv) < rsymm%ntrans) then
     if (parallel%iammaster) then
        write(7,*) 'ERROR: grid points in irreducible wedge ' &
             ,'have weights'
        write(7,*) 'different from the order of the reduced group'
        write(7,*) ' minimum and maximum values of weight =' &
             ,minval(ninv),maxval(ninv)
        write(7,*) ' order of the reduced group = ',rsymm%ntrans
     endif
     write(9,*) 'ERROR: grid points in irreducible wedge ' &
          ,'have weights'
     write(9,*) 'different from the order of the reduced group'
     write(9,*) ' minimum and maximum values of weight = ' &
          ,minval(ninv),maxval(ninv)
     write(9,*) ' order of the reduced group = ',rsymm%ntrans
     ierr = 272
  endif
  ierrvec = ierr
  call pisum(ierrvec,1,parallel%groups_num,parallel%gmaster_comm)
  ierr = ierrvec(1)
  if (ierr > 0) then
     ierr = 272
     return
  endif
  !
  !  Transform irows to show where the first row of each node
  !  begins. I.e., 3 nodes: first has 10 rows, second has 5 rows,
  !  third 6 rows irows(0) = 1, irows(1) = 11, irows(2) = 16,
  !  irows(3) = 22.
  !  Used for permuting the rows later.
  !
  irows(nnodes) = nwedge + 1
  do i = nnodes-1, 0, -1
     irows(i) = irows(i+1) - irows(i)
  enddo

  nef2=0
  !
  !  Permutation so that rows on the same processors appear in
  !  adjacent positions in the indexw numbering. Perform a
  !  "stiching" movement through space, in order to ensure that
  !  grid points with close (i,j,k) indices are also physically
  !  close.
  !
  j00 = -1
  k00 = -1
  do i = nxmin,nxmax
     j00 = -j00
     if (j00 == 1) then
        ny1 = nymin
        ny2 = nymax
     else
        ny1 = nymax
        ny2 = nymin
     endif
     do j = ny1, ny2, j00
        k00 = -k00
        if (k00 == 1) then
           nz1 = nzmin
           nz2 = nzmax
        else
           nz1 = nzmax
           nz2 = nzmin
        endif
        do k = nz1, nz2, k00
           if (grid%indexw(i,j,k) /= 0) then
              iold = grid%indexw(i,j,k)
              inode = Part(grid,nnodes,i,j,k,.true.)
              nef           = irows( inode )
              grid%indexw(i,j,k) = nef
              do ii = 1, rsymm%ntrans
                 grid%rindex(rinv(ii,iold)) = nef
              enddo
              irows( inode )= irows( inode ) + 1
              grid%kx(nef)  = i
              grid%ky(nef)  = j
              grid%kz(nef)  = k
           endif
        enddo
     enddo
  enddo
  deallocate(ninv, rinv)
  !
  !  Restoring irows so that they point at the first row of each
  !  processor
  !
  irows(nnodes) = grid%nwedge + 1
  do i = nnodes -1, 1,-1
     irows(i) = irows(i-1)
  enddo
  irows(0) = 1
  !
  !  Transfer irows to parallel structure
  !
  parallel%irows = irows
  !
  !  if in a periodic system, points must be "wrapped around" the
  !  boundary so that points on opposite sides of the box are
  !  actually recognized as neighbors. For that, we go over
  !  -nx-nord2 to -nx-1 and nx to nx+nord2-1 (and so on in the y
  !  and z directions), thus finding all the neighbors of points
  !  that have neighbors that are "wrapped" on the other side of the
  !  box. For each point, we go over 6 combinations in all,
  !  coresponding to the 6 faces of the cube. Every time i00, j00
  !  and k00 are flipped in sign accordingly so that all the 6
  !  combinations are covered. In any one of these combinations, the
  !  elements are mapped to the corresponding element on the
  !  other side of the boundary (i.e., shifted by 2*nx and so on).
  !

  grid_limits(1,1)=nxmin-nord2
  grid_limits(1,2)=nxmin-1
  grid_limits(1,3)=nxmax+1
  grid_limits(1,4)=nxmax+nord2
  grid_limits(1,5)=nxmin
  grid_limits(1,6)=nxmax

  grid_limits(2,1)=nymin-nord2
  grid_limits(2,2)=nymin-1
  grid_limits(2,3)=nymax+1
  grid_limits(2,4)=nymax+nord2
  grid_limits(2,5)=nymin
  grid_limits(2,6)=nymax

  grid_limits(3,1)=nzmin-nord2
  grid_limits(3,2)=nzmin-1
  grid_limits(3,3)=nzmax+1
  grid_limits(3,4)=nzmax+nord2
  grid_limits(3,5)=nzmin
  grid_limits(3,6)=nzmax

  gwrap(1)=grid%n1
  gwrap(2)=grid%n2
  gwrap(3)=grid%n3

  if (pbc%per >= 1) then
  !  ---- -x direction
     do i = grid_limits(1,1),grid_limits(1,2)
        do j = grid_limits(2,5),grid_limits(2,6)
           do k = grid_limits(3,5),grid_limits(3,6)
              grid%indexg(i,j,k) = grid%indexg(i+gwrap(1),j,k)
           enddo
        enddo
     enddo
  !  ---- +x direction
     do i = grid_limits(1,3),grid_limits(1,4)
        do j = grid_limits(2,5),grid_limits(2,6)
           do k = grid_limits(3,5),grid_limits(3,6)
              grid%indexg(i,j,k) = grid%indexg(i-gwrap(1),j,k)
           enddo
        enddo
     enddo
  endif
  if (pbc%per >= 2) then
  !  ---- -y direction
     do i = grid_limits(1,5),grid_limits(1,6)
        do j = grid_limits(2,1),grid_limits(2,2)
           do k = grid_limits(3,5),grid_limits(3,6)
              grid%indexg(i,j,k) = grid%indexg(i,j+gwrap(2),k)
           enddo
        enddo
     enddo
  !  ---- +y direction
     do i = grid_limits(1,5),grid_limits(1,6)
        do j = grid_limits(2,3),grid_limits(2,4)
           do k = grid_limits(3,5),grid_limits(3,6)
              grid%indexg(i,j,k) = grid%indexg(i,j-gwrap(2),k)
           enddo
        enddo
     enddo
  endif
  if (pbc%per >= 3) then
  !  ---- -z direction
     do i = grid_limits(1,5),grid_limits(1,6)
        do j = grid_limits(2,5),grid_limits(2,6)
           do k = grid_limits(3,1),grid_limits(3,2)
              grid%indexg(i,j,k) = grid%indexg(i,j,k+gwrap(3))
           enddo
        enddo
     enddo
  !  ---- +z direction
     do i = grid_limits(1,5),grid_limits(1,6)
        do j = grid_limits(2,5),grid_limits(2,6)
           do k = grid_limits(3,3),grid_limits(3,4)
              grid%indexg(i,j,k) = grid%indexg(i,j,k-gwrap(3))
           enddo
        enddo
     enddo
  endif

  do ii=1,3
    if(grid%lap_dir(ii) /= 0) then
      do jj=1,3
       if(grid%lap_neig(jj,ii)>0) then
         idx1(jj,1)=1
         idx1(jj,2)=3
       else
         if(grid%lap_neig(jj,ii)<0) then
           idx1(jj,1)=3
           idx1(jj,2)=1
         else
           idx1(jj,1)=5
           idx1(jj,2)=5
         endif
       endif
      enddo

     do i = grid_limits(1,idx1(1,1)) &
          , grid_limits(1,idx1(1,1)+1)
      do j = grid_limits(2,idx1(2,1)) &
           , grid_limits(2,idx1(2,1)+1)
        do k = grid_limits(3,idx1(3,1)) &
             , grid_limits(3,idx1(3,1)+1)

           grid%indexg(i,j,k)= &
            grid%indexg(i+gwrap(1)*grid%lap_neig(1,ii) &
                       ,j+gwrap(2)*grid%lap_neig(2,ii) &
                       ,k+gwrap(3)*grid%lap_neig(3,ii))
       enddo
      enddo
     enddo


    do i = grid_limits(1,idx1(1,2)) &
          , grid_limits(1,idx1(1,2)+1)
      do j = grid_limits(2,idx1(2,2)) &
           , grid_limits(2,idx1(2,2)+1)
        do k = grid_limits(3,idx1(3,2)) &
             , grid_limits(3,idx1(3,2)+1)
           grid%indexg(i,j,k)= &
              grid%indexg(i-gwrap(1)*grid%lap_neig(1,ii) &
                         ,j-gwrap(2)*grid%lap_neig(2,ii) &
                         ,k-gwrap(3)*grid%lap_neig(3,ii))
       enddo
      enddo
     enddo

     !the following if treats the case where all coordinates
     !of the neighbor are different than zero. In such a case
     !one has to treat not just the combined edge of all 3 indices 
     !but also the situation were one of the indices is inside the 
     !permitted grid and the other 2 are outside. We assume that all
     !neighbor coordinates are +/- 1. 

     if(sum(abs(grid%lap_neig(1:3,ii)))>2) then
       do jj=1,3
         idx2=idx1
         gwrap2=gwrap
         idx2(jj,1)=5
         idx2(jj,2)=5
         gwrap2(jj)=0

         do i = grid_limits(1,idx2(1,1)) &
          , grid_limits(1,idx2(1,1)+1)
          do j = grid_limits(2,idx2(2,1)) &
            , grid_limits(2,idx2(2,1)+1)
            do k = grid_limits(3,idx2(3,1)) &
               , grid_limits(3,idx2(3,1)+1)

            grid%indexg(i,j,k)= &
             grid%indexg(i+gwrap2(1)*grid%lap_neig(1,ii) &
                        ,j+gwrap2(2)*grid%lap_neig(2,ii) &
                        ,k+gwrap2(3)*grid%lap_neig(3,ii))
            enddo
          enddo
         enddo

         do i = grid_limits(1,idx2(1,2)) &
          , grid_limits(1,idx2(1,2)+1)
          do j = grid_limits(2,idx2(2,2)) &
            , grid_limits(2,idx2(2,2)+1)
            do k = grid_limits(3,idx2(3,2)) &
               , grid_limits(3,idx2(3,2)+1)

            grid%indexg(i,j,k)= &
             grid%indexg(i-gwrap2(1)*grid%lap_neig(1,ii) &
                        ,j-gwrap2(2)*grid%lap_neig(2,ii) &
                        ,k-gwrap2(3)*grid%lap_neig(3,ii))
            enddo
          enddo
         enddo

      enddo
     endif
   endif
  enddo

 
 
!
! Set the indexes for grid points outside the domain to be nwedge+1 (it was
! zero before). This is relevant along non-periodic directions only.
!
  do k = nzmin - nord2,nzmax + nord2
     do j = nymin - nord2,nymax + nord2
        do i = nxmin - nord2,nxmax + nord2
           if (grid%indexw(i,j,k) == 0) grid%indexw(i,j,k) = grid%nwedge + 1
           if (grid%indexg(i,j,k) == 0) grid%indexg(i,j,k) = grid%ndim + 1
        enddo
     enddo
  enddo

  grid%rindex(grid%ndim+1) = grid%nwedge + 1
  grid%rtrans(grid%ndim+1) = 1

  !  Finally, define grid shift in units of lattice vector
  if (pbc%per >= 1) then
     pbc%mx = minval(grid%fx)
     pbc%shift(1) = (grid%shift(1) + real(pbc%mx,dp)) * &
          grid%step(1) * twopi / pbc%box_size(1)
  else
     pbc%mx = 0
     pbc%shift(1) = zero
  endif
  if (pbc%per >= 2) then
     pbc%my = minval(grid%fy)
     pbc%shift(2) = (grid%shift(2) + real(pbc%my,dp)) * &
          grid%step(2) * twopi / pbc%box_size(2)
  else
     pbc%my = 0
     pbc%shift(2) = zero
  endif
  if (pbc%per >= 3) then
     pbc%mz = minval(grid%fz)
     pbc%shift(3) = (grid%shift(3) + real(pbc%mz,dp)) * &
          grid%step(3) * twopi / pbc%box_size(3)
  else
     pbc%mz = 0
     pbc%shift(3) = zero
  endif
  if (parallel%iammaster .and. pbc%is_on) then
     write(7,*) 'pbc shift indices = ',pbc%mx,pbc%my,pbc%mz
     write(7,'(a,3f10.4)')  &
          ' pbc shift [latt. vect. units] = ',pbc%shift/twopi
  endif

end subroutine grid_partition
!===============================================================
!
!  Simple partitioning routine.
!  Split the nwedge-long vector into nnodes equal parts.
!
!---------------------------------------------------------------
integer function Part(grid,nnodes,i,j,k,lprint)

  use constants
  use grid_module
  implicit none
  !
  !  Input/Output variables:
  !
  type (grid_data), intent(inout) :: grid
  integer, intent(in) :: nnodes
  integer, intent(in) :: i,j,k
  logical, intent(in) :: lprint

  integer nef, ii
!---------------------------------------------------------------

  !  Give the processor number for the current i,j,k
  !  For this routine partition according to rows:
  !  (ijk) -> Indexw(ijk)
  !
  nef = grid%indexw(i,j,k)
  do ii = 0, nnodes-1
     if (nef >= grid%ist(ii) .and. nef < grid%ist(ii+1)) then
        Part = ii
        return
     endif
  enddo

  Part = -1
  if (lprint) write(7,*) 'The row nef:',nef &
       ,'is outside bounds nwedge:',grid%nwedge

end function Part
!===============================================================
!
!  Determines if the point rr is outside the computational
!  domain.  The grid structure is needed since it contains
!  the type of domain shape and the necessary parameters.
!
!---------------------------------------------------------------
logical function outside_domain(grid,rr)

  use constants
  use grid_module
  implicit none
  !
  !  Input/Output variables:
  !
  type (grid_data), intent(in) :: grid
  real(dp), intent(in) :: rr(3)
  !
  !  Work variables:
  !
  real(dp) :: rad2, rr_rad2, coord, polar1, polar2
  !---------------------------------------------------------------

  ! generic return -> point rr is outside the domain
  outside_domain = .true.
  
  select case(grid%domain_shape)
  case (0)
     ! Sphere
     ! points inside the sphere obey 
     ! x**2 + y**2 + z**2 <= r**2

     ! grid%d_shape_param(1) is the radius
     rad2 = grid%d_shape_param(1)**2

     rr_rad2 = dot_product(rr, rr)

     outside_domain = (rr_rad2 > rad2)
  case (1)
     ! Ellipsoid
     ! points inside the sphere obey 
     ! (x/r_x)**2 + (y/r_y)**2 + (z/r_z)**2 <= 1

     outside_domain = ((rr(1) / grid%d_shape_param(1))**2 + &
        (rr(2) / grid%d_shape_param(2))**2 + &
        (rr(3) / grid%d_shape_param(3))**2 > 1.0d0)
  case (2)
     ! Cylinder
     ! points inside the cylinder obey
     ! x**2 + y**2 <= r**2 AND |z| <= length
     ! for a cylinder oriented along the z-axis

     ! Get the length coordinate and the polar coordinates
     select case(grid%i_shape_param(1))
     case (1)
        ! oriented along x-axis
        coord = rr(1)
        polar1 = rr(2)
        polar2 = rr(3)
     case (2)
        ! oriented along y-axis
        coord = rr(2)
        polar1 = rr(1)
        polar2 = rr(3)
     case (3)
        ! oriented along z-axis
        coord = rr(3)
        polar1 = rr(2)
        polar2 = rr(1)
     end select

     outside_domain = ((polar1**2 + polar2**2 > grid%d_shape_param(1)**2) &
        .or. coord**2 > (grid%d_shape_param(2)/2.0d0)**2)
  case (3)
     ! Box
     ! If the coordinate is outside the length/2 (in magnitude)...

     outside_domain = (rr(1)**2 > (grid%d_shape_param(1)/2.0d0)**2 .or. &
        rr(2)**2 > (grid%d_shape_param(2)/2.0d0)**2 .or. &
        rr(3)**2 > (grid%d_shape_param(3)/2.0d0)**2)
  end select

end function outside_domain
!===============================================================
