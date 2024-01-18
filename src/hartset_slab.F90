!===============================================================
!
!  Boundary conditions and builder of charge distribution for 
!  2-dimensional (slab) systems.
!
!  Generaly based on the calculation of the potential outside a sheet
!  of periodic charges, as described by Lenard-Jones & Dent, 
!  "Cohesion at a crystal surface", Trans. Faraday Soc, 1927
!
! author: Ayelet Mor (2007)
! 
!---------------------------------------------------------------
subroutine hartset_slab(grid,pbc,parallel,norder,lpole, &
     coeke,rho,brho)

  use constants
  use grid_module
  use pbc_module
  use parallel_data_module
  implicit none
  
  !
  !  Input/Output variables:
  !  -----------------------
  
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  number of neighbors used for derivation (on one side)
  integer, intent(in) :: norder
  !  order of multipole expansion
  integer, intent(in) :: lpole
  !  coefficients for taking a second derivative
  real(dp), intent(in) :: coeke(-norder:norder,3+grid%lap_dir_num)
  !  charge density on the 3-d grid
  real(dp), intent(in) :: rho(parallel%ldn)
  !  effective charge density = original charge density augmented by
  !  appropriate boundary terms, to be used as the source term in hpotcg
  real(dp), intent(out) :: brho(parallel%ldn)
  
  !
  !  Work variables:
  !  ---------------

  !  offset for each processor
  integer ioffset
  !  deta = det(latt_vec)
  real(dp) :: deta
  !  twopia = 2 * pi / deta
  real(dp) :: twopia
  ! nzmax = grid%nzmax - grid%norder
  integer nzmax
  ! nzmin = grid%nzmin + grid%norder
  integer nzmin
  !  nz = nzmax - nzmin + 1
  integer nz
  ! nzneighs = grid%nzmax - grid%nzmin + 1
  integer nzneighs
  !  lap_neibs = 2*( 3 + grid%lap_dir_num)
  integer lap_neibs

  !  Variables for each grid point (labeled by j)
  real(dp) :: &    
       rrho, &             ! = rho(j)
       rr(2)               ! vector of (u(j), v(j))

  ! index variables
  integer :: &
       izz, &              ! = z(j)/grid%step(3)
       izdiff              ! = abs(z(j) - Z(neighbor(j))/grid%step(3)
 
  !  Vectors of the fourier expansion
  !  klm(l,m) = l*b1 + m*b2  where b1,b2 are the reciprocal lattice vectors
  !  kklm(l,m) = |klm(l,m)|
  !  so that klm(l,m) = -klm(-l,-m); klm(-l,m) = -klm(l,-m)
  !
  !  in practice, since klm is used only to be multiplied by r,
  !  klm(l,m) is set to be (l,m) and r is (u,v) instead of the 
  !  correct form where klm = (l,m)*(b1,b2), r = (u,v)*(a1,a2)
  !  since b_i*a_j = 2pi*delta_ij

  real(dp) :: klm(-lpole:lpole, 0:lpole, 2)
  real(dp) :: kklm(-lpole:lpole, 0:lpole)

  !  Temporary storage of the sine and cosine functions
  real(dp) :: clm(-lpole:lpole, 0:lpole)
  real(dp) :: slm(-lpole:lpole, 0:lpole)

  !  Multipoles (one for the sine part and one for the cosine part)
  !  qslm(izz) = sum(rho(j) * slm(j,l,m)) on all points j with height zz
  !  qclm(izz) = sum(rho(j) * clm(j,l,m)) on all points j with height zz

  real(dp), allocatable :: qslm(:,:,:),qclm(:,:,:)

  !  Z axis exponents
  !  tlm(l,m,zdiff) = exp(-kklm(l,m)*zdiff)

  real(dp), allocatable :: tlm(:,:,:)

  !  Contribution to the boundary potential (the contribution should
  !  be real, since the summation is over real functions).
  real(dp) :: cpole

  !  Average charge density for each xy plane in the slab

  real(dp), allocatable :: sigma(:), sig_counter(:)

  !  other counters
  integer ii,ll,mm,zindex,jrow,jshell

  !  temporarly variables
  real(dp) :: tmp(3),rs(3)
  real(dp) :: tmp2
  real(dp) :: vect1(2), vect2(2)
  character(len=4) :: idstring
  integer :: memcount, alcstat

  !  constants
  real(dp), parameter :: pi8 = eight*pi
  real(dp), parameter :: pi2 = two*pi
  !---------------------------------------------------------------


  !
  ! initialize variables
  !
  write(idstring,'(I4.4)') parallel%iam
  !open(unit=15, file='cpole.dat.'//idstring, form='formatted')
  !open(unit=16, file='qclm.dat.'//idstring, form='formatted')

  lap_neibs = 2*( 3 + grid%lap_dir_num )
  deta = abs(pbc%latt_vec(1,1)*pbc%latt_vec(2,2) - &
         pbc%latt_vec(2,1)*pbc%latt_vec(1,2))
  twopia = pi2 / deta

  nzmax = grid%nzmax - grid%norder 
  nzmin = grid%nzmin + grid%norder + 1
  nz = nzmax - nzmin + 1
  nzneighs = grid%nzmax - grid%nzmin + 1
  
  ioffset = parallel%irows(parallel%group_iam) - 1
  
  !  Allocate arrays that need to be allocated
  memcount =(lpole-(-lpole)+1) * (lpole-0+1) * &
    ((grid%nzmax-grid%norder)-(grid%nzmin+grid%norder)+1)
  allocate(qslm(-lpole:lpole, 0:lpole, & 
           grid%nzmin+grid%norder:grid%nzmax-grid%norder), &
           stat=alcstat)
  call alccheck('qslm',memcount,alcstat)

  memcount =(lpole-(-lpole)+1) * (lpole-0+1) * &
    ((grid%nzmax-grid%norder)-(grid%nzmin+grid%norder)+1)
  allocate(qclm(-lpole:lpole, 0:lpole, & 
           grid%nzmin+grid%norder:grid%nzmax-grid%norder), &
           stat=alcstat)
  call alccheck('qclm',memcount,alcstat)

  memcount =(lpole-(-lpole)+1) * (lpole-0+1) * &
    (grid%nzmax-grid%nzmin+2)
  allocate(tlm(-lpole:lpole, 0:lpole, 0:(grid%nzmax-grid%nzmin+1)), &
           stat=alcstat)
  call alccheck('tlm',memcount,alcstat)

  memcount = ((grid%nzmax-grid%norder)-(grid%nzmin+grid%norder)+1)
  allocate(sigma(grid%nzmin+grid%norder:grid%nzmax-grid%norder), &
           stat=alcstat)
  call alccheck('sigma',memcount,alcstat)

  memcount = ((grid%nzmax-grid%norder)-(grid%nzmin+grid%norder)+1)
  allocate(sig_counter(grid%nzmin+grid%norder:grid%nzmax-grid%norder), &
           stat=alcstat)
  call alccheck('sig_counter', memcount, alcstat)
  
  !
  !  Zero out all the arrays
  !

  brho(:) = zero
  klm(:,:,:) = zzero
  kklm(:,:) = zzero
  clm(:,:) = zzero
  slm(:,:) = zzero
  qclm(:,:,:) = zzero
  qslm(:,:,:) = zzero
  tlm(:,:,:) = zzero
  sigma(:) = zero
  sig_counter(:) = zero

  !
  ! Define the klm and kklm tables.
  !
  
  vect1 = (/1,0/) * pi2/(sqrt(dot_product(pbc%latt_vec(:,1), pbc%latt_vec(:,1))))
  vect2 = (/0,1/) * pi2/(sqrt(dot_product(pbc%latt_vec(:,2), pbc%latt_vec(:,2))))

  do ll = -lpole, lpole
     do mm = 1, lpole
        klm(ll,mm,:) = ll*vect1 + mm*vect2 
        tmp = ll*pbc%bvec(:,1) + mm*pbc%bvec(:,2)
        kklm(ll,mm) = sqrt(dot_product(tmp,tmp))
     enddo
  enddo

  do ll = 1, lpole
     tmp = ll*pbc%bvec(:,1)
     tmp2 = sqrt(dot_product(tmp,tmp))
     klm(ll,0,:) = ll*vect1
     kklm(ll,0) = tmp2
  enddo


  !
  !  Build the qclm and qslm tables.
  !  add in the contribution of every grid point we have data for.
  !  build the tables for each izz height.

  do ii = parallel%irows(parallel%group_iam), &
          parallel%irows(parallel%group_iam+1)-1

     rrho = rho(ii-ioffset)

     ! Convert from 1d indexing to 3d indexing (u,v,z)
     rr(1) = (grid%shift(1) + grid%kx(ii)) * grid%step(1)
     rr(2) = (grid%shift(2) + grid%ky(ii)) * grid%step(2)
     izz = grid%kz(ii)

     ! Initiate the clm and slm tables for the given grid point
     call cos_lm(rr,lpole,klm,clm(-lpole:lpole, 0:lpole))
     call sin_lm(rr,lpole,klm,slm(-lpole:lpole, 0:lpole))


     ! Fill in the qslm & qclm tables, without the case of (ll,mm)=(0,0)
     do ll = -lpole, lpole
        do mm = 1, lpole
                qclm(ll,mm,izz) = qclm(ll,mm,izz) + rrho * clm(ll,mm)
                qslm(ll,mm,izz) = qslm(ll,mm,izz) + rrho * slm(ll,mm)
        enddo
     enddo

     do ll = 1, lpole
        qclm(ll,0,izz) = qclm(ll,0,izz) + rrho * clm(ll,0)
        qslm(ll,0,izz) = qslm(ll,0,izz) + rrho * slm(ll,0)
     enddo

  enddo   ! ii


  ! 
  ! Global sum for colecting the qclms and qslms from each processor.
  !  
  call psum(qclm,(2*lpole+1)*(lpole+1)*(nz), &
  parallel%group_size,parallel%group_comm)
  call psum(qslm,(2*lpole+1)*(lpole+1)*(nz), &
  parallel%group_size,parallel%group_comm)

  !
  ! Multiply each qslm, qclm by twopia/kklm(ll,mm)
  ! 

  do ll = -lpole, lpole
     do mm = 1, lpole
        tmp2 = one / kklm(ll,mm)
        qclm(ll,mm,:) = qclm(ll,mm,:) * tmp2
        qslm(ll,mm,:) = qslm(ll,mm,:) * tmp2
     enddo
  enddo

  do ll = 1, lpole
     tmp2 = one / kklm(ll,0)
     qclm(ll,0,:) = qclm(ll,0,:) * tmp2
     qslm(ll,0,:) = qslm(ll,0,:) * tmp2
  enddo

  qclm = qclm * twopia 
  qslm = qslm * twopia


  
  !
  ! Build the average charge density table

  do ii = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     izz = grid%kz(ii)
     rrho = rho(ii-ioffset)
     if (izz > (nzmin-1) .and. izz < (nzmax+1)) then
        sig_counter(izz) = sig_counter(izz) + 1
        sigma(izz) = sigma(izz) + rrho
     endif
  enddo

  !
  ! Global sum for colecting the sigmas from each processor.
  ! and calculate the average on each zz.

  call psum(sigma,nz,parallel%group_size,parallel%group_comm)
  call psum(sig_counter,nz,parallel%group_size,parallel%group_comm)

  ! The potential created by the average desity equals -2pi*avg_sigma*zdiff
  ! and that should be multiplied by the z axis grid step.

  do zindex = nzmin, nzmax
     if (sig_counter(zindex) > 0) then
        sigma(zindex) = -pi2*grid%step(3)*sigma(zindex)/sig_counter(zindex)
     endif
  enddo


  !
  ! Build the tlm table for zdiff exponents
  ! add in all the coefficients for the potential - 2*hcub
  
  do ll = -lpole, lpole
      do mm = 1, lpole
        do zindex = 1, nzneighs-1
                tlm(ll,mm,zindex) = exp(-kklm(ll,mm)*zindex*grid%step(3))
        enddo
        tlm(ll,mm,0) = 1
      enddo
  enddo

  do ll = 1, lpole
     do zindex = 1, nzneighs-1
        tmp2 = exp(-kklm(ll,0)*zindex*grid%step(3)) ! kklm(-ll,0) = kklm(ll,0)
        tlm(ll,0,zindex) = tmp2
        tlm(-ll,0,zindex) = tmp2
     enddo
     tlm(ll,0,0) = 1
     tlm(-ll,0,0) = 1
  enddo

  tlm = tlm * grid%hcub * two

  !
  ! Start with the boundary rho (brho) that equals 8*pi*rho
  ! (4pi because of the poisson's equation 
  ! and 2 to transfer from Hartree to Rydberg)
  do ii = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     brho(ii-ioffset) = rho(ii-ioffset)*pi8
  enddo

  !
  !  Theoreticaly, We need to go through each neighbour in all 
  !  necessary directions (+/-x, +/-y, +/-z) for the order needed.
  !  In practice, only the +/-z directions are of interest in this case. 
  !
  
  do ii = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     do jshell = 1, norder
        do jrow = 5, 6  
           !
           !  If this neighbor is outside the sphere then its neibs value
           !  will be greater than ndim. We can then find the position of
           !  this point and calculate its contribution to brho.
           !
           if (parallel%neibs(jrow+(jshell-1)*lap_neibs,ii-ioffset) &
                > parallel%nwedge) then
              rr(1) = (grid%shift(1) + grid%kx(ii)) * grid%step(1)
              rr(2) = (grid%shift(2) + grid%ky(ii)) * grid%step(2)
              izz = grid%kz(ii)

              select case(jrow)
              case (5)
                 izz = izz - jshell
              case (6)
                 izz = izz + jshell
              end select

              cpole = zero

              ! Initiate the clm and slm tables for the given grid point
              call cos_lm(rr,lpole,klm,clm(-lpole:lpole, 0:lpole))
              call sin_lm(rr,lpole,klm,slm(-lpole:lpole, 0:lpole))

              do zindex = nzmin, nzmax
                 izdiff = abs(izz - zindex)
                 cpole = cpole + sigma(zindex)*izdiff*grid%step(3)

                 do ll = -lpole, lpole
                    do mm = 1, lpole
                       tmp2 = clm(ll,mm)*qclm(ll,mm,zindex) + &
                              slm(ll,mm)*qslm(ll,mm,zindex)
                       cpole = cpole + tlm(ll,mm,izdiff)*tmp2
                    enddo
                 enddo

                 do ll = 1, lpole
                    tmp2 = clm(ll,0)*qclm(ll,0,zindex) + &
                           slm(ll,0)*qslm(ll,0,zindex)
                    cpole = cpole + tlm(ll,0,izdiff)*tmp2
                 enddo

              enddo     ! zindex       


              !
              !  Add the contributions into brho, but scale by 2 for Rydberg
              !  and by coeke (second derivative of the potential-
              !  to obtain the denity)
              !
              brho(ii-ioffset) = brho(ii-ioffset) - &
                   two*coeke(jshell,(jrow+1)/2)*cpole

             !rs(1) = (grid%shift(1) + grid%kx(ii)) * grid%step(1)
             !rs(2) = (grid%shift(2) + grid%ky(ii)) * grid%step(2)
             !rs(3) = (grid%shift(3) + grid%kz(ii)) * grid%step(3)
             ! write(15,106) rs(1), rs(2), rs(3), jshell,   cpole

           endif

        enddo           ! jrow 
     enddo              ! jshell
  enddo                 ! ii

  ! deallocate temp work arrays


  if(allocated(qslm)) deallocate(qslm)
  if(allocated(qclm)) deallocate(qclm)
  if(allocated(tlm)) deallocate(tlm)
  if(allocated(sigma)) deallocate(sigma)
  if(allocated(sig_counter)) deallocate(sig_counter)


!106 format(3(f15.2,1x),(i,1x),3(f15.9,1x))
!107 format(f15.2, 1x, 2(i,1x), 2(e15.7,1x))
!close(15)
!close(16)

end subroutine hartset_slab

!===============================================================
!
!  Functions for the calculation of cosine and sine functions
!  rr must be of length 2.
!
!

subroutine cos_lm(rr,n,vec_in,vec_out)
  use constants
  implicit none

  real(dp), intent(in) :: rr(2)
  integer, intent(in) :: n
  real(dp), intent(in) :: vec_in(-n:n,0:n,2)
  real(dp), intent(out) :: vec_out(-n:n,0:n)
  real(dp) :: mult
  integer :: ll, mm

  do ll = -n, n
     do mm = 1, n
        mult = dot_product(rr,vec_in(ll,mm,:))
        vec_out(ll,mm) = cos(mult)
     enddo
  enddo

  do ll = 1, n
     mult = dot_product(rr,vec_in(ll,0,:))
     vec_out(ll,0) = cos(mult)
  enddo

  vec_out(0,0) = 0
  return

end subroutine cos_lm
!==============================================================

subroutine sin_lm(rr,n,vec_in,vec_out)
  use constants
  implicit none

  real(dp), intent(in) :: rr(2)
  integer, intent(in) :: n
  real(dp), intent(in) :: vec_in(-n:n,0:n,2)
  real(dp), intent(out) :: vec_out(-n:n,0:n)
  real(dp) :: mult
  integer :: ll, mm

  do ll = -n, n
     do mm = 1, n
        mult = dot_product(rr,vec_in(ll,mm,:))
        vec_out(ll,mm) = sin(mult)
     enddo
  enddo

  do ll = 1, n
     mult = dot_product(rr,vec_in(ll,0,:))
     vec_out(ll,0) = sin(mult)
  enddo

  vec_out(0,0) = 0
  return
 
end subroutine sin_lm
!==============================================================
