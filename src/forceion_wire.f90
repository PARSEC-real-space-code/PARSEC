!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Computes the ionic energy and corresponding force in wire
!  geometry by doing a summation over the reciprocal space.
!
!  author: Murilo Tiago, UTexas, February 2007
!
!---------------------------------------------------------------
subroutine forceion_wire(clust,alatt,enuc,ipr,zv,ierr)

  use constants
  use cluster_module
  implicit none
  !
  !  Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(inout) :: clust
  !  lattice parameter (length of periodic cell along x direction)
  real(dp), intent(in) :: alatt
  !  print flag
  integer, intent(in) :: ipr
  !  Ionic energy
  real(dp), intent(out) :: enuc
  !  Electric charge for each chemical element
  real(dp), intent(in) :: zv(clust%type_num)
  !  Output error, 500 < ierr < 511
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  number of G vectors in the reciprocal space, one-dimensional
  integer :: ng
  !  2 * Pi / lattice_parameter
  real(dp) :: twopia
  !  2 * Pi * G / lattice_parameter
  real(dp) :: twopina
  !  2 * 4 * Z * Z' / lattice_parameter
  !  extra factor of 2 from hartrees -> rydbergs conversion
  real(dp) :: zovera
  !  counters
  integer :: ja, jap, ity, iat, ig
  !  vector quantities:
  !
  real(dp) :: &
       rr(3), &    ! rr = ratom(j) - ratom(j')
       rx, &       ! component of rr along the x axis
       rperp       ! normalized component of rr perpendicular to the x axis
  !  arguments and values of modified Bessel functions
  real(dp) :: arg_k, arg_0, bessk(0:1)
  !  complex phases and sum accumulators
  complex(dpc) :: phase0, expi, esum0, fsum_x, fsum_p
  !  temporary arrays for the ionic forces and coordinates
  real(dp), dimension(3,clust%atom_num) :: fsum, ratom
  !  temporary array for the electric charge on each ion
  real(dp) :: zz(clust%atom_num)
  !  external functions
  real(dp), external :: expint

  !  constants
  real(dp), parameter :: small = 1.d-20, r_tol = 2.5d-5
  real(dp), parameter :: eta = 1.d4, euler = 0.57721566490153286d0

  !---------------------------------------------------------------

  ng = nint ( alatt * 20.0 / two / pi / 1.d-5)
  twopia = twopi / alatt
  ierr = 0

  !  Store temporary arrays for atom coordinates and charges.
  ja = 0
  do ity = 1, clust%type_num
     do iat = 1, clust%natmi(ity)
        ja = ja + 1
        zz(ja) = zv(ity)
        ratom(:,ja) = (/ clust%xatm(ja), clust%yatm(ja), clust%zatm(ja) /)
     enddo
  enddo

  !  Zero-point energy (interaction between each atom and its images)
  !  contains only the log divergence, which canceled with the Hartree 
  !  term
  enuc = zero

  !  Start sum in g space.
  fsum(:,:) = zero
  do ja = 1, clust%atom_num
     do jap = 1, ja - 1
        zovera = eight*zz(ja)*zz(jap)/alatt
        rr = ratom(:,ja) - ratom(:,jap)
        rx = rr(1)
        rperp = sqrt( rr(2)*rr(2) + rr(3)*rr(3) )
        if (rperp < r_tol) then
           phase0 = zi * twopia * rx
           phase0 = exp( phase0 )
           expi = zone
           esum0 = log(eta*two*alatt) * half - euler / four
           fsum_x = zzero
           fsum_p = zzero
           arg_k = zero
           arg_0 = pi/alatt/eta
           arg_0 = arg_0 * arg_0
           twopina = zero
           do ig = 1, ng
              arg_k = arg_0 * real(ig,dp) * real(ig,dp)
              arg_k = expint(arg_k,ierr) * half
              expi = expi * phase0
              twopina = twopina + twopia
              esum0 = esum0 + expi * arg_k
              fsum_x = fsum_x + expi * twopina * arg_k
              if (arg_k < small) exit
           enddo
        else
           rr = rr / rperp

           phase0 = zi * twopia * rx
           phase0 = exp( phase0 )
           expi = zone
           esum0 = log(two*alatt/rperp) / two - euler/two
           fsum_x = zzero
           fsum_p = zone * half / rperp
           arg_k = zero
           arg_0 = twopia * rperp
           twopina = zero
           do ig = 1, ng
              arg_k = arg_k + arg_0
              call bessel_k(1,arg_k,bessk)
              expi = expi * phase0
              twopina = twopina + twopia
              esum0 = esum0 + expi * bessk(0)
              fsum_x = fsum_x + expi * twopina * bessk(0)
              fsum_p = fsum_p + expi * twopina * bessk(1)
              if (bessk(0) < small .and. bessk(1) < small) exit
           enddo
        endif

        if (ig >= ng) then
           write(7,*) ' ERROR in forceion_wire ! Radial distance between ', &
                'ions ',ja,' and ',jap,' is too small.'
           write(7,*) ' increase value of ng. ',ig,ng
           ierr = 502
        endif

        enuc = enuc + zovera * real(esum0,dp)
        fsum(1,ja) = fsum(1,ja) + zovera * aimag(fsum_x)
        fsum(1,jap) = fsum(1,jap) - zovera * aimag(fsum_x)
        fsum(2:3,ja) = fsum(2:3,ja) + zovera * real(fsum_p,dp) * rr(2:3)
        fsum(2:3,jap) = fsum(2:3,jap) - zovera * real(fsum_p,dp) * rr(2:3)
     enddo
  enddo
  clust%force = fsum
  !  If print flag on, print the ion-ion forces.
  if (ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion-ion interaction:'
     write(7,*) '================================'
     write(7,22)
     write(7,*)
     do ja = 1, clust%atom_num
        write(7,20) ja,clust%xatm(ja),clust%yatm(ja),clust%zatm(ja), &
             clust%force(:,ja)
     enddo
     write(7,*)
     write(7,'(a,f20.8,a)') ' Ion-Ion Energy = ',enuc,' [Ry]'
  endif

20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ')

end subroutine forceion_wire
!===============================================================
!
!  Computes the exponential integral E1(x), for x > 0
!
!  (C) Copr. 1986-92 Numerical Recipes Software #0).
!
!---------------------------------------------------------------------
function expint(x,ierr)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  real(dp), intent(in) :: x
  real(dp) :: expint
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  constants:
  integer, parameter :: maxit = 1000000
  real(dp), parameter :: eps = 1.d-20, fpmin = 1.d-39
  real(dp), parameter :: euler = 0.57721566490153286d0

  ! counters
  integer :: jj
  real(dp) :: del, fact
  !---------------------------------------------------------------------
  ierr = 0
  if (x <= zero ) then
     write(9,*) 'ERROR: bad arguments in expint'
     expint = zero
     ierr = 1
     goto 10
  else
     expint = -log(x) - euler
     fact = one
     do jj = 1, maxit
        fact = -fact*x/real(jj,dp)
        del = -fact/real(jj,dp)
        expint = expint + del
        if (abs(del) < abs(expint)*eps) goto 10
     enddo
     write(9,*) 'ERROR: series failed in expint'
     ierr = 2
  endif
10 continue

end function expint
!=====================================================================
