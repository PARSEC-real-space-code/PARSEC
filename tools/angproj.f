c     ===============================================================
c
c     Copyright (C) 2005 Finite Difference Research Group
c     This file is part of parsec, http://www.ices.utexas.edu/parsec/
c
c     This program projects DFT wave-functions into s, p, and d
c     partial waves by integrating them over spheres of radius rdh
c     around atom sites. At the end, the density of states (DOS) and
c     angle-projected DOS (angdos) are calculated using Gaussian
c     convolution. The total angdos (summed over all partial waves)
c     gives a measure of how much each wave-function is concentrated
c     around the atom sites (if total angdos = DOS, the wave-function
c     is totally concentrated within rdh of the atom sites).
c
c     INPUT: parsec.dat (DFT wave-functions from a parsec run)
c            angproj.in (input atomic coordinates and parameters, see example)
c
c     input parameters:
c     natom = number of atoms
c     pbcflag = 1 for periodic system, 0 otherwise
c     rdh = muffin-tin, ion radius (a.u.), should be similar to the pseudop.
c     cut-off radius or sligthly bigger.
c     eta = dispersion in Gaussian convolution (eV); for better DOS,
c     should be of order of the typical spacing of eigenvalues.
c
c     OUTPUT: angproj.dat (calculated projections)
c             angdos.dat (DOS and angdos)
c
c     each line in angproj.dat refers to an electronic state, columns
c     are:
c       column 1 : state index
c              2 : spin index
c              3 : energy, in eV
c              4 : electron occupancy
c              5 : muffin-tin projection, integral of wave-function squared
c                  over spheres around all atoms
c              6 : s projection, normalized over muffin-tin projection
c              7 : p projection, normalized
c              8 : d projection, normalized
c
c     columns in angdos.dat are:
c       column 1 : energy, in eV
c              2 : total DOS, in eV^-1, normalized so that its integral
c                  over energy is the total number of states
c              3 : LDOS, projected over muffin-tin spheres
c              4 : s component of LDOS
c              5 : p component of LDOS
c              6 : d component of LDOS
c
c     author : Murilo Tiago (Univ. of Minnesota, February 2005).
c
c     ---------------------------------------------------------------
c
c     Constants
c
      module constants
      integer, parameter :: dp = kind(1.0d0)
      integer, parameter :: dpc = kind((1.0d0,1.0d0))
      real(dp), parameter :: pi = 3.1415926535897932384626433832795d0
      real(dp), parameter :: ryd = 13.60569172d0
      integer, parameter :: nproj = 9 ! number of angular projections
      integer, parameter :: nw = 1000 ! number of mesh points in DOS, angdos
      integer, parameter :: npt = 3
      end module constants
c
c     ---------------------------------------------------------------
      program angproj
c
      use constants
      implicit none

      integer ::
     1     nspin,               ! number of spins
     2     nstate,              ! number of DFT wave-functions
     3     ndim,                ! number of grid points (irreducible wedge)
     4     ntrans,              ! number of symmetry operations
     5     natom,               ! number of atoms to project onto
     6     pbcflag,             ! = 1 for periodic boundary conditions,
                                ! = 0 otherwise.
     7     ngl,                 ! number of radial points in Gauss-Legendre 
                                ! integration (10 is generally a good value)
     8     ncplx,               ! real(0) or complex(1) wave-functions
     9     nkpt,                ! number of k-points
     1     nreplica             ! number of atom replicas (periodic systems only)

      real(dp) ::
     1     rdh,                 ! "radius" of atomic spheres
     2     sigma,               ! dispersion of radial shell (see subroutine angular)
     3     hh(3),               ! grid spacing for x,y,z directions
     4     rsize,               ! size of grid for x,y,z directions
     5     shift(3),            ! grid shift from parsec.dat
     6     g_amp(0:nproj),      ! angular amplitudes
     7     eta,                 ! dispersion in Gaussian convolution (angdos)
     8     emin,emax,           ! min,max value of energy in DOS, angdos
     9     angdos(0:nproj),     ! angular density of states
     1     avec(3,3),           ! unit lattice vectors
     2     anorm(3,3)           ! normalizedunit lattice vectors

      integer, allocatable ::
     1     rgrid(:,:),          ! coordinates of grid points, in units of hh
     2     istate(:,:),         ! number of wave-functions per spin component
     3     irep(:,:,:),         ! irreducible representation of wave-functions
     4     chi(:,:),            ! character table
     5     indx(:,:,:)          ! indices of printed wave-functions

      real(dp), allocatable ::
     1     ratm(:,:),           ! atomic coordinates, in a.u.
     2     trans(:,:,:),        ! rotation matrices from symmetry operations
     2     eigen(:,:,:),        ! DFT eigenvalues, converted to eV
     3     occup(:,:,:),        ! occupancy factors of DFT eigenstates
     4     wfn(:),              ! wave-functions
     5     proj_at(:,:,:,:,:),  ! projection of each eigenstate onto each 
                                ! atom and orbital
     6     proj(:,:,:,:),       ! summed projection for all atoms
     7     rgl(:),              ! radial positions in Gauss-Leg. integral
     8     wgl(:),              ! radial weights in Gauss-Legendre integral
     9     pgl(:,:),            ! integration array
     1     kw(:)                ! weights of k-points

      complex(dpc), allocatable ::
     1     zwfn(:),             ! wave-functions
     2     zpgl(:,:)            ! integration array

c     counters
      integer ii, jj, is, ikp, jr, iat, ip, jl, itrans, icell1,
     1     icell2, icell3
      real(dp) :: rr(3), rt(3), dist, tsum, h32, pr_t, pr_s, pr_p, pr_d,
     1     omega, dos, fac, fac2, mtmp(3,3), rmin(3), dmin

c     ---------------------------------------------------------------

c
c     Read angproj.in
c
      open(10,file='angproj.in',status='old')
      read(10,*) natom, rdh, eta

      allocate(ratm(3,natom))
      do jj=1,natom
         read(10,*) (ratm(ii,jj),ii=1,3)
      enddo
      close(10)
c
c     Start reading parsec.dat
c
      open(20,file='parsec.dat',form='unformatted')

      read(20)
      read(20) nspin, ncplx, pbcflag
      if (pbcflag .eq. 1) then
         write(6,*) ' angle-resolved DOS for periodic system'
         read(20) hh(:)
         read(20) avec(:,:)
         do ii = 1, 3
            read(20)
         enddo
         read(20) nkpt
         allocate(kw(nkpt))
         do ii = 1, 4
            read(20)
         enddo
         read(20) (kw(ii), ii=1,nkpt)
      else
         write(6,*) ' angle-resolved DOS for confined system'
         read(20) hh(1),rsize
         hh(2:3) = hh(1)
         nkpt = 1
         allocate(kw(nkpt))
         kw = 1.d0
         avec = 0.d0
         avec(1,1) = rsize + 2.d0 * hh(1)
         avec(2,2) = rsize + 2.d0 * hh(1)
         avec(3,3) = rsize + 2.d0 * hh(1)
      endif
ccm
      if (nkpt .eq. 0) then
         nkpt = 1
         if (allocated(kw)) deallocate(kw)
         allocate(kw(nkpt))
         kw = 1.d0
      endif
ccm
      if (ncplx .eq. 1) then
          write(6,*) ' wave-functions are complex'
      else
          write(6,*) ' wave-functions are real'
      endif
      read(20)
      read(20) shift
      do ii = 1, 7
         read(20)
      enddo
c
c     Renormalize lattice vectors
c
      anorm = avec
      do ii = 1, 3
         fac = sqrt( dot_product(anorm(:,ii),anorm(:,ii)) )
         anorm(:,ii) = anorm(:,ii)/fac
      enddo
      mtmp = anorm
      call mtrxin(mtmp,fac,fac2)

      h32 = sqrt( hh(1) * hh(2) * hh(3) * fac )
      sigma = minval( hh )
c
c     search for maximum number of wave-functions
      do is = 1,nspin
         do ikp = 1, nkpt
            read(20) nstate
            if (is .eq. 1 .and. ikp .eq. 1) then
               allocate(irep(nstate,nkpt,nspin))
               allocate(eigen(nstate,nkpt,nspin))
               allocate(occup(nstate,nkpt,nspin))
               allocate(istate(nkpt,nspin))
               allocate(indx(nstate,nkpt,nspin))
            endif
            read(20) (irep(ii,ikp,is),ii=1,nstate)
            read(20) (eigen(ii,ikp,is),ii=1,nstate)
            read(20) (occup(ii,ikp,is),ii=1,nstate)
         enddo
         read(20)
         read(20)
      enddo
      eigen = eigen * ryd

c     number of points in radial integration, npt should be small for
c     fast calculations but no less than about 2
      ngl = npt * nint( rdh/sigma )
c
      allocate(rgl(ngl))
      allocate(wgl(ngl))
      if (ncplx .eq. 0) then
         allocate(pgl(ngl,0:nproj))
      else
         allocate(zpgl(ngl,0:nproj))
      endif

      call mygauleg(0,rdh,rgl,wgl,ngl)

      write(6,14) ' number of atoms = ',natom
      write(6,13) ' atomic radius = ',rdh,' a.u.'
      write(6,14) ' number of points in radial integration = ',ngl
      write(6,13) ' dispersion in radial integrals = ',sigma,' a.u.'
      write(6,13) ' dispersion in convoluted DOS = ',eta,' eV'
 13   format(a,f8.4,a)
 14   format(a,i8)
c
c     Initialize arrays
c
      allocate(proj_at(nstate,nkpt,nspin,natom,0:nproj))
      proj_at = 0.d0
      allocate(proj(nstate,nkpt,nspin,0:nproj))
      proj = 0.d0
c
c     Start calculation
c

c     Define parameters for periodic boundary conditions.
      if (pbcflag .eq. 1) then
         nreplica = 1
      else
         nreplica = 0
      endif

      if (natom .eq. 0) then
         istate = nstate
         do ii = 1, nstate
            indx(ii,:,:) = ii
         enddo
         close(20)
         goto 21
      endif
      do is = 1, nspin
         do ikp = 1, nkpt
            read(20) ndim, ntrans
            read(20)
            allocate(trans(3,3,ntrans))
            read(20) (trans(:,:,jj),jj=1,ntrans)
            read(20)
            read(20)
            allocate(chi(ntrans,ntrans))
            read(20) ((chi(ii,jj),ii=1,ntrans),jj=1,ntrans)
            allocate(rgrid(3,ndim))
            read(20) ((rgrid(ii,jj),ii=1,3),jj=1,ndim)
            read(20) istate(ikp,is)
            if (istate(ikp,is) .eq. 0) then
               deallocate(trans)
               deallocate(chi)
               deallocate(rgrid)
               cycle
            endif
            read(20) (indx(ii,ikp,is),ii=1,istate(ikp,is))
            if (ncplx .eq. 0) then
               allocate(wfn(ndim))
            else
               allocate(zwfn(ndim))
            endif
            do jj = 1,istate(ikp,is)
               if (ncplx .eq. 0) then
                  read(20) (wfn(jr), jr=1,ndim)
               else
                  read(20) (zwfn(jr), jr=1,ndim)
               endif
               do iat = 1, natom
c     initialize radial integration
                  if (ncplx .eq. 0) then
                     pgl = 0.d0
                  else
                     zpgl = 0.d0
                  endif

                  do jr = 1, ndim
                     do itrans = 1, ntrans
c     Get the position of current atom with respect to each grid
c     point. Look for the atom closer to this grid point (taking
c     cell periodicity into account if system is periodic).
                        rr = rgrid(:,jr) + shift
                        do ii = 1, 3
                           rr(ii) = hh(ii) * rr(ii)
                        enddo
                        rr = matmul(anorm,rr)
                        rt = matmul(rr,trans(:,:,itrans))
                        dmin = 9.d9
                        do icell1 = -nreplica,nreplica
                           do icell2 = -nreplica,nreplica
                              do icell3 = -nreplica,nreplica
                                 rr = ratm(:,iat) - rt
     1                                + avec(:,1)*real(icell1,dp)
     2                                + avec(:,2)*real(icell2,dp)
     3                                + avec(:,3)*real(icell3,dp)
                                 dist = DOT_PRODUCT(rr,rr)
                                 if (dist .lt. dmin) then
                                    rmin = rr
                                    dmin = dist
                                 endif
                              enddo
                           enddo
                        enddo
                        rr = rmin
                        dist = (sqrt(DOT_PRODUCT(rr,rr)) - rdh)/sigma

c     if this point is too far away, skip it
                        if (dist .gt. 5.d0) cycle

c     get angular amplitudes at this grid point and sum over the grid
c     pgl holds the total amplitude (summed over all angular momenta)
                        do jl = 1, ngl
                           call angular(rr,rgl(jl),sigma,g_amp)
                           if (ncplx .eq. 0) then
                           pgl(jl,1:nproj) = pgl(jl,1:nproj) +
     1                          g_amp(1:nproj) * wfn(jr) * h32 *
     2                          chi(irep(jj,ikp,is),itrans)
                           pgl(jl,0) = pgl(jl,0)+g_amp(0)*wfn(jr) **2
                        else
                           zpgl(jl,1:nproj) = zpgl(jl,1:nproj) +
     1                          g_amp(1:nproj) * zwfn(jr) * h32 *
     2                          chi(irep(jj,ikp,is),itrans)
                           zpgl(jl,0) = zpgl(jl,0)+g_amp(0)*zwfn(jr) **2
                        endif
                        enddo
                     enddo      ! itrans = 1, ntrans

                  enddo
c     radial integration
                  tsum = 0.d0
                  if (ncplx .eq. 0) then
                     do jl = 1, ngl
                        tsum = tsum + wgl(jl)*(rgl(jl))**2*pgl(jl,0)
                     enddo
                  else
                     do jl = 1, ngl
                        tsum = tsum + wgl(jl)*(rgl(jl))**2
     1                       *abs( zpgl(jl,0) )
                     enddo
                  endif
                  proj_at(jj,ikp,is,iat,0) = tsum
                  do ip = 1, nproj
                     tsum = 0.d0
                     if (ncplx .eq. 0) then
                        do jl = 1, ngl
                           tsum = tsum + wgl(jl)*
     1                          (rgl(jl)*pgl(jl,ip))**2
                        enddo
                     else
                        do jl = 1, ngl
                           tsum = tsum + wgl(jl)*
     1                          abs( (rgl(jl)*zpgl(jl,ip))**2 )
                        enddo
                     endif
                     proj_at(jj,ikp,is,iat,ip) = tsum
                  enddo

               enddo            ! iat = 1, natom

            enddo               ! jj = 1,istate(is)
            deallocate(trans)
            deallocate(chi)
            deallocate(rgrid)
            if (ncplx .eq. 0) then
               deallocate(wfn)
            else
               deallocate(zwfn)
            endif
         enddo                  ! ikp = 1, nkpt
      enddo                     ! is = 1, nspin
      close(20)
c
c     Sum over atoms
c
      do is = 1, nspin
         do ikp = 1, nkpt
            do jj = 1, istate(ikp,is)
               do ip = 0, nproj
                  proj(jj,ikp,is,ip) = sum( proj_at(jj,ikp,is,:,ip) )
               enddo
            enddo
         enddo
      enddo
c
c     Print out projections
C     (note: projections are not weighted by occupancy factor!)
c
      open(30,file='angproj.dat')

      do is = 1, nspin
         do ikp = 1, nkpt
            do jj = 1, istate(ikp,is)
               ii = indx(jj,ikp,is)
               pr_t = proj(jj,ikp,is,0)
               pr_s = proj(jj,ikp,is,1)/proj(jj,ikp,is,0)
               pr_p = sum( proj(jj,ikp,is,2:4) )/proj(jj,ikp,is,0)
               pr_d = sum( proj(jj,ikp,is,5:9) )/proj(jj,ikp,is,0)
               write(30,'(3(i4,1x),f12.6,1x,5(f10.4))')
     1              ii,ikp,is,eigen(ii,ikp,is),
     2              occup(ii,ikp,is),pr_t,pr_s,pr_p,pr_d
            enddo
         enddo
      enddo
      close(30)

 21   continue
c
c     Calculate angle-resolved density of states using convolution of
c     original density of states with a Gaussian function.
c
      emin = minval(eigen) - 10.0*eta
      emax = maxval(eigen) + 10.0*eta

      if (nspin .eq. 1) then
         open(41,file='angdos.dat',form='formatted')
      else
         open(41,file='angdos_up.dat',form='formatted')
         open(42,file='angdos_down.dat',form='formatted')
      endif

      do is = 1, nspin
         ip = 40 + is
         do ii = 1, nw
            omega = emin + (emax - emin)*real(ii,dp)/real(nw,dp)
            angdos = 0.d0
            dos = 0.d0
            do ikp = 1, nkpt
               do jj = 1, istate(ikp,is)
                  fac = omega - eigen(indx(jj,ikp,is),ikp,is)
                  fac2 = exp( -fac*fac/(2.d0*eta**2) )
                  fac2 = fac2/sqrt(pi*2.d0)/eta
                  dos = dos + fac2
                  angdos(:) = angdos(:) + fac2 * proj(jj,ikp,is,:)
     1                 * kw(ikp)
               enddo
            enddo
            pr_t = angdos(0)
            pr_s = angdos(1)
            pr_p = sum( angdos(2:4) )
            pr_d = sum( angdos(5:9) )
            write(ip,'(8(f11.5))') omega,dos,pr_t,pr_s,pr_p,pr_d
         enddo
      enddo

      close(41)
      if (nspin .eq. 2) close(42)

      end program angproj

c     ===============================================================
c
c     For a given vector rr(1:3), calculates amplitudes of angular
c     functions, weighted along the radial direction by a Gaussian-like
c     delta function centered at distance rshell and with dispersion
c     sigma. Angular dependence of output amplitudes is:
c     g_amp(0) : 1.0 (used for total amplitude)
c          (1) : s component
c          (2) : p_x component
c          (3) : p_y component
c          (4) : p_z component
c          (5) : d_xy component
c          (6) : d_(x^2 - y^2) component
c          (7) : d_xz component
c          (8) : d_yz component
c          (9) : d_(3z^2 - r^2) component
c
c     ---------------------------------------------------------------
      subroutine angular(rr,rshell,sigma,g_amp)

      use constants
      implicit none
c
c     Input/Output variables:
c
c     position
      real(dp), intent(in) :: rr(3)
c     central point and dispersion of Gaussian function
      real(dp), intent(in) :: rshell,sigma
c     amplitude of angular projections (delta function * spherical
c     harmonic) at this position
      real(dp), intent(out) :: g_amp(0:nproj)
c
c     Work variables:
c
      real(dp) :: dist, rg, alp, pref_0, pref_1, pref_2, pref_g,
     1     rnorm(3)
      real(dp), parameter :: rtol = 1.d-4 ! tolerance in distance

c     ---------------------------------------------------------------
c
c     parameters
c
      pref_0 = sqrt(1.d0/(4.d0*pi))
      pref_1 = sqrt(3.d0/(4.d0*pi))
      pref_2 = sqrt(15.d0/(4.d0*pi))
      pref_g = sqrt(2.d0 * pi) * sigma
c
c     normalize position vector
      dist = sqrt( DOT_PRODUCT(rr,rr) )
c
c     if distance is too short, angles can not be defined
      if ( dist .lt. rtol ) then
         g_amp = 0.d0
         return
      endif
      rnorm = rr/dist

c     the radial part of atomic wave-function for s,p and d
      alp = (dist - rshell)**2/2.d0/sigma**2
      rg = exp(-alp) / dist**2 / pref_g
c     partial waves
      g_amp(0) = rg
      g_amp(1) = rg*pref_0
      g_amp(2) = rg*pref_1*rnorm(1)
      g_amp(3) = rg*pref_1*rnorm(2)
      g_amp(4) = rg*pref_1*rnorm(3)
      g_amp(5) = rg*pref_2*rnorm(1)*rnorm(2)
      g_amp(6) = rg*pref_2*(rnorm(1)**2 - rnorm(2)**2)/(2.d0)
      g_amp(7) = rg*pref_2*rnorm(1)*rnorm(3)
      g_amp(8) = rg*pref_2*rnorm(2)*rnorm(3)
      g_amp(9) = rg*pref_2*(3.d0*rnorm(3)**2 - 1.d0)/(2.d0*sqrt(3.d0))

      end subroutine angular

c     ===============================================================
c
c     Given the lower and upper limits of integration x1 and x2, and
C     given n, this routine returns arrays x(1:n) and w(1:n),
C     containing the abscissas and weights of the Gauss-Legendre
C     n-point quadrature formula.
c     from Numerical Recipes in Fortran
c
c     ---------------------------------------------------------------
      subroutine mygauleg(x1, x2, x, w, n)  

      use constants
      implicit none  
c
c     Input/Output variables:
c
      real(dp), intent(in) :: x1, x2  
      integer, intent(in) :: n  
      real(dp), intent(out) :: w(n), x(n)  
c
c     Work variables:
c
      real(dp) :: p1, p2, p3, pp, xl, xm, z, z1, dj
      integer :: i, j, m  
      real(dp), parameter :: eps = 3.0d-14  

c     ---------------------------------------------------------------

      m = (n + 1) / 2  
      xm = 0.5d0 * (x2 + x1)  
      xl = 0.5d0 * (x2 - x1)  
      do i = 1, m  
         z = cos(pi * (real(i, dp) - 0.25d0) / (real(n, dp) + 0.5d0))
         do
            p1 = 1.d0
            p2 = 0.d0
            do j = 1, n  
               p3 = p2  
               p2 = p1
               dj = real(j, dp)
               p1 = ((2.d0*dj - 1.d0)*z*p2 - (dj - 1.d0)*p3) / dj
            end do
            pp = real(n, dp) * (z * p1 - p2) / (z * z - 1.d0)  
            z1 = z  
            z = z1 - p1 / pp  
            if (abs(z - z1) .lt. eps) exit
         end do
         x(i) = xm - xl * z  
         x(n + 1 - i) = xm + xl * z  
         w(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)  
         w(n + 1 - i) = w(i)  
      end do

      end subroutine mygauleg
c     ===============================================================
c
c     Inverts 3x3 matrix m corresponding to a symmetry operation,
c     storing the result in m. It also calculates determinant and
c     trace of the input matrix. Matrix inversion is aborted if
c     det<del.
c
c     ---------------------------------------------------------------
      subroutine mtrxin(m,det,tr)

      use constants
      implicit none
c
c     Input/Output variables:
c
      real(dp), intent(inout) :: m(3,3)
      real(dp), intent(out) :: det,tr
c
c     Work variables:
c
      real(dp) :: a(3,3),del,x
      integer i,j
c     ---------------------------------------------------------------
c
c     compute matrix of cofactors
c
      a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
      a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
      a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
      a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
      a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
      a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
      a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
      a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
      a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
c
c     compute determinant
c
      det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
      tr = m(1,1) + m(2,2) + m(3,3)
      del = 1.0d-05
      if (abs(det) < del) stop 501
c
c     form mi
c
      do i=1,3
         do j=1,3
            x = a(i,j)/det
            m(i,j) = x
c            if (x < 0.0d0) m(i,j) = x - del
c            if (x > 0.0d0) m(i,j) = x + del
         enddo
      enddo

      end subroutine mtrxin
c     ===============================================================

