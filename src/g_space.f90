!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Calculates g-vectors (kgv) for lattice points inside a G-grid.
! G-vectors are ordered by stars depending on symmetry and pairs of
! (+G,-G) vectors are put together. Notice that
! the G-grid is defined as (-kmax,kmax), and the number of grid points
! along each direction is odd by construction. But the FFT box is
! defined as (m,m + N - 1), where N is *even*. In the parent
! routine, kmax is defined so that 2*kmax + 1 >= N.
!
! Adapted from plane-wave programs written by S. Froyen and
! J. L. Martins.
!
!---------------------------------------------------------------
subroutine g_space(pbc,symm,kmax,ipr,ierr)

  use constants
  use pbc_module
  use symmetry_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  ! symmetry operations
  type (symmetry), intent(in) :: symm
  ! half-size of reciprocal space
  integer, intent(inout) :: kmax(3)
  ! print flag
  integer, intent(in) :: ipr
  ! error flag, 260 < ierr < 266
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  integer ns,ng
  integer i,j,k,l,istar,istop,ir,itest,iprot,iavail,ncount
  integer kk(3)
  integer, dimension (:), allocatable :: izstar
  ! number of g-vectors per star (transfered to pbc%mstar)
  integer, dimension(:), allocatable :: tmstar
  ! kinetic energy for each star (transfered to pbc%ek)
  real(dp), dimension(:), allocatable :: tek

  integer ktran(3,symm%ntrans)

  real(dp) :: gmod,rgm,gmr,dkk(3),dkt(3)

  real(dp) :: fi,diff,bdot(3,3)
  complex(dpc) :: phtmp
  complex(dpc) :: expt(symm%ntrans)
  real(dp), dimension (:), allocatable :: gm

  real(dp), parameter :: delta = 1.0d-8

  !---------------------------------------------------------------

  bdot = pbc%bdot

  ! g = k1*b1+k2*b2+k3*b3  

  ng=(2*kmax(1)+1)*(2*kmax(2)+1)*(2*kmax(3)+1)

  if(ng <= 1) then
     write(7,56) ng
     ierr = 261
     return
  endif

  if (.not. allocated (gm)) allocate (gm(ng))
  allocate(tek(ng))
  allocate(tmstar(ng))
  call pbc_set_ng (ng,pbc)

  ng = 0
  do i=-kmax(1),kmax(1)
     dkk(1) = real(i,dp)
     do j=-kmax(2),kmax(2)
        dkk(2) = real(j,dp)
        do k=-kmax(3),kmax(3)
           dkk(3) = real(k,dp)
           call matvec3('N',bdot,dkk,dkt)
           gmod = dot_product(dkk,dkt)
           ng = ng + 1
           gm(ng) = gmod
           pbc%kgv(1,ng) = i
           pbc%kgv(2,ng) = j
           pbc%kgv(3,ng) = k
        enddo
     enddo
  enddo

  ! Sort the g-vectors according to length and collect
  l=ng/2+1
  ir=ng
20 continue

  if(l > 1)then
     l=l-1
     rgm=gm(l)
     kk(:)=pbc%kgv(:,l)
  else
     rgm=gm(ir)
     kk(:)=pbc%kgv(:,ir)
     gm(ir)=gm(1)
     pbc%kgv(:,ir)=pbc%kgv(:,1)
     ir=ir-1
     if(ir == 1)then
        gm(1)=rgm
        pbc%kgv(:,1)=kk(:)
        goto 25
     
     endif
  endif
  i=l
  j=l+l

21 continue

  if(j <= ir)then
     if(j < ir)then
        if(gm(j) < gm(j+1))j=j+1
     endif
     if(rgm < gm(j))then
        gm(i)=gm(j)
        pbc%kgv(:,i)=pbc%kgv(:,j)
        i=j
        j=j+j
     else
        j=ir+1
     endif

     go to 21

  endif
  gm(i)=rgm
  pbc%kgv(:,i)=kk(:)

  go to 20

25 continue
  !
  ! Sort g vectors of same length by the kgv(1,i), then kgv(2,i),
  ! then the kgv(3,i) value. this makes the kgv sort machine
  ! independent.
  !
  istar = 2

  do i=2,ng
     if( i == ng) then
        diff = 1.0
     else
        diff = gm(i+1)-gm(istar)
     endif
     if(diff > delta) then
        istop = i
        do j=istar,istop-1
           do k=j+1,istop
              if (pbc%kgv(1,k) < pbc%kgv(1,j)) then
                 rgm = gm(j)
                 kk(:) = pbc%kgv(:,j)
                 gm(j) = gm(k)
                 pbc%kgv(:,j) = pbc%kgv(:,k)
                 gm(k) = rgm
                 pbc%kgv(:,k)=kk(:)
              elseif (pbc%kgv(1,k) == pbc%kgv(1,j)) then
                 if (pbc%kgv(2,k) < pbc%kgv(2,j)) then
                    rgm = gm(j)
                    kk(:) = pbc%kgv(:,j)
                    gm(j) = gm(k)
                    pbc%kgv(:,j) = pbc%kgv(:,k)
                    gm(k) = rgm
                    pbc%kgv(:,k)=kk(:)
                 elseif (pbc%kgv(2,k) == pbc%kgv(2,j)) then
                    if (pbc%kgv(3,k) < pbc%kgv(3,j)) then
                       rgm = gm(j)
                       kk(:) = pbc%kgv(:,j)
                       gm(j) = gm(k)
                       pbc%kgv(:,j) = pbc%kgv(:,k)
                       gm(k) = rgm
                       pbc%kgv(:,k)=kk(:)
                    endif
                 endif
              endif
           enddo
        enddo
        istar = istop + 1
     endif
  enddo
  !
  ! collects g-vectors into stars by symmetry
  ! first star is (0 0 0)
  !
  ns = 1
  tek(1) = zero
  pbc%phase(1) = zone
  pbc%conj(1) = one
  pbc%inds(1) = 1
  tmstar(1) = 1

  iprot = 1
  iavail = 2
  itest = 2

30 continue

  if(itest > ng) then
     diff = one
  else
     diff = abs(gm(itest)-gm(iprot))
  endif
  if(diff > delta) then
     ! new star
     iprot = iavail
     ns = ns + 1
     tek(ns) = gm(iprot) / two
     phtmp = zzero
     ncount = 0
     do i=1,symm%ntrans
        do j = 1, 3
           ktran(j,i) = dot_product(symm%gmtrx(j,:,i),pbc%kgv(:,iprot))
        enddo
        fi = dot_product(symm%tnp(:,i),ktran(:,i)*twopi)
        expt(i) = cmplx(cos(fi),sin(fi),dp)
        if ((ktran(1,i)-pbc%kgv(1,iprot) == 0) .and. &
             (ktran(2,i)-pbc%kgv(2,iprot) == 0) .and. &
             (ktran(3,i)-pbc%kgv(3,iprot) == 0)) then
           phtmp = phtmp + expt(i)
           ncount = ncount + 1
        endif
     enddo
     pbc%phase(iprot) = phtmp/real(ncount,dp)
     pbc%conj(iprot) = one
     tmstar(ns) = 1
     pbc%inds(iprot) = ns

     ! searches for the inverse of star member
     do i=iprot+1,ng
        if ((pbc%kgv(1,i)+pbc%kgv(1,iprot) == 0) .and. &
             (pbc%kgv(2,i)+pbc%kgv(2,iprot) == 0) .and. &
             (pbc%kgv(3,i)+pbc%kgv(3,iprot) == 0)) then
           if(i /= iprot+1) then
              kk(:) = pbc%kgv(:,iprot+1)
              pbc%kgv(:,iprot+1) = pbc%kgv(:,i)
              pbc%kgv(:,i) = kk(:)
              gmr = gm(iprot+1)
              gm(iprot+1) = gm(i)
              gm(i) = gmr
           endif
           !
           ! verifies if the inverse is related by symmetry
           !
           ncount = 0
           phtmp = zzero
           do j=1,symm%ntrans
              if ((ktran(1,j)-pbc%kgv(1,iprot+1) == 0) .and. &
                   (ktran(2,j)-pbc%kgv(2,iprot+1) == 0) .and. &
                   (ktran(3,j)-pbc%kgv(3,iprot+1) == 0)) then
                 phtmp = phtmp + expt(j)
                 ncount = ncount + 1
              endif
           enddo

           if(ncount /= 0) then
              pbc%phase(iprot+1) = phtmp/real(ncount,dp)
              pbc%conj(iprot+1) = one
           else
              pbc%phase(iprot+1) = pbc%phase(iprot)
              pbc%conj(iprot+1) = -one
           endif

           tmstar(ns) = tmstar(ns) + 1
           pbc%inds(iprot+1) = ns

           goto 31

        endif
     enddo

31   continue

     iavail = iprot + 2
     itest = iprot + 2

  else
     !
     ! checks if it is related by symmetry
     !
     ncount = 0
     phtmp = zzero
     do i=1,symm%ntrans
        if ((ktran(1,i)-pbc%kgv(1,itest) == 0) .and. &
             (ktran(2,i)-pbc%kgv(2,itest) == 0) .and. &
             (ktran(3,i)-pbc%kgv(3,itest) == 0)) then
           phtmp = phtmp + expt(i)
           ncount = ncount + 1
        endif
     enddo

     if(ncount /= 0) then
        !
        ! it is related by symmetry
        !
        if(itest /= iavail) then
           kk(:) = pbc%kgv(:,iavail)
           pbc%kgv(:,iavail) = pbc%kgv(:,itest)
           pbc%kgv(:,itest) = kk(:)
           gmr = gm(iavail)
           gm(iavail) = gm(itest)
           gm(itest) = gmr
        endif
        pbc%phase(iavail) = phtmp/real(ncount,dp)
        pbc%conj(iavail) = one
        tmstar(ns) = tmstar(ns) + 1
        pbc%inds(iavail) = ns
        !
        ! searches for the inverse of star member
        !
        do i=iavail+1,ng
           if ((pbc%kgv(1,i)+pbc%kgv(1,iavail) == 0) .and. &
                (pbc%kgv(2,i)+pbc%kgv(2,iavail) == 0) .and. &
                (pbc%kgv(3,i)+pbc%kgv(3,iavail) == 0)) then
              if(i /= iavail+1) then
                 kk(:) = pbc%kgv(:,iavail+1)
                 pbc%kgv(:,iavail+1) = pbc%kgv(:,i)
                 pbc%kgv(:,i) = kk(:)
                 gmr = gm(iavail+1)
                 gm(iavail+1) = gm(i)
                 gm(i) = gmr
              endif
              !
              ! verifies if the inverse is related by symmetry
              !
              ncount = 0
              phtmp = zzero
              do j=1,symm%ntrans
                 if ((ktran(1,j)-pbc%kgv(1,iavail+1) == 0) .and. &
                      (ktran(2,j)-pbc%kgv(2,iavail+1) == 0) .and. &
                      (ktran(3,j)-pbc%kgv(3,iavail+1) == 0)) then
                    phtmp = phtmp + expt(j)
                    ncount = ncount + 1
                 endif
              enddo

              if(ncount /= 0) then
                 pbc%phase(iavail+1) = phtmp/real(ncount,dp)
                 pbc%conj(iavail+1) = one
              else
                 pbc%phase(iavail+1) = pbc%phase(iavail)
                 pbc%conj(iavail+1) = -one
              endif
              tmstar(ns) = tmstar(ns) + 1
              pbc%inds(iavail+1) = ns

              goto 32

           endif
        enddo

32      continue

        if(itest == iavail) then
           itest = itest +2
        else
           itest = itest +1
        endif
        iavail = iavail +2

     else
        !
        ! increases test counter
        !
        itest = itest + 1
     endif

  endif

  if(iavail <= ng) goto 30
  !
  ! reset kmax 
  !
  kmax = maxval(pbc%kgv,dim=2)
  !
  ! checks stars with zero phase factors
  !
  if (.not. allocated (izstar)) allocate (izstar (1:ns))
  ncount = 1
  izstar(1) = 1
  do i=2,ns
     if(abs(pbc%phase(ncount+1)-one) > delta) then
        do j=ncount+1,ncount+tmstar(i)
           pbc%phase(j) = zzero
        enddo
        izstar(i) = 0
     else
        izstar(i) = 1
     endif
     ncount = ncount + tmstar(i)
  enddo
  !
  ! Save information of stars
  !
  call pbc_set_nstar (pbc,ns)
  pbc%mstar(1:pbc%nstar) = tmstar(1:pbc%nstar)
  call dcopy(pbc%nstar,tek,1,pbc%ek,1)

  write(7,50) ng,ns,(kmax(i),i=1,3)
  if (ipr >= 5) then
     !
     ! print g-vectors
     !
     istop = 0
     do i=1,pbc%nstar
        write(7,52) i,pbc%mstar(i),pbc%ek(i),izstar(i)
        istar = istop+1
        istop = istar+pbc%mstar(i)-1
        do j=istar,istop
           write(7,54) j,pbc%inds(j),(pbc%kgv(k,j),k=1,3), &
                real(pbc%phase(j)),aimag(pbc%phase(j)),pbc%conj(j)
        enddo
     enddo
  endif

  deallocate(gm)
  deallocate(izstar)
  deallocate(tmstar)
  deallocate(tek)

50 format(/,1x,i8,' g-vectors are set up in ',i7,' stars', &
        ' - kmax = ',3i4)
52 format(/,' star no',i4,'  with',i3,' members - kin energy =', &
        f12.8,4x,i3,//, &
        4x,'i',3x,'inds',5x,'kx',3x,'ky',3x,'kz',9x,'phase')
54 format(i5,2x,i5,2x,3i5,f10.5,3x,f10.5,5x,f5.2)
56 format('  *** stopped in gspace   no. of g-vec.= ',i6)

end subroutine g_space
!===============================================================
! M. Jain:
! The remaining routines here are needed for the BerkeleyGW 
! interface. If the compilation of these takes too much time
! then it is easy to wrap them in ifdef BGW flags.
!===============================================================
subroutine generate_gspace(gs, pbc, symms)  

  use constants
  use pbc_module
  use symmetry_module
  use gspace_module

  implicit none
  !

  !      INPUT:
  !      ------
  !
 !integer, intent(in) :: ioutput                               ! output flag
  type(pbc_data), intent(in) :: pbc
  type(symmetry), intent(in) :: symms  
  !
  !      INPUT/OUTPUT:
  !      -------------
  !
  type(gspace), intent(inout) :: gs  
  !
  !
  !     The routine to collect the gvectors into stars was written by JLM.
  !
  !     All the other stuff was written by  Bernd Pfrommer, UC Berkeley
  !
  !
  !
  !     ------------------------- actions --------------------------------
  !
  !     *  generates the gspace according to the variable
  !                gs%gmax
  !                gs%rk(3)  (the k-point)
  !     *  if
  !           gs%igvec    generate the gvectors (h,k,l)
  !
  !
  !     this subroutine allocates memory, and attaches it to various
  !     structures:
  !     *  if
  !
  !           gs%igvec    ==>  gs%gvec
  !
  !     ------------------ local variables ------------------------------
  !
  real(dp) :: qk(3), gmod, gmax2, rgm, rkmod, fi, gmr, diff, t0
  complex(dp) :: ph, ex(48), str
  integer :: &
       ns, &            ! number of stars
       i, j, k, &       ! g-space integers
       l, ir, ih, &
       ncount, &        ! counts number of symmops leaving G invariant
       ktran(3, 48), &  ! symmetry transformed vector
       ngt, &           ! number of vectors in total gspace
       ng, &            ! number of vectors on myproc
       myblocks, kstart, nrods, &
       ntotrod, &       ! total number of rods everywhere
       ib, &            ! counts blocks
       ip, &            ! counts planes within one block
       iprot, iavail, itest, &
       kk(3), kmax(3), nx, ny, p, iproc, lorder, nt, ja, jmax, &
       igv(3), ili, &
       sflag            ! detect if inside sphere
  !
  !     variables for the generic gspace loop
  !
  integer :: il, fftn(4), ffth(4), iorder, irod, gv(4), gv3,icomp1  
  real(dp), parameter :: delta = 1.0d-6
  integer, allocatable :: gtvec(:,:) 
  real(dp), allocatable :: ekt(:)  
  logical :: iflag
  !!integer, external :: findvec  !never used, why kept it here? (it failed the g95 compiler)
  !real(dp), external :: gimmetime  
  !     -----------------------------------------------------------------
  !     calculate length of k-vector

!  t0 = gimmetime()  

  rkmod = sqrt(dot_product(gs%rk, matmul(pbc%bdot, gs%rk)))  
  !     determine boundaries kmax such that grid holds the sphere snugly
  call findbound(gs%gmax + rkmod, pbc%bvec, kmax(1), kmax(2), kmax(3))
  gs%kmax(:) = kmax(:)
  !
  !
  !     Pick nice FFT grid if it has not already been picked
  !
  !
  if (gs%fftsize(1) <= 0 .or. gs%fftsize(2) <= 0 .or. gs%fftsize(3) <= 0) then
     gs%fftsize(1:3) = 2 * kmax(1:3) + 1
     call adjustfft (gs%fftsize, symms)
  end if
  gs%ig0 = -1                         ! will be set to index of G=0
  gmax2 = gs%gmax * gs%gmax  
  !
  !     ------- count the number of elements in total and small gspace. --
  !
  gs%length = 0  
  do i = -kmax(1), kmax(1)  
     qk(1) = real(i, dp) + gs%rk(1)  
     do j = -kmax(2), kmax(2)  
        qk(2) = real(j,dp) + gs%rk(2)  
        do k = -kmax(3), kmax(3)  
           qk(3) = real(k, dp) + gs%rk(3)  
           !              calculate norm with special metric
           gmod = (qk(1) * pbc%bdot(1, 1) + qk(2) * pbc%bdot(2, 1) + &
                qk(3) * pbc%bdot(3, 1)) * qk(1) + &
                (qk(1) * pbc%bdot(1, 2) + qk(2) * pbc%bdot(2, 2) + &
                qk(3) * pbc%bdot(3, 2)) * qk(2) + &
                (qk(1) * pbc%bdot(1, 3) + qk(2) * pbc%bdot(2, 3) + &
                qk(3) * pbc%bdot(3, 3)) * qk(3)
           if (gmod <= gmax2) then                  ! we are inside the sphere
              gs%length = gs%length + 1
           end if
        end do                                      ! loop along 3rd dimension
     end do      ! 2nd dim
  end do         ! 1st dim
  !
  !     allocate arrays
  !
  ! array for index mapping
  if (gs%imap) then
     allocate(gs%indexg(-kmax(1):kmax(1),-kmax(2):kmax(2),-kmax(3):kmax(3)))
     gs%indexg = 0
  else
     nullify(gs%indexg)
  endif
  if (gs%igvec) then
     allocate(gs%gvec(3, gs%length))  
  else
     nullify(gs%gvec)
  end if
  !
  !
  !     ---------------- now really generate the gspace--------------
  !
  ng = 0  
  ngt = 0  
  do i = -kmax(1), kmax(1)  
     qk(1) = real(i, dp) + gs%rk(1)  
     do j = -kmax(2), kmax(2)  
        qk(2) = real(j, dp) + gs%rk(2)  
        do k = -kmax(3), kmax(3)  
           qk(3) = real(k, dp) + gs%rk(3)  
           !              calculate norm with special metric
           gmod = (qk(1) * pbc%bdot(1, 1) + qk(2) * pbc%bdot(2, 1) + &
                qk(3) * pbc%bdot(3, 1)) * qk(1) + &
                (qk(1) * pbc%bdot(1, 2) + qk(2) * pbc%bdot(2, 2) + &
                qk(3) * pbc%bdot(3, 2)) * qk(2) + &
                (qk(1) * pbc%bdot(1, 3) + qk(2) * pbc%bdot(2, 3) + &
                qk(3) * pbc%bdot(3, 3)) * qk(3)
           if (gmod <= gmax2) then  !     -------- we are inside the sphere
              ng = ng + 1  
              if (i == 0 .and. j == 0 .and. k == 0) gs%ig0 = ng  
                 ! i hold the G=0 vector
              if (gs%imap) gs%indexg(i,j,k) = ng
              if (gs%igvec) then  
                 gs%gvec(1, ng) = i  
                 gs%gvec(2, ng) = j  
                 gs%gvec(3, ng) = k  
                 !                     write(9,*) ng, nx,ny,gs%gvec(:,ng)
              end if
           end if
        end do    ! loop along 3rd dimension
     end do     ! 2nd dim
  end do        ! 1st dim

  !call myflush(9)  
  !call myflush(7)  

  return  

940 format(' TIME FOR GSPACE SETUP:',f12.3)  

end subroutine generate_gspace


!
subroutine adjustfft(nfft, symms)

  use symmetry_module
  implicit none           

  !
  !     1996 Bernd Pfrommer
  !     
  !     

  !
  !     INPUT:
  !     -----

  type(symmetry) :: symms       ! the symmetry operations of the crystal
  !
  !     INPUT/OUTPUT:
  !     ------------

  integer :: nfft(3)            ! fft grid dimensions
  !
  !     DESCRIPTION
  !     -----------
  !
  !     Increases the fft grid in nfft(3) to fit a suitable value, e.g.
  !     a power of two, three ...
  !
  !     Also checks to make sure the resulting grid is compatible with
  !     the nonprimitive translations, e.g. if there is a translation
  !     of 1/3, the grid must be a multiple of 3.
  !
  !
  !     ------------------- local arrays ---------------------------

  integer :: i, j, k, ns, sym_k
  real(dp) :: tau
  integer :: fracprim(4), fracfac(3)

  logical :: ifit, fracpresent(4)




  ! M. Jain:
  ! Find the Lowest common multiple of all the integers
  ! that will ensure that the fractional translations will
  ! take a point on the FFT grid to another on the FFT grid

  do i = 1, 3
     fracprim = (/ 1, 2, 3, 4 /) 
     fracpresent = (/ .false., .false., .false., .false. /) 
     do sym_k = 1, symms%ntrans
        if (abs(symms%tnp(i,sym_k)) .gt. 1.0d-5) then
           do j = 1, 4
               if (abs(abs(nint(symms%tnp(i,sym_k)*fracprim(j))) - one) &
                   .lt. 1.0d-5) then
                  fracpresent(j) = .true.
               endif
           enddo
        endif
     enddo
     fracfac(i) = 1
     do j = 1, 4
        if (fracpresent(j)) then
           fracfac(i) = fracfac(i)*fracprim(j)
        endif
     enddo
     if (fracpresent(2) .and. fracpresent(4)) then
         fracfac(i) = fracfac(i)/2
     endif
  enddo

  

  do i=1,3
    ! redundant check
    ifit = .false.
    do while( .not.(ifit) )
      nfft(i) = fastnum( nfft(i) )
      ! check if nfft(i) is compatible with symmetry
      if (mod(nfft(i),fracfac(i)) .eq. 0) ifit=.true. 
      ! if not increase nfft(i)
      if( .not.(ifit) ) nfft(i)=nfft(i)+1
    enddo
  enddo
 

  return

110 format(' *** ERROR IN ADJUSTFFT: ', &
       'FFT GRIDSIZE EXCEEDS IN DIRECTION',i3, &
       /' FFT GRIDSIZE = ', 3i6, &
       /' MAX GRIDSIZE = ', i6, &
       /'     THIS CAN BE DUE TO STRANGE NONPRIMITIVE', &
       /'     TRANSLATIONS IN THE SYMMETRY OPERATIONS', &
       /'     OR JUST BECAUSE THE SYSTEM IS TOO BIG.', &
       /'     ADD LARGER FFT SIZES IN ADJUSTFFT OR', &
       /'     IMPROVE THE COORDINATES/LATTICE VECTORS.')


  contains

  function fastnum( nfft_in )
  integer :: fastnum
  integer,intent(in) :: nfft_in

  integer :: nfft_out
  integer :: n, i
  !integer :: prim(6), primpower(6)
  integer :: prim(3), primpower(6)
  logical :: notfound

!  - best grids are given by : 2^a*3^b*5^c*7^d*11^e*13^f
!        a,b,c,d are arbitary, e,f are 0,1

  ! allowed prime factors
  !prim = (/ 2, 3, 5, 7, 11, 13 /)
  prim = (/ 2, 3, 5/)

  nfft_out=nfft_in
  notfound=.true.
  do while( notfound )

    notfound=.false.
    n=nfft_out
    primpower=0

    ! generate prime power expansion
    do i=1,size(prim)

      do while( mod(n,prim(i))==0 )

        n=n/prim(i)
        primpower(i)=primpower(i)+1

      enddo

    enddo

    n=1
    do i=1,size(prim)
      n=n*prim(i)**primpower(i)
    enddo

    !if( n/=nfft_out .or. primpower(5) > 1 .or. primpower(6) > 1 ) then
    if( n/=nfft_out) then
      notfound=.true.
      nfft_out=nfft_out+1
    endif

  enddo

  ! final result
  fastnum = nfft_out

  end function fastnum


end subroutine adjustfft


!
subroutine generate_k_gspaces(k_gspace_irk, maxlength, irk, fftsize, &
     pbc, kpt, kpwt, syms, emax)

  use constants
  use electronic_struct_module
  use pbc_module
  use symmetry_module
  use gspace_module

  implicit none
  !
  !
  !
  !     INPUT:
  !     -----
  !
  real(dp), intent(in) :: kpt(3) ! kpoint in reciprocal lattice units
  real(dp), intent(in) :: kpwt ! weight of this kpoint 
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  type(symmetry), intent(in) :: syms  
  real(dp), intent(in) :: emax          ! energy cutoff for the wave functions!
  real(dp), intent(in) :: fftsize(1:3)
  integer, intent(in) :: irk
  !
  !     OUTPUT:
  !     ------
  !
  type(gspace), intent(inout) :: k_gspace_irk
  integer, intent(inout) :: maxlength   ! maximum length of all gspaces generated
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     generates the small gspaces for the various kpoints
  !
  !
  !     --------  local variables ---------------------------------
  !
  integer :: mtrxdim, i  

  if (irk == 1) then
     write(7, 10)  
     maxlength = 0  
  endif
  k_gspace_irk%gmax = sqrt(emax)  
  k_gspace_irk%rk(:) = kpt(:)
  k_gspace_irk%name = 'k_gspace'  
  k_gspace_irk%fftsize(:) = fftsize(:)
  call generate_gspace(k_gspace_irk, pbc, syms)  
  if (k_gspace_irk%length > maxlength) maxlength = k_gspace_irk%length
  mtrxdim = k_gspace_irk%length  
  write(7, 20) irk, kpwt,(kpt(i), i = 1, 3), mtrxdim  

  call myflush(7)
  
10 format(/' HAMILTONIAN MATRIX DIMENSIONS:',/1x,27('-'), &
       &           /12x,'WEIGHT', 4x,'K-POINT',30x,'DIMENSION')

20 format(i4,4x,4f10.4,12x,i6)  

end subroutine generate_k_gspaces


!
subroutine generate_potential_gspace(pot_gspace, pbc, syms, emax)
  
  use constants
  use pbc_module
  use symmetry_module
  use gspace_module
  !
  implicit none   
  !
  !     INPUT:
  !     -----
  !
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  type(symmetry), intent(in) :: syms  
  real(dp), intent(in) :: emax                   ! energy cutoff for the wavefunctions!
  !
  !     OUTPUT:
  !     ------
  !
  type(gspace), intent(inout) :: pot_gspace    ! potential gspace
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     generates a large potential gspace with all the features:
  !     stars, gvectors, kinetic energy, structure factors
  !
  !
  !     ------------------------------------------------------------------
  !
  pot_gspace%gmax = two * sqrt(emax)          ! 4*emax
  pot_gspace%rk = (/ zero, zero, zero /)    ! no shift
  pot_gspace%name = 'potential'  
  pot_gspace%fftsize = (/ 0, 0, 0 /)

  call generate_gspace(pot_gspace, pbc, syms)  
  write(7, 10) pot_gspace%length, pot_gspace%fftsize(1:3)  
  call myflush(7)  

10 format(' SIZE OF POTENTIAL GSPACE = ', i8, &
       &      ' WITH FFT GRID ', 3i6)

end subroutine generate_potential_gspace


!
!=======================================================================
!
function scp(v1, v2)

  use constants
  implicit none

  real(dp) :: scp                    !p    scalar product without metric
  real(dp), intent(in) :: v1(3), v2(3)  

  scp = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
  
  return  

end function scp
!
!=======================================================================
!

subroutine findbound(qmax, b, maxi, maxj, maxk)  
  !p
  !p    adapted from subroutine ban 3.2.93, Bernd Pfrommer
  !p
  !p    ! bi,bj,bk should form a right-handed system !
  !p
  !p    subroutine findbound finds boundaries of the space built up by
  !p    the vectors bi,bj,bk if a sphere of
  !p    radius qmax should be placed inside that space. In maxi, maxj, max
  !p    integers values are returned representing the coordinates along bi
  !p    that are necessary to involve the qmax-sphere entirely.
  !p
  !     -------------------- documentation of variables ------------------
  !p    qmax      radius of involved sphere
  !p    b(i,k)    ith component of recip. lattice vector k.
  !p    maxi,maxj,maxk integer variables used to return required minimum
  !p              boundaries.
  !p    volp      volume produkt bi*(bj x bk)
  !p    proi,proj,prok  see documentation in source
  !p    bj_times_bk, bk_times_bj, bi_times_bk  absolute value of cross product
  !p

  use constants
  implicit none  
  real(dp), intent(in) :: qmax                ! radius of involved sphere
  integer, intent(out) :: maxi, maxj, maxk    ! return values: boundaries
  real(dp) :: scp                             ! standard scalar product
  real(dp), intent(in) :: b(3, 3)             ! reciprocal lattice vectors
  !
  !p    ------- local variables -------------------
  !
  real(dp) :: bi(3), bj(3), bk(3)      ! rec lattice vectors (or any vector)
  integer :: i  
  real(dp) :: volp                     ! volume product
  real(dp) :: proi, proj, prok         ! see below

  real(dp) :: bj_times_bk, bk_times_bi, bi_times_bj

! abs. value of cross prod.  
  do i = 1, 3  
     bi(i) = b(i, 1)  
     bj(i) = b(i, 2)  
     bk(i) = b(i, 3)  
  end do

  !p    volume product bi*(bj x bk)

  volp = bi(1) * (bj(2) * bk(3) - bj(3) * bk(2)) + &
       bi(2) * (bj(3) * bk(1) - bj(1) * bk(3)) + &
       bi(3) * (bj(1) * bk(2) - bj(2) * bk(1))
  !
  !p    project bi on normal vector to the plane (bk,bj) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bk and bj.
  !p    proi=bi*(bj x bk)/|bj x bk|
  !
  bj_times_bk = sqrt(scp(bj, bj) * scp(bk, bk) - scp(bj, bk)**2)
  if (bj_times_bk <= 1.0d-6) then  
     write (7, *) 'error in findbound: bj,bk parallel'  
     stop
  else  
     proi = abs(volp / bj_times_bk)  
     maxi = int(qmax / proi)  
  end if
  !
  !p    project bj on normal vector to the plane (bk,bi) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bk and bi.
  !p    proj=bj*(bk x bi)/|bk x bi|
  !
  bk_times_bi = sqrt(scp(bk, bk) * scp(bi, bi) - scp(bk, bi)**2)
  if (bk_times_bi <= 1.0d-6) then  
     write(7, *) 'error in findbound: bk,bi parallel'  
     stop
  else  
     proj = abs(volp / bk_times_bi)  
     maxj = int(qmax / proj)  
  end if
  !
  !p    project bk on normal vector to the plane (bi,bj) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bi and bj.
  !p    prok=bk*(bi x bj)/|bi x bj|
  !
  bi_times_bj = sqrt(scp(bi, bi) * scp(bj, bj) - scp(bi, bj)**2)
  if (bi_times_bj <= 1.0d-6) then  
     write(7, *) 'error in findbound: bi,bj parallel'  
     stop
  else  
     prok = abs(volp / bi_times_bj)  
     maxk = int(qmax / prok)  
  end if

end subroutine findbound
