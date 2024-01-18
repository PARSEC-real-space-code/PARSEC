c     ===============================================================
c
c     Copyright (C) 2005 Finite Difference Research Group
c     This file is part of parsec, http://www.ices.utexas.edu/parsec/
c
c     Reads file parsec.dat and produces output to be visualizes with
c     OpenDX (or DataExplorer). Output functions are: self-consistent
c     potential, density, or wave-functions.
c
c     input file : parsec.dat
c     output file: [ out_file ]
c
c     Murilo L Tiago, UTexas-Austin, August 2005, mtiago@ices.utexas.edu
c
c     ---------------------------------------------------------------
c
c     Constants
c
      module constants
      integer, parameter :: dp = kind(1.0d0)
      integer, parameter :: dpc = kind((1.0d0,1.0d0))
      end module constants
c
c     ---------------------------------------------------------------
      program plotdx

      use constants
      implicit none

      integer, parameter :: nmax = 100

      character (len=200) :: outfile
      character (len=26) :: label
      integer :: ipr, ispflag, wfnflag,neig, pbcflag, ncplx
      integer :: ndim, ndim_full, ntrans, nspin, nstate, nkpt
      integer, dimension(nmax) :: indxeig, indxspn, indxkpt
      integer :: ii, jj, kk, is, ikp, jr, nw(3), rmin(3), rmax(3)
      real(dp) :: hh(3), rsize, shift(3), vec(3), spin
      real(dp) :: avec(3,3), anorm(3,3), mtmp(3,3), dvol, rtmp
      integer, allocatable :: rgrid(:,:), rgrid_full(:,:),
     1     irep(:,:,:), chi(:,:), istate(:,:), indxprint(:),
     2     indxfound(:), imtrx(:,:,:)
      real(dp), allocatable :: wfn(:),trans(:,:,:),
     1     plotbox(:,:,:), kw(:), kpt(:,:)
      complex(dpc), allocatable :: zwfn(:)
c
      open(20,file='parsec.dat',form='unformatted',status='old')
      ipr = 0
      ispflag = 0
      wfnflag = 0
      pbcflag = 0

      write(6,*) 'Output file ? '
      read(5,'(a)') outfile
      open(30,file=trim(outfile),form='formatted')
      write(6,*) 'Reading input file parsec.dat'
      write(6,*) 'writing output to file ',trim(outfile)

      write(6,*) 'Plot density (1), potential (2), or ',
     1     'wave-functions (3) ? '
      read(5,*) ipr
      select case(ipr)
      case (1)
         write(6,*) 'Printing electron density [ 1/bohr^3 ]'
      case (2)
         write(6,*) 'Printing self-consistent potential [Ry]'
      case (3)
         write(6,*) 'Printing wave-functions '
         write(6,*)
     1        'Print real part (1), imaginary part(2) or '
     2        ,'absolute value squared (3) ?'
         read(5,*) wfnflag
      end select

      read(20) label
      read(20) nspin, ncplx, pbcflag
      write(6,*) ' parsec.dat file created on ',label
      if (pbcflag .eq. 1) then
         read(20) hh(:)
         read(20) avec(:,:)
         do ii = 1, 3
            read(20)
         enddo
         read(20) nkpt
         do ii = 1, 3
            read(20)
         enddo
         if (nkpt .ne. 0) then
            allocate(kw(nkpt))
            allocate(kpt(3,nkpt))
            read(20) ((kpt(ii,jj),ii=1,3),jj=1,nkpt)
            read(20) (kw(ii), ii=1,nkpt)
            write(6,*) ' K-point coordinates and weights : '
            do ii = 1, nkpt
               write(6,'(3f10.4,3x,f10.4)') kpt(:,ii),kw(ii)
            enddo
         else
            read(20)
            read(20)
         endif
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
      if (nkpt .eq. 0) then
         nkpt = 1
         allocate(kw(nkpt))
         kw = 1.d0
         allocate(kpt(3,nkpt))
         kpt = 0.d0
      endif

      if (ncplx .eq. 1) then
          write(6,*) ' wave-functions are complex'
      else
          write(6,*) ' wave-functions are real'
      endif

      ispflag = 1
      if (ipr .eq. 1 .and. nspin .eq. 2) then
         write(6,*) ' System is spin polarized'
         write(6,*) ' Plot up component(1),',
     1        ' down component (2), total charge density (3) or '
         write(6,*) 'spin density (4)?'
         read(5,*) ispflag
      endif
      if (ipr .eq. 2 .and. nspin .eq. 2) then
         write(6,*) ' System is spin polarized'
         write(6,*) ' Plot up component(1) or', 
     1        ' down component (2) of potential?'
         read(5,*) ispflag
      endif

      if (ipr .eq. 3 .and. nspin .eq. 1) then
         if (nkpt .gt. 1) then
            write(6,*)
     1           ' Which eigenstates/k-points will be printed out ?'
            write(6,*) ' End list with 0 0. Output is added up '
         else
            write(6,*) ' Which eigenstates will be printed out ?'
            write(6,*) ' End list with 0. Output is added up '
         endif
         neig = 0
         do
            neig = neig + 1
            if (neig .gt. nmax) then
               write(6,*)
     1              'ERROR: too many eigenstates. Increase nmax'
               write(6,*) ' STOP'
               stop
            endif
            if (nkpt .gt. 1) then
               read(5,*) indxeig(neig), indxkpt(neig)
            else
               read(5,*) indxeig(neig)
               indxkpt(neig) = 1
            endif
            indxspn(neig) = 1
            if (indxeig(neig) .eq. 0) exit
            do ii = 1, neig-1
               if ((indxeig(ii) .eq. indxeig(neig)) .and. 
     1              (indxkpt(ii) .eq. indxkpt(neig))) then
                  neig = neig - 1
                  exit
               endif
            enddo
         enddo
         neig = neig - 1
      endif
      if (ipr .eq. 3 .and. nspin .eq. 2) then
         if (nkpt .gt. 1) then
            write(6,*)
     1           ' Which eigenstates/k-points/spins '
     2           ,'will be printed out ?'
            write(6,*) ' End list with 0 0 0. Output is added up '
         else
            write(6,*) ' Which eigenstates/spins will '
     1           ,'be printed out ?'
            write(6,*) ' End list with 0 0. Output is added up '
         endif
         neig = 0
         do
            neig = neig + 1
            if (neig .gt. nmax) then
               write(6,*)
     1              'ERROR: too many eigenstates. Increase nmax'
               write(6,*) ' STOP'
               stop
            endif
            if (nkpt .gt. 1) then
               read(5,*) indxeig(neig), indxspn(neig), indxkpt(neig)
            else
               read(5,*) indxeig(neig), indxspn(neig)
               indxkpt(neig) = 1
            endif
            do ii = 1, neig-1
               if ((indxeig(ii) .eq. indxeig(neig)) .and. 
     1              (indxkpt(ii) .eq. indxkpt(neig)) .and. 
     2              (indxspn(ii) .eq. indxspn(neig))) then
                  neig = neig - 1
                  exit
               endif
            enddo
            if (indxeig(neig) .eq. 0) exit
         enddo
         neig = neig - 1
      endif

c     Build full real-space grid. Assume tnp = 0.d0.
      read(20)
      read(20) shift
      read(20) ndim, ntrans
      allocate(imtrx(3,3,ntrans))
      read(20) imtrx
      allocate(trans(3,3,ntrans))
      read(20) trans
      read(20)
      read(20)
      allocate(chi(ntrans,ntrans))
      read(20) ((chi(ii,jj),ii=1,ntrans),jj=1,ntrans)
      ndim_full = ndim * ntrans

      allocate(rgrid(3,ndim))
      allocate(rgrid_full(3,ndim_full))
      read(20) ( (rgrid(ii,jr),ii=1,3), jr=1,ndim )
      call unfoldgrid(rgrid,rgrid_full,ndim,ntrans,shift,trans)
      rmin = minval(rgrid_full,2)
      rmax = maxval(rgrid_full,2)
      nw = rmax - rmin + 1
      deallocate(trans)

      call dx_latticevec(avec)
c
c     Renormalize lattice vectors
c
      anorm = avec
      do ii = 1, 3
         rtmp = sqrt( dot_product(anorm(:,ii),anorm(:,ii)) )
         anorm(:,ii) = anorm(:,ii)/rtmp
      enddo
      mtmp = anorm
      call mtrxin(mtmp,dvol,rtmp)

      dvol = hh(1) * hh(2) * hh(3) * dvol

      allocate(plotbox(nw(1),nw(2),nw(3)))
      plotbox = 0.0
      allocate(wfn(ndim))
c
      do is = 1, nspin
         do ikp = 1, nkpt
            read(20) nstate
            if (is .eq. 1 .and. ikp .eq. 1) then
               write(6,*) nstate,' eigenstates in file'
               allocate(irep(nstate,nkpt,nspin))
               allocate(istate(nkpt,nspin))
               allocate(indxprint(nstate))
            endif

            read(20) (irep(jj,ikp,is), jj=1,nstate)
            read(20)
            read(20)
         enddo

c     potential, trivial representation
         chi = 1
         read(20) (wfn(jr), jr=1,ndim)
         if (ipr .eq. 2 .and. is .eq. ispflag) call
     1        addbox(ndim,ntrans,nw,rmin,imtrx,
     2        rgrid,chi(1,1),shift,wfn,plotbox)
c     charge density, trivial representation
         read(20) (wfn(jr), jr=1,ndim)
         if (ipr .eq. 1 .and. ispflag .eq. is) call
     1        addbox(ndim,ntrans,nw,rmin,imtrx,
     2        rgrid,chi(1,1),shift,wfn,plotbox)
         if (ipr .eq. 1 .and. ispflag .eq. 3) call
     1        addbox(ndim,ntrans,nw,rmin,imtrx,
     2        rgrid,chi(1,1),shift,wfn,plotbox)
         spin = 3.d0 - 2.d0*is
         if (ipr .eq. 1 .and. ispflag .eq. 4) call
     1      addbox(ndim,ntrans,nw,rmin,imtrx,
     2        rgrid,chi(1,1),shift,spin*wfn,plotbox)

      enddo
      deallocate(chi,imtrx,rgrid)
      if (ipr .ne. 3) goto 90

      deallocate(wfn)
      allocate(indxfound(neig))
      indxfound = 0

      do is = 1, nspin
         do ikp = 1, nkpt
            read(20) ndim, ntrans
            allocate(imtrx(3,3,ntrans))
            read(20) imtrx
            read(20)
            read(20)
            read(20)
            allocate(chi(ntrans,ntrans))
            read(20) ((chi(ii,jj),ii=1,ntrans),jj=1,ntrans)
            allocate(rgrid(3,ndim))
            read(20) ( (rgrid(ii,jr),ii=1,3), jr=1,ndim )
            read(20) istate(ikp,is)
            if (istate(ikp,is) .eq. 0) then
               deallocate(imtrx)
               deallocate(chi)
               deallocate(rgrid)
               cycle
            endif
            read(20) (indxprint(ii),ii=1,istate(ikp,is))

            allocate(wfn(ndim))
            if (ncplx .ne. 0) allocate(zwfn(ndim))

            if (wfnflag .eq. 3) irep = 1

c     wave-functions, must add characters from the corresponding representation
            do jj = 1, istate(ikp,is)
               if (ncplx .eq. 0) then
                  read(20) (wfn(jr), jr=1,ndim)
               else
                  read(20) (zwfn(jr), jr=1,ndim)
               endif

               do ii = 1,neig
                  if (indxeig(ii) .eq. indxprint(jj) 
     1                 .and. indxspn(ii) .eq. is
     1                 .and. indxkpt(ii) .eq. ikp) then
                     indxfound(ii) = 1
                     if (ncplx .eq. 0) then
                        select case (wfnflag)
                        case(2)
                           wfn = 0.d0
                        case(3)
                           do jr = 1, ndim
                              wfn(jr) = wfn(jr) * wfn(jr)
                           enddo
                        end select
                     else
                        select case (wfnflag)
                        case(1)
                           wfn = real(zwfn,dp)
                        case(2)
                           wfn = imag(zwfn)
                        case(3)
                           do jr = 1, ndim
                              wfn(jr) = abs( zwfn(jr) )**2
                           enddo
                        end select
                     endif
                     call addbox(ndim,ntrans,nw,rmin,imtrx,rgrid,
     1                    chi(1:ntrans,irep(indxprint(jj),ikp,is)),
     2                    shift,wfn,plotbox)
                  endif
               enddo
            enddo               ! jj = 1, istate(ikp,is)
            deallocate(imtrx, chi, rgrid)
            deallocate(wfn)
            if (ncplx .eq. 1) deallocate(zwfn)
         enddo                  ! ikp = 1, nkpt
      enddo                     ! is = 1, nspin
      close(20)
      do ii = 1, neig
         if (indxfound(ii) .ne. 1) then
            write(6,*) ' WARNING!!!!! eigenvector # ',indxeig(ii)
     >           ,' at k-point ',indxkpt(ii),' spin ',indxspn(ii)
     >           ,' not found!'
         endif
      enddo

 90   continue

      write(6,*) ' Done! sum = ',sum(plotbox)*dvol

      write(30, 5)              ! write the header for the array
      write(30, 10) nw(1) * nw(2) * nw(3)

 5    format(/'# this is the object defining the data values')  
 10   format('object 1 class array type float rank 0 items ',i8,
     1     ' data follows')
c
      write(30,'(5g15.7)') (((plotbox(ii,jj,kk),ii=1,nw(1))
     1     ,jj=1,nw(2)),kk=1,nw(3))
c
      write(30, 20) nw(3), nw(2),  nw(1)
 20   format(/'# this is the object defining the grid connections',
     1     /'object 2 class gridconnections counts ',3i4)
                                !     write header for gridpositions
      write(30, 30) nw(3), nw(2), nw(1)
 30   format(//'# this is the object defining the grid ',
     1     /'object 3 class gridpositions counts ',3i4 )

      do ii = 1, 3
         vec(ii) = ( minval(rgrid_full(ii,:)) + shift(ii) ) * hh(ii)
      enddo
      vec = matmul(anorm,vec)
      write(30,'(a,3f20.10)') 'origin',vec

      do ii = 1, 4
         anorm(:,ii) = anorm(:,ii) * hh(ii)
      enddo
      write(30,50) anorm(:,3)
      write(30,50) anorm(:,2)
      write(30,50) anorm(:,1)
 50   format('delta ', 3f20.10)  
                                !     write header for field object
      write(30, 60)  
 60   format(/'# this is the collective object, one for each grid',
     1     /'object 4 class field',
     2     /'component "positions"   value 3',
     3     /'component "connections" value 2',
     4     /'component "data"        value 1')
      write(30,*)
      write(30,'(2a)') '# Date label : ',label
      if (ipr .eq. 1) then
         if (nspin .eq. 2 .and. ispflag .eq. 4) then
            write(30,'(a)') '# spin density'
         elseif (nspin .eq. 2 .and. ispflag .eq. 3) then
            write(30,'(a)') '# electron density, total'
         elseif (nspin .eq. 2 .and. ispflag .eq. 2) then
            write(30,'(a)') '# electron density, spin up component'
         elseif (nspin .eq. 2 .and. ispflag .eq. 1) then
            write(30,'(a)') '# electron density, spin down component'
         endif
      elseif (ipr .eq. 2) then
         write(30,'(a)') '# self-consistent potential'
         if (nspin .eq. 2 .and. ispflag .eq. 1) then
            write(30,'(a)') '# spin up component'
         elseif (nspin .eq. 2 .and. ispflag .eq. 2) then
            write(30,'(a)') '# spin down component'
         endif
      elseif (ipr .eq. 3) then
         write(30,'(a)') '# electron wave-functions'
         if (wfnflag .eq. 1) then
            write(30,'(a)') '# real part'
         elseif (wfnflag .eq. 2) then
            write(30,'(a)') '# real part'
         elseif (wfnflag .eq. 3) then
            write(30,'(a)') '# absolute value squared'
         endif
         if (nkpt .eq. 1) then
            if (nspin .eq. 1) then
               write(30,'(a)') '# printed eigenstates : '
               write(30,70) (indxeig(ii),ii=1,neig)
            elseif (nspin .eq. 2) then
               write(30,'(a)') '# printed eigenstates (ii,isp): '
               write(30,75) (indxeig(ii),indxspn(ii),ii=1,neig)
            endif
         else
            if (nspin .eq. 1) then
               write(30,'(a)') '# printed eigenstates (ii,ikp): '
               write(30,80) (indxeig(ii),indxkpt(ii),ii=1,neig)
            elseif (nspin .eq. 2) then
               write(30,'(a)') '# printed eigenstates (ii,ikp,isp): '
               write(30,85) (indxeig(ii),indxkpt(ii),indxspn(ii),ii=1
     1              ,neig)
            endif
         endif
      endif
 70   format('# ',10i4)
 75   format(10('# ',3('( ',i4,' , ',i1,' ) '),/))
 80   format(10('# ',3('( ',i4,' , ',i4,' ) '),/))
 85   format(10('# ',3('( ',i4,' , ',i4,' , ',i1,' ) '),/))

      close(30)

      end program plotdx
c     ===============================================================
      subroutine unfoldgrid(rgrid,rgrid_full,ndim,ntrans,shift,trans)

      use constants
      implicit none

      integer, intent(in) :: ndim,ntrans
      real(dp), intent(in) :: shift(3),trans(3,3,ntrans)
      integer, intent(in) :: rgrid(3,ndim)
      integer, intent(out) :: rgrid_full(3,ntrans*ndim)

      integer ir, jr, kt, ii
      real(dp) :: rvec(3), tvec(3)

      do ir = 1, ndim
         do kt = 1, ntrans
            jr = kt + (ir-1)*ntrans
            rvec = rgrid(:,ir) + shift
            do ii = 1, 3
               tvec(ii) = dot_product(rvec,trans(:,ii,kt))
            enddo
            tvec = tvec - shift
            rgrid_full(:,jr) = nint(tvec)
         enddo
      enddo

      end subroutine unfoldgrid
c     ===============================================================
      subroutine addbox(ndim,ntrans,nw,rmin,imtrx,rgrid,chi,shift,
     1     wfn,plotbox)

      use constants
      implicit none

      integer, intent(in) :: ndim,ntrans,nw(3),rmin(3)
      integer, intent(in) :: imtrx(3,3,ntrans),rgrid(3,ndim),
     1     chi(ntrans)
      real(dp), intent(in) :: shift(3),wfn(ndim)
      real(dp), intent(inout) :: plotbox(nw(1),nw(2),nw(3))

      integer ir, kt, i1, i2, i3, ii
      real(dp) :: wfn_tmp
      real(dp) :: rtmp(3),gridpt(3),rmtrx(3,3,ntrans)

      rmtrx = imtrx
      do ir = 1, ndim
         do kt = 1, ntrans
            gridpt = real(rgrid(:,ir),dp) + shift
            do ii = 1, 3
               rtmp(ii) = dot_product(rmtrx(ii,:,kt),gridpt)
            enddo
            gridpt = rtmp - shift
            i1 = nint(gridpt(1)) - rmin(1) + 1
            i2 = nint(gridpt(2)) - rmin(2) + 1
            i3 = nint(gridpt(3)) - rmin(3) + 1
            if (i1 .lt. 1 .or.i1 .gt. nw(1)) then
              write(6,*) ' ERROR i1 ',i1,ir,kt,rgrid(:,ir)
            endif
            if (i2 .lt. 1 .or.i2 .gt. nw(2)) then
              write(6,*) ' ERROR i2 ',i2,ir,kt,rgrid(:,ir)
            endif
            if (i3 .lt. 1 .or.i3 .gt. nw(3)) then
              write(6,*) ' ERROR i3 ',i3,ir,kt,rgrid(:,ir)
            endif
            wfn_tmp = wfn(ir) * real(chi(kt),dp)
            plotbox(i1,i2,i3) = plotbox(i1,i2,i3) + wfn_tmp
         enddo
      enddo

      end subroutine addbox
c     ===============================================================
      subroutine dx_latticevec(avec)
c
c     writes files LATTICE_VEC.dx and UCELL_FRAME.dx for the
c     IBM data explorer, a visualization software package.
c     1996 Bernd Pfrommer, UC Berkeley.
c
      use constants
      implicit none

      integer :: i, j
      real(dp) :: avec(3,3)
 
      open(19, file = 'LATTICE_VEC.dx', form= 'formatted')
      write(19, 10)  
      write(19, 20) 1, 3  
      write(19, 30) ((avec(i, j), i = 1, 3), j = 1, 3)
      write(19, 20) 2, 3  
      write(19, 35)  
      write(19, 40) 3, 1, 2  
      write(19, 50)  
      close(19)
  
      open(19, file= 'UCELL_FRAME.dx',form = 'formatted')
      write(19, 60)  
      write(19, 15) 3, 12  
      write(19,'(2i3)') 0, 1, 0, 2, 0, 3, 1, 4, 1, 5, 3, 5, 3, 6, 2,
     1     6, 2, 4, 7, 5, 7, 6, 7, 4
      write(19, *) 'attribute "element type" string "lines"'  
      write(19, 20) 4, 8  
      write(19, 30) 0.0,0.0,0.0
      write(19, 30) ((avec(i, j), i = 1, 3), j = 1, 3)
      write(19, 30) (avec(i, 1) + avec(i, 2), i = 1, 3)  
      write(19, 30) (avec(i, 1) + avec(i, 3), i = 1, 3)  
      write(19, 30) (avec(i, 2) + avec(i, 3), i = 1, 3)  
      write(19, 30) (avec(i, 1) + avec(i, 2) + avec(i, 3), i = 1, 3)
      write(19, 70) 5, 12  
      write(19, '(12(''1.0''/))')  
      write(19, *) 'attribute "dep" string "connections"'  
      write(19, 45) 6, 5, 4, 3  
      write(19, 50)  
      close(19)  

 10   format(2('#',/),'#    LATTICE VECTOR INFO:',2(/'#'))  
 15   format('object ',i2,' class array type '
     1     ,'int rank 1 shape 2 items ',i4,' data follows')
 20   format('object ',i2,' class array type '
     1     ,'float rank 1 shape 3 items ',i4,' data follows')
 30   format(3(f15.8))  
 35   format(3('0 0 0'/))  
 40   format('object ',i3,' class field',/'component "data" value '
     1     ,i2,/'component "positions" value ',i2)
 45   format('object ',i3,' class field',/'component "data" value '
     1     ,i2,/'component "positions" value ',i2,/'component
     2     "connections" value ',i2)

 50   format('end')  
 60   format(2('#',/),'#   UNIT CELL FRAME:',2(/'#'))  

 70   format('object ',i4,' array type float rank 0 items ',i4
     1     ,' data follows')

      end subroutine dx_latticevec
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

