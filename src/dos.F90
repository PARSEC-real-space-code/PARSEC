!===============================================================
!     1996 Bernd Pfrommer
!
!     Based on subroutines by Alberto Garcia, Paul Delaney, and
!     Michel Cote.
!
! CHANGES
! 
!  01/2007 imported to Parsec, and modified accordingly.
!   Or Cohen (changed to dos_graph and added dos)
!  *****************************
!  Notice that this function completely messes up with the arrays
!  of PARSEC, so it better be used last.
!  *****************************
!  03/2007 Added Ylm DOS generation
!---------------------------------------------------------------
  subroutine dos(clust,elec_st,grid,pbc,parallel,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pbc_module
  use parallel_data_module
  implicit none  


  !     INPUT:
  !     -----
  !
  !     the cluster
      type (cluster), intent(in) :: clust
  !     electronic structure
      type (electronic_struct), intent(inout) :: elec_st
  !     grid related data
      type (grid_data), intent(in) :: grid
  !     pbc related data
      type (pbc_data), intent(in) :: pbc
  !     parallel computation related data
      type (parallel_data), intent(in) :: parallel
  !     error flag
      integer, intent(out) :: ierr

      
  !
  !     The weights for different angular momenta, different bands,
  !     kpoints, and spins. A fancy interpolation between the function
  !     values at different kpoints is done inside the subroutine
  !     spec_tet.
  !
  !
  !   DESCRIPTION:
  !   -----------
  !
  !   computes DOS and angular momentum resolved DOS using a
  !   tetrahedron method. The data must be given on a Monkhorst-Pack grid
  !
  !   -------------------local variables ------------------------------
  !
  integer :: n, i, j, k, l, nx, ny, nz, isp, ie, it, iylm, iat, b, i1, &
       j1, k1, ibn, nmesh, nkpt, kidx(4), itet, ntetra, na, nmytet, ikp, &
       itetstart, itetend, neqtet, idum, mnatmi, irp, nn, nylm,  &
       atom_type(clust%atom_num), atom_num(clust%atom_num)
  real(dp) :: e(0:3), f(0:3), t0, t1, ecorn(4), fcorn(4), &
       vol, mineigen, maxeigen, delta, fav, t, ebound(2)
  real(dp), allocatable :: et(:,:,:,:), ylmdos(:,:,:,:), fulldos(:,:), dummy(:), &
  ylmwgt(:,:,:,:,:,:)
  integer, allocatable :: itetdef(:,:)  
  real(dp) :: tetrylmdos, dos_factor
  character(len=10),  dimension (16) :: ylmname=(/'s         ', &
  'px        ','py        ','pz        ','dxy       ', &
  'dzy       ','dxz       ','dz2       ','dx2-y2    ', & 
  'fz3       ','fxz2      ','fyz2      ','fz(x2-y2) ', &
  'fxyz      ','fx(3y2-x2)','fy(3x2-y2)'/)

  !
  !     --------------------------------------------------------------
  !
  call mysecond(t0)  
  write(7, 110)  
  nkpt = elec_st%nkpt
  if (.not. (pbc%is_on)) then
     nkpt = 1
  else if (nkpt < 1) then  
     write(7, *) 'DOS: kpoint map is not set up!'  
     ierr=-9987
     stop 
  end if

  
  nmesh = pbc%dos_pnum  
  allocate(fulldos(nmesh, elec_st%nspin),stat=ierr)  
  call alccheck('fulldos',nmesh*elec_st%nspin,ierr)
  fulldos = zero
  allocate(dummy(nmesh), stat =ierr)  
  call alccheck('dummy',nmesh,ierr)
  dummy = zero


  !
  !
  !     rewrite the eigenvalues in eV relative to fermi level
  !
  allocate(et(elec_st%nrep,elec_st%nstate, nkpt, elec_st%nspin), stat=ierr)  
  call alccheck('et',elec_st%nrep*elec_st%nstate*nkpt*elec_st%nspin,ierr)  

  !
  ! compute the Ylm eigenvalues -> ylmwgt=Ylm*(eigenfunctions)
  ! nylm = number of Ylm function used
  !
  if (pbc%ylmdos_l >= 0 .and. pbc%is_on) then
     nylm = 0 
     do  n= 0, pbc%ylmdos_l
        nylm=nylm+n*2+1
     end do        
     allocate(ylmdos(nmesh, nylm, clust%atom_num, elec_st%nspin), stat=ierr)  
     call alccheck('ylmdos',nmesh*nylm*clust%atom_num*elec_st%nspin,ierr)
     allocate(ylmwgt(nylm, clust%atom_num,elec_st%nrep,elec_st%nstate, &
      nkpt, elec_st%nspin), stat=ierr)  
     call alccheck('ylmwgt',nylm*clust%atom_num*elec_st%nrep*elec_st%nstate &
     *nkpt*elec_st%nspin,ierr) 
     call ylm_eig(clust,elec_st,grid,pbc,parallel,ylmwgt,nkpt,nylm,ierr)
     ylmdos = zero
  else 
     nylm=0
  end if
  
  do isp = 1, elec_st%nspin       
     do ikp = 1, nkpt
      do irp = 1, elec_st%nrep  
            do ibn = 1, elec_st%nstate             
              et(irp, ibn, ikp, isp) = (elec_st%eig(irp,ikp,isp)%en(ibn) &
              - elec_st%efermi) * rydberg
            end do
        end do
     end do
  end do
  !
  !     determine lowest and highest eigenvalues
  !
  mineigen = et(1, 1, 1, 1)  
  maxeigen = et(1, 1, 1, 1)
  do isp = 1, elec_st%nspin  
     do ikp = 1, nkpt  
      do irp = 1, elec_st%nrep  
          if (et(irp,1, ikp, isp) < mineigen) mineigen = et(irp, 1, ikp, isp)  
          if (et(irp,elec_st%eig(irp,ikp,isp)%nec, ikp, isp) > maxeigen) maxeigen = &
              et(irp,elec_st%eig(irp,ikp,isp)%nec, ikp, isp) !??
        end do
     end do
  end do

  delta = abs(maxeigen - mineigen) / real(nmesh, dp)
  ebound(1) = mineigen - 1.0d1 * delta  
  ebound(2) = maxeigen + 1.0d1 * delta  
  delta = abs(ebound(2) - ebound(1)) / real(nmesh, dp)    
  if (.not. (pbc%is_on)) then
     nylm = 0     
     do isp = 1, elec_st%nspin  
       do irp = 1, elec_st%nrep          
          do ibn = 1,elec_st%nstate
            n = int((et(irp,ibn, 1, isp)-ebound(1)) / delta +1)            
            fulldos(n,isp) = fulldos(n,isp) + elec_st%eig(irp,1,isp)%occ(ibn) 
          end do        
       end do
     end do 
  else
     ! number of tetrahedra
     ntetra = 6 * elec_st%mpgrid(1) * elec_st%mpgrid(2) * elec_st%mpgrid(3)
     vol = twopi**3 / pbc%vcell / ntetra                           
     allocate(itetdef(4, ntetra), stat =ierr)  
     call alccheck('itetdef',4*ntetra,ierr)
!     allocate(kmapf(nx*ny*nz), stat=ierr)  
!     call alccheck('kmapf',nx*ny*nz,ierr)
     call setup_tetrahedra(elec_st%mpgrid, elec_st%kptgrid, itetdef, pbc%bvec)  !?? 
     
     neqtet = 0  
     idum = 0  
     do itet = 1, ntetra                  ! loop through all tetrahedra
       kidx = itetdef(:, itet)  
       do isp = 1, elec_st%nspin  
          do irp = 1, elec_st%nrep  
           do ibn = 1, elec_st%ntotal(irp)
                ecorn(:) = et(irp,ibn, kidx(:), isp)               ! energy at corners
                fcorn(:) = one
                call spec_tet(vol, ecorn, fcorn, nmesh, ebound, fulldos(1, isp), &
                dummy(1), neqtet)
                if (nylm > 0) then  
                  do iat = 1,clust%atom_num
                       do iylm = 1, nylm  
                         ecorn(:) = et(irp, ibn, kidx(:), isp)   ! energy at corners
                         fcorn(:) = ylmwgt(iylm, iat, irp, ibn, kidx(:), isp)
                         call spec_tet(vol, ecorn, fcorn, nmesh, ebound, &
                              dummy(1), ylmdos(1, iylm, iat, isp), idum)
                      end do
                   end do
                end if
          end do 
          end do
       end do
     end do
     ibn = 0
     do irp = 1, elec_st%nrep
       ibn = ibn + elec_st%ntotal(irp)
     end do   
     neqtet = neqtet / (elec_st%nspin * ibn)  !??
     write(7, 210) ntetra  
     write(7, 200) neqtet
  end if  

  !
  !     ------------ print out ---------------------------------------
  !
  
  open(11, file = 'DOS', status = 'unknown', form = 'formatted')  
  write(11, *) 'Density Of States - Units:  states/eV/unit cell vs. eV'
  write(11, '(a15,2f12.6)') 'Fermi levels: ', elec_st%efermi , &
       elec_st%efermi

  if (elec_st%nspin==1) then 
   write(11, '(a17)') 'Tot. DOS: '  
  else
   write(11, '(a17)') 'Spin DOS: '  
  end if
  if (pbc%is_on) then  
     dos_factor = pbc%vcell / twopi**3 / elec_st%nspin
  else
     dos_factor = one
  end if
  do isp = 1, elec_st%nspin  
     do ie = 1, nmesh  
        write(11, * ) ebound(1) + ie * delta, fulldos(ie, isp) * dos_factor             
     end do
     write(11, * ) '&'  
  end do
  if (elec_st%nspin > 1) then  
     write(11, '(a100)') 'Sum of spin DOS: '  
     do ie = 1, nmesh  
        write(11, * ) ebound(1) + ie * delta, fulldos(ie, 1) * dos_factor + &
             fulldos(ie, 2) * dos_factor
     end do
     write(11, * ) '&'  
  end if
  !
  !     ------------------------------------------------------------------
  !
  !
  !     print Ylm resolved DOS if it is requested
  !  
  if (nylm > 0) then  
     do isp = 1, elec_st%nspin  
        do iat = 1, clust%atom_num  
           do iylm = 1, nylm
              write(11, 124) trim(clust%name(clust%atype(iat))) , iat, iylm, isp 
              write(11, 126) trim(clust%name(clust%atype(iat))) , iat, trim(ylmname(iylm)), isp
              do ie = 1, nmesh  
                 write(11, *) ebound(1) + ie * delta, &
                      ylmdos(ie, iylm, iat, isp) * dos_factor
              end do
              write(11, *) '&'  
           end do
        end do
     end do
  end if
 
  close(11)  

  write(7, 987)  
  call mysecond(t1)  
  write (7, 988) t1 - t0  

  deallocate(itetdef)  
  
  deallocate(fulldos)  
  deallocate(et)
  deallocate(dummy)
  if (pbc%ylmdos_l >= 0 .and. pbc%is_on) then
     deallocate(ylmwgt)
     deallocate(ylmdos)
  end if
    
  call myflush(7)
  
  return  

110 format(/' Generation of DOS : ',/1x,17('-'))  
124 format('LDOS : Atom Type = ',a2,', Atom Num =', i2, ', Ylm Num = ',i2,', ', 'Spin = ',i2)
126 format('Graph Name : ',a, i1,'_',a,'_sp',i1)
200 format(' Number of tetrahedra will eqaul corners:',i10)  

210 format(' Tot. number of tetrahedra :             ',i10)  

987 format(/' The DOS data has been writted to the file "DOS"'/)  

988 format(/' TIME FOR DENSITY OF STATES: ',f14.4)  

end subroutine dos
!
!     =================================================================
!
subroutine setup_tetrahedra(grid, k_index, itetra, bvec)  
  !
  !     1996 Bernd Pfrommer
  !     originated from a code by Paul Delaney
  !
  use constants
  implicit none  
  !
  !     INPUT:
  !     -----
  !  
  integer, intent(in) :: grid(3), &  ! the Monkhorst pack k-grid size kx,ky,kz
       k_index(0:grid(1)-1, 0:grid(2)-1, 0:grid(3)-1)         ! map to irred k
  real(dp), intent(in) :: bvec(3, 3)  ! the reciprocal lattice vectors, as col
                                                               ! in the matrix
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: &
       itetra(4, 6 * grid(1) * grid(2) * grid(3))     ! index into irred. k-pt
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     1) Finds optimal 6-tetrahedral disection of the monkhorst-pack
  !     grid, such that the edgelengths of the tetrahedra are as
  !     regular as possible
  !
  !     2) With this, generates the index array itetra which maps each
  !     corner of the 6*nx*ny*nz tetrahedra to an irreducible kpoint
  !
  !
  !
  !     -------------------------- local variables ------------------
  !
  integer :: ntetra, i, j, k, it, iv, idx, idy, idz, ivert, tset(3, 4, 6), &
       ivert2, ivert1, goodi, goodj, goodk
  real(dp) :: thisedgemin, thisedgemax, cart(3, 4), edge, edgemin, edgemax
  ! tolerance in the finding the right k point
  real(dp), parameter :: tol_bz = 1.d-8
 
   !
  !
  !

!!$  integer, parameter :: initialtetraoffset(3, 4, 6) = (/ &
!!$       0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, &
!!$       0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, &
!!$       1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, &
!!$       1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1 /)

  ! This way in F90 - CJP 16/04/00

  integer, parameter :: initialtetraoffset(3,4,6) = reshape(source=(/ &
       0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, &
       0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, &
       1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, &
       1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1 /),shape=(/3,4,6/))

  edgemin = zero
  edgemax = 1.0d10  

  !     First set up the matrix of components of the reciprocal lattice
  !     vectors along the Cartesian basis. bcomp(i,j) holds the dot
  !     product of dual basis vector b_i with Cartesian basis vector
  !     e_j.
  !     Now go through all the four long diagonals, specified by the
  !     corner (i,j,k)

  k = 0  
  do i = 0, 1  
     do j = 0, 1  
        do ntetra = 1, 6  
           do ivert = 1, 4  
              call mirrorijk(initialtetraoffset(1, ivert, ntetra), &
                   tset (1, ivert, ntetra), i, j, k)
           end do
        end do
        thisedgemin = 1.0d10  
        thisedgemax = zero
        do ntetra = 1, 6  
           do ivert = 1, 4  
              cart(:, ivert) = matmul(bvec, real(tset(:, ivert, ntetra), dp))  
           end do
           do ivert1 = 1, 3  
              do ivert2 = ivert1 + 1, 4  
                 edge = (cart(1, ivert1) - cart(1, ivert2))**2 + &
                      (cart(2, ivert1) - cart(2, ivert2))**2 + &
                      (cart(3, ivert1) - cart(3, ivert2))**2
                 edge = sqrt(edge)  
                 thisedgemax = max(thisedgemax, edge)  
                 thisedgemin = min(thisedgemin, edge)  
              end do
           end do
        end do
        if (thisedgemax < edgemax) then  
           goodi = i  
           goodj = j  
           goodk = k  
           edgemax = thisedgemax  
           edgemin = thisedgemin  
        end if
     end do
  end do
  do ntetra = 1, 6  
     do ivert = 1, 4  
        call mirrorijk(initialtetraoffset(1, ivert, ntetra), &
             tset (1, ivert, ntetra), goodi, goodj, goodk)
     end do
  end do
  !      write(0,*)
  !      write(0,*) 'Choosing long-diagonal based at i,j,k=',
  !     $     goodi,goodj,goodk
  !      write(0,*) 'The longest edge length isp ',edgemax,' a.u.'
  !      write(0,*) 'The shortest isp ',edgemin,' a.u.'
  !     ------------------------------------------------------------------
  !     ----------- now use the optimial tetrahedron to set up the -------
  !     -------------------------- index array ---------------------------
  !
  !      write(7,*) 'TETRAHEDRON GENERATION:'

  ntetra = 0  
  do i = 0, grid(1) - 1  
     do j = 0, grid(2) - 1  
        do k = 0, grid(3) - 1  
           do it = 1, 6  
              !                  write(7,*) 'tetrahedron:',it
              do iv = 1, 4  
                 idx = mod(i + tset(1, iv, it), grid(1))  
                 idy = mod(j + tset(2, iv, it), grid(2))  
                 idz = mod(k + tset(3, iv, it), grid(3))  
                 itetra(iv, ntetra + it) = abs(k_index(idx, idy, idz))  
                 !               write(7,*) idx,idy,idz,itetra(iv,ntetra+it)
              end do
           end do
           ntetra = ntetra + 6  
        end do
     end do
  end do

end subroutine setup_tetrahedra
!
!     =================================================================
!
subroutine mirrorijk(in, out, i, j, k)  
  !
  !     Do mirrors in i around i=.5 if i=1; similarly for j,k
  !
  !     1995/96 by Paul Delaney
  !
  integer, intent(in) :: i, j, k  
  integer, intent(in) :: in(3)
  integer, intent(out) :: out (3)  
  out (1) = in (1)  
  out (2) = in (2)  
  out (3) = in (3)  
  if (i == 1) out(1) = 1 - out(1)  
  if (j == 1) out(2) = 1 - out(2)  
  if (k == 1) out(3) = 1 - out(3)
  
  return  

end subroutine mirrorijk

subroutine ylm_eig(clust,elec_st,grid,pbc,parallel,ylmwgt,nkpt,nylm,ierr)
  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pbc_module
  use parallel_data_module
  implicit none  


  !     INPUT:
  !     -----
  !
  !     the cluster
      type (cluster), intent(in) :: clust
  !     electronic structure
      type (electronic_struct), intent(in) :: elec_st
  !     grid related data
      type (grid_data), intent(in) :: grid
  !     pbc related data
      type (pbc_data), intent(in) :: pbc
  !     parallel computation related data
      type (parallel_data), intent(in) :: parallel      
  !     number of kpoint used
      integer, intent(in) :: nkpt
  !     number of Ylm function used
      integer, intent(in) :: nylm
  !     error flag
      integer, intent(out) :: ierr 

  real(dp), intent(inout) :: ylmwgt(nylm, clust%atom_num,elec_st%nrep, &
      elec_st%nstate,nkpt, elec_st%nspin)

  !   DESCRIPTION:
  !   -----------
  !
  !   Computes the dot product of the Ylm functions of each atom and the
  !   different eigenfunctions
  !
  !   -------------------local variables ------------------------------  
  !
  
  real(dp) :: ylm(16), uvec(3), xvec(3), rr(3), xx, &
     yy, zz, rtmp, rsize, rnorm, rcomp, sigma, fac1, fac2, ylmsum
  real(dp), allocatable :: rgl(:), & ! radial positions in Gauss-Leg. integral
      wgl(:),    &         ! radial weights in Gauss-Legendre integral
      ylmwgl(:,:,:,:,:,:)  ! temporal array to store the ylm weight for
                           ! every radial point  
  complex(dp), allocatable :: zylmwgl(:,:,:,:,:,:) ! same, but complex
  complex(dp) :: ztmp
  integer :: igr,irp,isp,ikp,iat,iylm, ibn, icellx, icelly, &
             icellz, ngl, maxngl, igl
  
  
  
  !  Spherical harmonics coefficients - for s, p, and d
  !  formulae, see, e.g., pp. 11-12 in W. A. Harrison, "Electronic
  !  Structure and the Properties of Solids: The Physics of the
  !  Chemical Bond". For f formulae, see, e.g., 
  !  S. F. A. Kettle, "Physical Inorganic Chemistry", p. 244.
  !  Notice that they *have* the factor sqrt( (2*l+1)/4*pi ).
  !
  !  1/(2*sqrt(pi)) for s
  real(dp), parameter :: c0  = 0.282094791773878143474039725779d0
  !  sqrt(3/(4*pi)), for px, py, pz 
  real(dp), parameter :: c1  = 0.488602511902919921586384622837d0
  !  sqrt(15/(4*pi)), dxy, dzy, dzx 
  real(dp), parameter :: c21 = 1.092548430592079070543385705802d0
  !  sqrt(5/(16*pi)), for dz^2 
  real(dp), parameter :: c22 = 0.315391565252520006030893690294d0
  !  sqrt(7/(16*pi)), for f_z(5z^2-3r^2)
  real(dp), parameter :: c31 = 0.373176332590115391414395913197d0
  !  sqrt(21/(32*pi)), for f_x(5z^2-r^2),f_y(5z^2-r^2)
  real(dp), parameter :: c32 = 0.457045799464465736158020696916d0
  !  sqrt(105/(16*pi)), for f_z(x^2-y^2)
  real(dp), parameter :: c33 = 1.445305721320277027694690077198d0
  !  sqrt(105/(4*pi)), for f_xyz
  real(dp), parameter :: c34 = 2.890611442640554055389380154398d0
  !  sqrt(35/(32*pi)), for f_x(x^2-3y^2),f_y(y^2-3x^2)
  real(dp), parameter :: c35 = 0.590043589926643510345610277541d0  
   
  ylmwgt=zero 
  ierr=0
  
  maxngl=0
  sigma = minval(grid%step)


  ! notice that the ngl here must be the same as set within the loop
  do iat = 1,clust%atom_num
     maxngl=max(maxngl,max(10,nint(clust%ylmdos_cutoff(clust%atype(iat))*5.0/sigma)))
  end do
  allocate(rgl(maxngl))
  allocate(wgl(maxngl))
   
  if (elec_st%cplx) then  
     allocate(zylmwgl(nylm, maxngl,elec_st%nrep,elec_st%nstate, &
     nkpt, elec_st%nspin), stat=ierr) 
     call alccheck('zylmwgl',nylm*maxngl*elec_st%nrep*elec_st%nstate &
     *nkpt*elec_st%nspin,ierr)  
  else
     allocate(ylmwgl(nylm, maxngl,elec_st%nrep,elec_st%nstate, &
     nkpt, elec_st%nspin), stat=ierr) 
     call alccheck('ylmwgl',nylm*maxngl*elec_st%nrep*elec_st%nstate &
     *nkpt*elec_st%nspin,ierr)  
  end if

  do iat = 1,clust%atom_num
  do isp = 1, elec_st%nspin  
  do ikp = 1, nkpt  
  do irp = 1, elec_st%nrep
  do ibn = 1, elec_st%nstate   
  do iylm=1,nylm
      ylmwgt(iylm,iat,irp,ibn,ikp,isp) = zero
  end do
  end do
  end do
  end do
  end do
  end do


              
  call determinant(pbc%avec_norm,fac1)  
  fac1 = sqrt(grid%step(1) * grid%step(2) * grid%step(3) * fac1) !integration factor, volume of integration step    
  fac2 = sqrt(two * pi) * sigma  ! Gaussian normalization
  do iat = 1,clust%atom_num
     ! setting the weights for the Gauss-Legendre radial grid
     ngl = max(10,nint(clust%ylmdos_cutoff(clust%atype(iat))*5.0/sigma))
     call mygauleg(0.d0,clust%ylmdos_cutoff(clust%atype(iat)),rgl,wgl,ngl)
     if (elec_st%cplx) then
          zylmwgl=zero
     else 
          ylmwgl=zero
     end if     
     do igl = 1 , ngl
         do igr = 1 , grid%ndim
         ! shifting the box to make sure that no point that is in the cutoff
         ! sphere is missed
         do icellx = -1, 1
         do icelly = -1, 1
         do icellz = -1, 1         
           uvec(1)=(grid%shift(1)+grid%fx(igr))*grid%step(1)+icellx*pbc%box_size(1)
           uvec(2)=(grid%shift(2)+grid%fy(igr))*grid%step(2)+icelly*pbc%box_size(2)
           uvec(3)=(grid%shift(3)+grid%fz(igr))*grid%step(3)+icellz*pbc%box_size(3)
           call matvec3('N',pbc%avec_norm,uvec,xvec)        
           !calculating the angular coefficients
           rr(1)=xvec(1)-clust%xatm(iat)
           rr(2)=xvec(2)-clust%yatm(iat)
           rr(3)=xvec(3)-clust%zatm(iat)
           rsize=sqrt(dot_product(rr,rr))
           !if the distance is larger the defined cutoff skip the point
           if (rsize > clust%ylmdos_cutoff(clust%atype(iat))) cycle           
           rcomp = exp(-((rsize - rgl(igl))**2/two/(sigma**2))) / (rsize**2) /fac2
           !computing a radial compunent, made of (r-rgl)^2 gaussian and 1/r^2 components
           !the exponent is equivelant to a delta function.
           rnorm=one/rsize
           xx=rr(1)*rnorm
           yy=rr(2)*rnorm
           zz=rr(3)*rnorm
           ylm(1)=c0 !compting the 
           if (nylm>=4) then
              ylm(2) = rcomp*c1*xx
              ylm(3) = rcomp*c1*yy
              ylm(4) = rcomp*c1*zz
           end if
           if (nylm>=9) then
              ylm(5) = rcomp*c21*xx*yy
              ylm(6) = rcomp*c21*zz*yy
              ylm(7) = rcomp*c21*xx*zz
              ylm(8) = rcomp*c22*(three*zz*zz-one)
              ylm(9) = rcomp*c21*half*(xx*xx-yy*yy)
           end if
           if (nylm==16) then
              ylm(10) = rcomp*c31*zz*(five*zz*zz-three)
              ylm(11) = rcomp*c32*xx*(five*zz*zz-one)
              ylm(12) = rcomp*c32*yy*(five*zz*zz-one)
              ylm(13) = rcomp*c33*zz*(xx*xx-yy*yy)
              ylm(14) = rcomp*c34*xx*yy*zz
              ylm(15) = rcomp*c35*xx*(xx*xx-three*yy*yy)
              ylm(16) = rcomp*c35*yy*(yy*yy-three*xx*xx)                   
           end if
           !calculating the (eigenfunction)*(ylm), summing
           !over all grid points within the cutoff radius
           if (elec_st%cplx .and. parallel%mxwd /= 2) then 
              do isp = 1, elec_st%nspin  
              do ikp = 1, nkpt  
              do irp = 1, elec_st%nrep
              do ibn = 1, elec_st%eig(irp,ikp,isp)%nec   
              do iylm=1,nylm
                  zylmwgl(iylm,igl,irp,ibn,ikp,isp)=zylmwgl(iylm,igl,irp,ibn,ikp,isp)+&
                  elec_st%eig(irp,ikp,isp)%zwf(grid%rindex(igr),ibn) * ylm(iylm) * fac1
              end do
              end do
              end do
              end do
              end do
           else if (elec_st%cplx .and. parallel%mxwd == 2) then
              do ikp = 1, nkpt  
              do irp = 1, elec_st%nrep
              do ibn = 1, elec_st%eig(irp,ikp,isp)%nec   
              do iylm=1,nylm
                  zylmwgl(iylm,igl,irp,ibn,ikp,1)=zylmwgl(iylm,igl,irp,ibn,ikp,1)+&
                  elec_st%eig(irp,ikp,1)%zwf(grid%rindex(igr),ibn) * ylm(iylm) & 
                  * fac1 * (half-elec_st%magmom(ibn,ikp))
                  zylmwgl(iylm,igl,irp,ibn,ikp,2)=zylmwgl(iylm,igl,irp,ibn,ikp,2)+&
                  elec_st%eig(irp,ikp,1)%zwf(grid%rindex(igr+parallel%mydim),ibn) * ylm(iylm) &
                  * fac1 * (half+elec_st%magmom(ibn,ikp))
              end do
              end do
              end do
              end do
           else          
              do isp = 1, elec_st%nspin  
              do ikp = 1, nkpt  
              do irp = 1, elec_st%nrep
              do ibn = 1, elec_st%eig(irp,ikp,isp)%nec   
              do iylm=1,nylm
                  ylmwgl(iylm,igl,irp,ibn,ikp,isp)=ylmwgl(iylm,igl,irp,ibn,ikp,isp)+&
                  elec_st%eig(irp,ikp,isp)%wf(grid%rindex(igr),ibn) * ylm(iylm) * fac1
              end do
              end do
              end do
              end do
              end do
           end if
        end do
        end do
        end do
        end do !end of 4 real space grid loops  
     end do    ! end of radial grid loop  
       
     if (elec_st%cplx) then 
        do igl = 1, ngl
        do isp = 1, elec_st%nspin  
        do ikp = 1, nkpt  
        do irp = 1, elec_st%nrep
        do ibn = 1, elec_st%eig(irp,ikp,isp)%nec  
        do iylm=1,nylm
            ylmwgt(iylm,iat,irp,ibn,ikp,isp)=ylmwgt(iylm,iat,irp,ibn,ikp,isp)+&
             abs((rgl(igl)*zylmwgl(iylm,igl,irp,ibn,ikp,isp))**2)*wgl(igl)                          
        end do
        end do
        end do
        end do
        end do
        end do
     else 
        do igl = 1, ngl
        do isp = 1, elec_st%nspin  
        do ikp = 1, nkpt  
        do irp = 1, elec_st%nrep
        do ibn = 1, elec_st%eig(irp,ikp,isp)%nec   
        do iylm=1,nylm
            ylmwgt(iylm,iat,irp,ibn,ikp,isp)=ylmwgt(iylm,iat,irp,ibn,ikp,isp)+&
            (rgl(igl)*ylmwgl(iylm,igl,irp,ibn,ikp,isp))**2*wgl(igl)                          
        end do
        end do
        end do
        end do
        end do
        end do
     end if
 !    do igl = 1, ngl
 !    do isp = 1, elec_st%nspin  
 !    do ikp = 1, nkpt  
 !    do irp = 1, elec_st%nrep
 !    do ibn = 1, elec_st%nstate 
 !    ylmsum=sum(ylmwgt(1:nylm,iat,irp,ibn,ikp,isp))
 !    ylmsum=1/ylmsum
 !    do iylm=1,nylm
 !        ylmwgt(iylm,iat,irp,ibn,ikp,isp)=ylmwgt(iylm,iat,irp,ibn,ikp,isp)*ylmsum                                   
 !    end do
 !    end do
 !    end do
 !    end do
 !    end do
 !    end do     

  end do      ! end of atoms loop    

  if (elec_st%cplx) then 
    deallocate(zylmwgl) 
  else 
    deallocate(ylmwgl)
  end if 

  deallocate(rgl)
  deallocate(wgl)


  end subroutine ylm_eig


!     ===============================================================
!
!     Given the lower and upper limits of integration x1 and x2, and
!     given n, this routine returns arrays x(1:n) and w(1:n),
!     containing the abscissas and weights of the Gauss-Legendre
!     n-point quadrature formula.
!     from Numerical Recipes in Fortran
!
!     ---------------------------------------------------------------
subroutine mygauleg(x1, x2, x, w, n)  
    use constants
    implicit none  
!
!     Input/Output variables:
!
    real(dp), intent(in) :: x1, x2  
    integer, intent(in) :: n  
    real(dp), intent(out) :: w(n), x(n)  
!
!     Work variables:
!
    real(dp) :: p1, p2, p3, pp, xl, xm, z, z1, dj
    integer :: i, j, m  
    real(dp), parameter :: eps = 3.0d-14  

!     ---------------------------------------------------------------

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
!
!     =================================================================
!      
!     Computes a 3x3 matrix determinant

subroutine determinant(m,det)
    use constants
    implicit none  
    real(dp), intent(in) :: m(3,3)
    real(dp), intent(out) :: det
    
    real(dp) :: a(3,3)
    integer i,j
    !     ---------------------------------------------------------------
    !
    !     compute matrix of cofactors
    !
    a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
    a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
    a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
    a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
    a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
    a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
    a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
    a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
    a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
    !
    !     compute determinant
    !
    det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
end subroutine determinant
!     ===============================================================      
