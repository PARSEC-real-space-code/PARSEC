subroutine calc_bands(band_st,elec_st,pbc,pot,u_pot,solver,parallel, &
           grid, nloc_p_pot, clust, p_pot, ipr,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use potential_module
  use non_local_psp_module
  use eigen_solver_module
  use parallel_data_module
  use pseudo_potential_module
  use bandstruc_module
  use pbc_module
  use grid_module
  
  implicit none

  !
  ! Input/Output variables:
  !
  ! bands structure
  type (bandstruc), intent(inout) :: band_st
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !grid
  type (grid_data), intent(in) :: grid
  ! pbc
  type (pbc_data), intent(inout) :: pbc
  ! potential related data
  type (potential), intent(in) :: pot
  ! on-site Coulomb interaction related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! non local pseudo potential part
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot
  ! pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! cluster data
  type (cluster), intent(in) :: clust
  ! printout flag
  integer, intent(in) :: ipr
  ! error flag, 400 < ierr < 421
  integer, intent(out) :: ierr

  !
  ! Work variables:
  !
  integer :: istate,isp,kplp,scf_nkpt,i,offset,j,lastkpt,spnum,jj,nrep,nmax
  real(dp), dimension(3)  :: vec, ktmp
  real(dp) :: min_length  
  !


!until we figure out what to do with symmetry
nrep = 1

!conversion of kpoints of start and end of lines to
!cartesian coordinates of inverse Bohr radius
!by multiplying by reciprocal lattice verctors

  do i = 1, band_st%nlines
     call matvec3('N',pbc%bvec,band_st%blines(i)%start,vec)
     band_st%blines(i)%startb= vec
     call matvec3('N',pbc%bvec,band_st%blines(i)%end,vec)
     band_st%blines(i)%endb= vec
  enddo
 
!find the length of the shortest line 

  do i = 1, band_st%nlines
     vec=band_st%blines(i)%endb-band_st%blines(i)%startb
     band_st%blines(i)%length= sqrt(dot_product(vec,vec))
  enddo
  
  min_length=band_st%blines(1)%length
  if (band_st%nlines >1) then
     do i=2,band_st%nlines
        if (band_st%blines(i)%length < min_length) min_length=band_st%blines(i)%length
     enddo
  endif

!calculate the number of kpoints for each line
!proportional to the length of the line
!for the last line the end point is also included

 do i = 1, band_st%nlines-1
    band_st%blines(i)%nkpt = band_st%blines(i)%length/min_length*band_st%npoints
 enddo 
 band_st%blines(band_st%nlines)%nkpt = band_st%blines(band_st%nlines)%length/min_length*band_st%npoints +1

!create list of kpoints for each line

  do i = 1, band_st%nlines
     allocate(band_st%blines(i)%kpts(3,band_st%blines(i)%nkpt))
     band_st%blines(i)%kpts = zero
     if (i< band_st%nlines) then 
        do j = 1,band_st%blines(i)%nkpt
        band_st%blines(i)%kpts(:,j)= band_st%blines(i)%startb + &
        (j-1)*(band_st%blines(i)%endb-band_st%blines(i)%startb)/band_st%blines(i)%nkpt
        enddo     
     else
        do j = 1,band_st%blines(i)%nkpt
        band_st%blines(i)%kpts(:,j)= band_st%blines(i)%startb + &
        (j-1)*(band_st%blines(i)%endb-band_st%blines(i)%startb)/(band_st%blines(i)%nkpt-1)
        enddo
     endif
  enddo

!devide each line to batches of kpoints and find eigenvalues for each batch      

  spnum = elec_st%nspin/elec_st%mxwd  
  scf_nkpt = elec_st%nkpt

! some array initializations...

 if (associated(solver%eig_init)) deallocate(solver%eig_init)
 allocate(solver%eig_init(nrep,scf_nkpt,spnum))
 solver%eig_init(:,:,:) = .false.
 if (associated(solver%nconvt)) deallocate(solver%nconvt)
 allocate(solver%nconvt(nrep,scf_nkpt,spnum))
 solver%nconvt(:,:,:) = 0
 if (associated(solver%eval_loc)) deallocate(solver%eval_loc)
 nmax = elec_st%nstate + solver%nadd + solver%winsize
 allocate(solver%eval_loc(nmax,nrep,spnum))
 if (solver%name == ARPACK) then
    if (associated(solver%zres_res)) deallocate(solver%zres_res)
    allocate(solver%zres_res(parallel%ldn*parallel%mxwd,nrep,scf_nkpt,spnum))
    solver%zres_res(:,:,:,:) = zzero
    if (associated(solver%info)) deallocate(solver%info)
    allocate(solver%info(nrep,scf_nkpt,spnum))
    solver%info = 0
 endif

!main loop- i is loop over lines, j is loop over batches of k-points in each line
  do i = 1, band_st%nlines

     elec_st%kpts = zero
     solver%kecoe1(:,:,:) = zero

     allocate(band_st%blines(i)%eigs(band_st%blines(i)%nkpt,spnum,elec_st%nstate))
     band_st%blines(i)%eigs=zero
     offset=0
     pbc%nkpt = elec_st%nkpt
     nloc_p_pot%nkpt = elec_st%nkpt

     do while (offset < band_st%blines(i)%nkpt) 
        
        lastkpt = min(offset + scf_nkpt, band_st%blines(i)%nkpt)
       
        elec_st%kpts = zero
        
        do j= offset+1, lastkpt
           elec_st%kpts(:,j-offset)= band_st%blines(i)%kpts(:,j)
           pbc%kpts(:,j-offset)=elec_st%kpts(:,j-offset)
           nloc_p_pot%kpts(:,j-offset)= elec_st%kpts(:,j-offset)
        enddo

        do kplp = 1, elec_st%nkpt
           do jj =0,grid%norder
           call matvec3('T',grid%grad_bvec_norm,elec_st%kpts(:,kplp),ktmp)
           solver%kecoe1(kplp,1:3,jj) = -2*zi*grid%coe1(jj,1:3)*ktmp
           enddo
        enddo
    
        solver%eig_init(:,:,:) = .false.

        call nonloc(clust,grid,p_pot,nloc_p_pot,pbc,parallel,ierr)

        call eigval(elec_st,pot,u_pot,solver,parallel,ipr,ierr)

        do j= offset+1, lastkpt
           do isp=1, spnum           
            band_st%blines(i)%eigs(j,isp,1:elec_st%nstate) & 
               = elec_st%eig(1,j-offset,isp)%en
           enddo
        enddo
      
       offset = lastkpt;
            
     enddo
  enddo

!write results to file
!eigenvalues are shifted so that Ef=0

if (parallel%iammaster) then

   open(88, file='bands.dat',form='formatted',status='unknown')


        ! Number of spins, number of states for each k-points, 
        ! number of band segments, Fermi energy (in Ry)
        write(88,*) '# of spins, # of states, # of segments, Ef (Ry)'
        write(88,*) spnum, elec_st%nstate, band_st%nlines, elec_st%efermi

        ! Number of k-points of each segment
        write(88,*) 'idx of segment, # of k-points'
        do jj = 1, band_st%nlines
                write(88,*) jj, band_st%blines(jj)%nkpt
        enddo

        ! For each spin, print out one line of (spin, # seg., # kp, kp coord)
        ! Then the eigenvalues (in Ry)
        write(88,*) 'spin, segment, kpt, kpt coord. (1/bohr), &
                &followed by eigs (Ry)'

        do isp = 1, spnum    
        do i = 1, band_st%nlines
        do j = 1, band_st%blines(i)%nkpt
                write(88,'(3i6,3(3x,f12.8))') isp, i, j, &
                        band_st%blines(i)%kpts(:,j)

                do istate = 1, elec_st%nstate
                        write(88,*) &
                        band_st%blines(i)%eigs(j,isp,istate)
                enddo
        enddo
        enddo
        enddo




   !do i=1, band_st%nlines
   !   do j=1,band_st%blines(i)%nkpt
   !        write(88,'(3i6,3(3x,f12.8))') isp,i,j, band_st%blines(i)%kpts(:,j)
   !        do istate = 1, elec_st%nstate
   !            do isp = 1,spnum
   !             write (88,'(i6,2(3x,f15.6))') &
   !              istate,band_st%blines(i)%eigs(j,isp,istate),&
   !              band_st%blines(i)%eigs(j,isp,istate)*rydberg
   !            enddo
   !        enddo
   !   enddo
   !enddo

   close(88)

   open(88, file="bands_plot.dat",form='formatted',status='unknown')
   do istate = 1, elec_st%nstate
      jj = 0
      do i=1, band_st%nlines
          do j=1,band_st%blines(i)%nkpt
              jj = jj + 1
              write(88,'(i6,3x,2(3x,f12.8))') jj,&
               (band_st%blines(i)%eigs(j,isp,istate)*rydberg,isp=1,spnum)
          enddo
      enddo
      write(88,'(1x)')
   enddo
   close(88)

endif

end subroutine calc_bands 

