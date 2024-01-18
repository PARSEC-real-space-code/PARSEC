!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine minimizes total energy using either one of 2
! schemes:
! (1) steepest descent, moving the atoms in the direction of the
! force acting on them; or 
! (2) BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm,
! implemented with the Limited Memory BFGS routines by 
! J. Nocedal.
! If instead the user desires to manually give a sequence of
! coordinates to be used, those coordinates are read from 
! manual.dat. This is useful for checking energy and force minima
! along pre-defined paths of atom motion.
!
!---------------------------------------------------------------
subroutine domove(clust,move,pbc,etot,bindev,ipr,ierr)

  use constants
  use cluster_module
  use movement_module
  use pbc_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! movement data
  type (movement), intent(inout) :: move
  ! periodic boundary conditions
  type (pbc_data), intent(in) :: pbc
  ! total energy in Rydberg, total energy per atom in eV
  real(dp), intent(in) :: etot,bindev
  ! printout flag (for BFGS)
  integer, intent(in) :: ipr
  ! error flag (for BFGS), 540 < ierr < 551
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! counters
  integer i,j,ity
  ! largest force component found, magnitude of forces
  real(dp) :: fmax, flen
  ! number of movable atoms
  integer nmovable
  ! negative of atomic forces of the movable atoms
  real(dp), allocatable :: mforce(:,:)
  !
  ! steepest descent variables:
  ! increments of positions
  real(dp) :: disp(3)
  ! scaling coefficient in formula of position increments
  real(dp) :: coeff
  ! max step(au), min step (au), absolute value of minimum force
  real(dp) :: stepmax, stepmin, abs_fmin
  !
  ! BFGS variables:
  integer, dimension(:), allocatable :: nbound
  integer nparam,iprint
  real(dp) :: xtol

  !---------------------------------------------------------------

  nmovable = sum(clust%mvat)
  nparam = 3*nmovable
  xtol = 1.d-15
  if (move%name == BFGS) then
     if (ipr >= 3) then
        iprint = 100
     elseif (ipr >= 4) then
        iprint = 101
     else
        iprint = 1
     endif
     allocate(nbound(nparam))
     if (move%stepmax > zero) then
        nbound = 2
     else
        nbound = 0
     endif
  endif
  !
  ! Select the forces on movable atoms, and compute the maximum
  ! component (note the negative sign in front of clust%force)
  !
  allocate(mforce(3,nmovable))

  j = 0
  fmax = zero
  do i = 1, clust%atom_num
     if (clust%mvat(i) == 1) then
        j = j + 1
        mforce(:,j) = -clust%force(:,i)
        flen = sqrt(dot_product(mforce(:,j),mforce(:,j)))
        if (flen > fmax) fmax = flen
     endif
  enddo
  !
  ! report total energy and max force belonging to the input
  ! coordinates
  !
  write(66,30) bindev, fmax
30 format('Total Energy = ',f12.5,1x,'[eV/at],' &
        ,1x,'Max force =',f9.5,1x,'[Ry/au]')
  write(66,*) 
  !
  ! If maximal force is below the user-supplied required minimal
  ! force then set 'minimum found' flag, report to files, and
  ! return.
  !
  if(fmax < move%fmin) then
     move%min_frc_found = .true.
     write(7,*) ' Local minimum of total energy was found '
     write(66,*) 'Local minimum of total energy was found '
     deallocate(mforce)
     if ( allocated(nbound) ) deallocate(nbound)
     return
  endif
  !
  ! In first iteration, initialize the coordinates of movable atoms
  ! and other parameters.
  !
  if (move%num == 0) then
     j = 0
     allocate(move%rcur(3,nmovable))
     do i = 1, clust%atom_num
        if (clust%mvat(i) == 1) then
           j = j + 1
           move%rcur(1,j) = clust%xatm(i)
           move%rcur(2,j) = clust%yatm(i)
           move%rcur(3,j) = clust%zatm(i)
        endif
     enddo
     ! for BFGS
     if (move%name == BFGS) then
        allocate(move%lower(3*nmovable))
        allocate(move%upper(3*nmovable))
        allocate(move%ibfgs(44+9*nmovable))
        allocate(move%wbfgs((2*move%mmax + 5)*3*nmovable + 11 &
             *move%mmax*move%mmax + 8*move%mmax))
        allocate(move%lbfgs(4))
        allocate(move%dbfgs(29))

        if (move%is_restart) then
           write(7,*) ' RESTART LBFGS FROM PREVIOUS RUN'
           open(45,file='relax_restart.dat',form='formatted',status='old')
           do j=1,clust%atom_num + 4
              read(45,*)
           enddo
           read(45,*) move%mmax
           read(45,'(60a)') move%cbfgs(1:60)
           read(45,'(60a)') move%cbfgs(61:120)
           do j=1,4
              read(45,*) move%lbfgs(j)
           enddo
           do j=1,3*nmovable
              read(45,*) move%lower(j)
              read(45,*) move%upper(j)
           enddo
           do j=1,44+9*nmovable
              read(45,*) move%ibfgs(j)
           enddo
           do j=1,size(move%wbfgs)
              read(45,*) move%wbfgs(j)
           enddo
           do j=1,29
              read(45,*) move%dbfgs(j)
           enddo
           close(45)
        else
           if (move%stepmax > zero) then
              do i =1, nmovable
                 do j = 1,3
                    move%lower(3*(i-1) + j) = move%rcur(j,i) - move%stepmax
                    move%upper(3*(i-1) + j) = move%rcur(j,i) + move%stepmax
                 enddo
              enddo
           else
              move%lower(:) = zero
              move%upper(:) = zero
           endif
           !
           ! Initialize internal BFGS parameters
           !
           move%cbfgs = 'START'
! #ifdef AJB_DEBUG
!      write(7,*)
!      write(7,*) '       now calling setulb I'
!      write(7,*)
! #endif
           call setulb(nparam,move%mmax, &
                move%rcur,move%lower,move%upper,nbound,etot,mforce &
                ,xtol,xtol,move%wbfgs,move%ibfgs(45) &
                ,move%cbfgs(1:60),iprint,move%cbfgs(61:120) &
                ,move%lbfgs,move%ibfgs(1),move%dbfgs)
        endif
     endif
  endif
  !
  ! Move the atoms according to the desired scheme:
  !
  select case (move%name)
     ! perform simple movement of atoms (steepest descent)
  case (STEEPDESC)

     ! The movement is related to the force component by  dx=stepmax
     ! *fx/(const+fx), and the same for dy, dz. For small f, dx->f
     ! /const. The choice of the constant made here guarantees that
     ! for f=fmin the movement is stepmin - the smallest movement
     ! allowed by the user. For large f, dx->stepmax, and hence
     ! stepmax - the largest movement allowed by the user, even if the
     ! force is huge. stepmax controls how fast the motion  proceeds.
     ! If it's too small, many steps will be needed. It it's too large,
     ! the coordinates will oscillate around the true minimum
     ! energy positions.  The best values for stepmin and stepmax
     ! strongly depend on how "hard" or "soft" the system under study
     ! is.
     stepmax = move%stepmax
     stepmin = move%stepmin
     abs_fmin = abs(move%fmin)
     !
     ! steepest descent with "adjustable" movement
     !
     coeff = log((stepmax-stepmin)/(stepmax-50*stepmin))/(1000*move%fmin)
     do i = 1, nmovable
        do j=1,3
           if (mforce(j,i) == zero) then
              disp(j) = zero
           else
              disp(j) =  -(mforce(j,i)/abs(mforce(j,i)))* &
                   (stepmin+(stepmax-stepmin)* &
                   (1-exp(-coeff*(abs(mforce(j,i))-abs_fmin))))
           endif
        enddo
        move%rcur(:,i) = move%rcur(:,i) + disp(:)
     enddo

     ! perform BFGS minimization of energy
  case (BFGS)
     do
! #ifdef AJB_DEBUG
!      write(7,*) '       now calling setulb II'
! #endif
        call setulb(nparam,move%mmax,move%rcur,move%lower, &
             move%upper,nbound,etot,mforce,xtol,xtol,move%wbfgs, &
             move%ibfgs(45),move%cbfgs(1:60),iprint, &
             move%cbfgs(61:120),move%lbfgs,move%ibfgs(1), &
             move%dbfgs)
! #ifdef AJB_DEBUG
!      write(7,*) '       finished setulb II'
! #endif
        ! BFGS routines need forces/energy at current coordinates
        if (move%cbfgs(1:2) == 'FG') then
           exit
        endif
        if (move%cbfgs(1:4) == 'CONV') then
           write(7,*) 'LBFGS minimization converged '
           write(7,*) move%cbfgs(1:60)
           write(66,*) 'Local minimum of total energy was found'
           move%done = .true.
           write(7,*) 'LBFGS is terminating'
           write(66,*) 'LBFGS is terminating'
           return
        endif
        if (move%cbfgs(1:4) == 'WARN') then
           write(7,*) 'LBFGS minimization warning'
           write(7,*) move%cbfgs(1:60)
           move%done = .true.
           write(7,*) 'LBFGS is terminating'
           write(66,*) 'LBFGS is terminating'
           return
        endif
        ! Error exits
        if (move%cbfgs(1:5) == 'ERROR') then
           write(7,*) 'ERROR in LBFGS routines : '
           write(7,*) move%cbfgs(1:60)
           move%done = .true.
           write(7,*) 'LBFGS is terminating'
           write(66,*) 'LBFGS is terminating'
           ierr = 541
           return
        endif
        if (move%cbfgs(1:4) == 'ABNO') then
           write(7,*) 'Abnormal exit from LBFGS routines : '
           write(7,*) move%cbfgs(1:60)
           move%done = .true.
           write(7,*) 'LBFGS is terminating'
           write(66,*) 'LBFGS is terminating'
           ierr = 542
           return
        endif
     enddo
     deallocate(nbound)

     ! no miminization, move atoms according to user-supplied
     ! coordinates
  case (MANUAL)
     ! Unless the results reported were for the final step, 
     ! read coordinates for next position from 'manual.dat'
     if (move%num < move%mxmove) then
        do i=1,nmovable
           read(67,*) move%rcur(:,i)
        enddo
        ! unless reading the last set of coordinates, skip empty line 
        ! between coordinate sets
        if (move%num < move%mxmove-1) read(67,*)
     endif
  end select
  !
  ! Update atom coordinates in clust structure (notice that only
  ! the coordinates of movable atoms are modified).
  !
  j = 0
  do i = 1, clust%atom_num
     if (clust%mvat(i) == 1) then
        j = j + 1
        clust%xatm(i) = move%rcur(1,j)
        clust%yatm(i) = move%rcur(2,j)
        clust%zatm(i) = move%rcur(3,j)
     endif
  enddo

  if (pbc%is_on) then
     !
     ! If this is a periodic system:
     ! It is possible that some of the atoms have "moved outside of
     ! the cell" because of the minimization. For the electronic
     ! structure problem, we need to "wrap" them back into the
     ! supercell.
     !
     call back_to_cell(pbc,clust%atom_num, &
          clust%xatm,clust%yatm,clust%zatm)
  endif
  !
  ! Save temporary information for a relaxation restart.
  !
! #ifdef AJB_DEBUG
      write(7,*) '       writing stuff to files'
! #endif
  open(45,file='relax_restart.dat',form='formatted')
  write(45,*) ' Atomic coordinates follow (all atoms)'
  do j=1,clust%atom_num
     write(45,*) clust%xatm(j),clust%yatm(j),clust%zatm(j)
  enddo
  write(45,*) 'Iteration number, number of movable atoms'
  write(45,*) move%num,nmovable
  if (move%name == BFGS) then
     ! if BFGS is done, write out info for a BFGS restart.
     write(45,*) ' LBFGS work variables follow'
     write(45,*) move%mmax,move%stepmax
     write(45,'(60a)') move%cbfgs(1:60)
     write(45,'(60a)') move%cbfgs(61:120)
     do j=1,4
        write(45,*) move%lbfgs(j)
     enddo
     do j=1,3*nmovable
        write(45,*) move%lower(j)
        write(45,*) move%upper(j)
     enddo
     do j=1,44+9*nmovable
        write(45,*) move%ibfgs(j)
     enddo
     do j=1,size(move%wbfgs)
        write(45,*) move%wbfgs(j)
     enddo
     do j=1,29
        write(45,*) move%dbfgs(j)
     enddo
! not flushing because file is being closed now
!     call myflush(45) !or is it 46?
     close(45)
! #ifdef AJB_DEBUG
      write(7,*) '       finished writing stuff to relax_restart '
! #endif
  endif
  !
  ! Report new coordinates to file, unless this was the last step.
  !
  if (move%num < move%mxmove) then
     write(66,32) move%num+1
     j=0
     do ity=1,clust%type_num
        write(66,*) clust%name(ity),'    ',clust%natmi(ity)
        do i=1,clust%natmi(ity)
           j = j + 1
           write(66,34) clust%xatm(j), clust%yatm(j), clust%zatm(j)
        enddo
     enddo
! #ifdef AJB_DEBUG
      write(7,*) '       finished writing stuff to atom_cor '
! #endif
  call myflush(66)
  endif

32 format('Step #',1x,i3)
34 format(3(3x,f10.6))

  deallocate(mforce)

end subroutine domove
!===============================================================
