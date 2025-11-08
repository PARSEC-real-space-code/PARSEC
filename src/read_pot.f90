subroutine read_eam_param(nty,Dm,am,rm,ar,br,Fe,ge,hg,rcut,rr,iunit)

  use constants  

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,1), &
                           am(nty,nty,1), &
                           rm(nty,nty,1)
  real(dp), intent(out) :: ar(nty),br(nty) 
  real(dp), intent(out) :: Fe(nty),ge(nty) 
  real(dp), intent(out) :: rr(nty)
  real(dp), intent(out) :: hg,rcut 
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo

 !!read global r_morse parameters
 !read(iunit,*) !blank
 !read(iunit,*) !global N
 !!number of interaction potential
 !do i = 1, nty
 !   do j = 1, nty
 !      if (nty == 1) then
 !         read(iunit,*) dummy, rm(i,j,1)
 !      else if(i .lt. j) then
 !         read(iunit,*) dummy, rm(i,j,1)
 !      endif
 !   end do
 !end do

  !read other parameters
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rcut!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           !read(iunit,*) !h 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rcut!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           !read(iunit,*) !h 
           write(*,*) Dm(i,j,1), am(i,j,1) 
        endif
     end do
  end do
  do i = 1, nty
     read(iunit,*) 
     read(iunit,*) !type
     read(iunit,*) !cutoff
     read(iunit,*) !#rmin
     read(iunit,*) dummy,ar(i)
     read(iunit,*) dummy,br(i)
     read(iunit,*) dummy,rr(i) 
     !read(iunit,*) !h 
  end do
  do i = 1, nty
     read(iunit,*) 
     read(iunit,*) !type
     read(iunit,*) !cutoff
     read(iunit,*) !#rmin
     read(iunit,*) dummy,Fe(i)
     read(iunit,*) dummy,ge(i)
     read(iunit,*) !F1_emb 
  end do

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia,1) = Dm (ia, ja,1)
        am (ja, ia,1) = am (ia, ja,1)
        rm (ja, ia,1) = rm (ia, ja,1)
     end do
  end do

end subroutine read_eam_param
!===============================================================
!
! This subroutine reads the classical potential parameters from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine read_morse_param(mol_dynamic,nty,Dm,am,rm,rc,hs,iunit)

  use constants  
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,1), &
                           am(nty,nty,1), &
                           rm(nty,nty,1), &
                           hs(nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
           write(*,*) Dm(i,j,1), am(i,j,1) 
        endif
     end do
  end do

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia,1) = Dm (ia, ja,1)
        am (ja, ia,1) = am (ia, ja,1)
        rm (ja, ia,1) = rm (ia, ja,1)
        hs (ja, ia  ) = hs (ia, ja  )
     end do
  end do

end subroutine read_morse_param
!===============================================================
!
! This subroutine reads the classical potential parameters from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine read_doublemorse_param(mol_dynamic,nty,Dm,am,rm,delta,rc,hs,iunit)

  use constants  
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,2), &
                           am(nty,nty,2), &
                           rm(nty,nty,2), &
                           delta(nty,nty), &
                           hs(nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,Dm(i,j,2)
           read(iunit,*) dummy,am(i,j,2)
           read(iunit,*) dummy,rm(i,j,2)
           read(iunit,*) dummy,delta(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,Dm(i,j,2)
           read(iunit,*) dummy,am(i,j,2)
           read(iunit,*) dummy,rm(i,j,2)
           read(iunit,*) dummy,delta(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        endif
     end do
  end do

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia, :) = Dm (ia, ja, :)
        am (ja, ia, :) = am (ia, ja, :)
        rm (ja, ia, :) = rm (ia, ja, :)
        delta(ja, ia) = delta(ia, ja)
        hs(ja, ia) = hs(ia, ja)
     end do
  end do

end subroutine read_doublemorse_param
!
subroutine read_ljmorse_param(mol_dynamic,nty,Dm,am,rm,elj,slj,delta,rc,hs,iunit)

  use constants
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic
  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,1), &
                           am(nty,nty,1), &
                           rm(nty,nty,1), &
                           delta(nty,nty), &
                           elj(nty,nty), &
                           slj(nty,nty), &
                           hs(nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,elj(i,j)
           read(iunit,*) dummy,slj(i,j)
           read(iunit,*) dummy,delta(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,elj(i,j)
           read(iunit,*) dummy,slj(i,j)
           read(iunit,*) dummy,delta(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        endif
     end do
  end do

  ! set parameters for opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia, :) = Dm (ia, ja, :)
        am (ja, ia, :) = am (ia, ja, :)
        rm (ja, ia, :) = rm (ia, ja, :)
        elj(ja, ia) = elj(ia, ja)
        slj(ja, ia) = slj(ia, ja)
        delta(ja, ia) = delta(ia, ja)
        hs(ja, ia) = hs(ia, ja)
     end do
  end do

end subroutine read_ljmorse_param
!===============================================================
!
! This subroutine reads the classical potential parameters from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine read_mtheta_param(mol_dynamic,nty,Dm,am,rm,rc,hs,gm,lm,iunit)
  use constants  
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,1), &
                           am(nty,nty,1), &
                           rm(nty,nty,1), &
                           hs(nty,nty), &
                           gm(nty,nty,nty), &
                           lm(nty,nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j, k
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  write(*,*) "Reading Morse-theta potential params" 
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
           write(*,*) Dm(i,j,1), am(i,j,1) 
        endif
     end do
  end do

  write(*,*) "Read Gammas" 
  !Read Gamma
  read(iunit,*)
  read(iunit,*) !type
  read(iunit,*) dummy, rc!cutoff
  read(iunit,*) !#rmin
  do i = 1, nty
    do j = 1, nty
      do k = j, nty
           read(iunit,*) dummy,gm(i,j,k)
      enddo
    enddo
  enddo

  write(*,*) "Read Lambdas" 
  !Read Lambda
  read(iunit,*)
  read(iunit,*) !type
  read(iunit,*) dummy, rc!cutoff
  read(iunit,*) !#rmin
  do i = 1, nty
    do j = 1, nty
      do k = j, nty
           read(iunit,*) dummy,lm(i,j,k)
      enddo
    enddo
  enddo

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia,1) = Dm (ia, ja,1)
        am (ja, ia,1) = am (ia, ja,1)
        rm (ja, ia,1) = rm (ia, ja,1)
        hs (ja, ia  ) = hs (ia, ja  )
     end do
  end do
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           gm(i,k,j) = gm(i,j,k)
           lm(i,k,j) = lm(i,j,k)
      enddo
    enddo
  enddo
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           write(*,*) "ijk and gamma"
           write(*,'(3i, 1f15.10)') i,j,k,gm(i,j,k)
           write(*,'(3i, 1f15.10)') i,j,k,lm(i,j,k)
      enddo
    enddo
  enddo


end subroutine read_mtheta_param
!===============================================================
!
! This subroutine reads the classical potential parameters from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine read_ljmt_param(mol_dynamic,nty,Dm,am,rm,elj,slj,rc,hs,gm,lm,iunit)
  use constants  
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: Dm(nty,nty,1), &
                           am(nty,nty,1), &
                           rm(nty,nty,1), &
                           elj(nty,nty), &
                           slj(nty,nty), &
                           hs(nty,nty), &
                           gm(nty,nty,nty), &
                           lm(nty,nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j, k
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  write(*,*) "Reading Morse-theta potential params" 
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,elj(i,j)
           read(iunit,*) dummy,slj(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,Dm(i,j,1)
           read(iunit,*) dummy,am(i,j,1)
           read(iunit,*) dummy,rm(i,j,1)
           read(iunit,*) dummy,elj(i,j)
           read(iunit,*) dummy,slj(i,j)
           if(mol_dynamic%sc) read(iunit,*) dummy,hs(i,j) 
           write(*,*) Dm(i,j,1), am(i,j,1) 
        endif
     end do
  end do

  write(*,*) "Read Gammas" 
  !Read Gamma
  read(iunit,*)
  read(iunit,*) !type
  read(iunit,*) dummy, rc!cutoff
  read(iunit,*) !#rmin
  do i = 1, nty
    do j = 1, nty
      do k = j, nty
           read(iunit,*) dummy,gm(i,j,k)
      enddo
    enddo
  enddo

  write(*,*) "Read Lambdas" 
  !Read Lambda
  read(iunit,*)
  read(iunit,*) !type
  read(iunit,*) dummy, rc!cutoff
  read(iunit,*) !#rmin
  do i = 1, nty
    do j = 1, nty
      do k = j, nty
           read(iunit,*) dummy,lm(i,j,k)
      enddo
    enddo
  enddo

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        Dm (ja, ia,1) = Dm (ia, ja,1)
        am (ja, ia,1) = am (ia, ja,1)
        rm (ja, ia,1) = rm (ia, ja,1)
        elj(ja, ia) = elj (ia, ja)
        slj(ja, ia) = slj (ia, ja)
        hs (ja, ia  ) = hs (ia, ja  )
     end do
  end do
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           gm(i,k,j) = gm(i,j,k)
           lm(i,k,j) = lm(i,j,k)
      enddo
    enddo
  enddo
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           write(*,*) "ijk and gamma"
           write(*,'(3i, 1f15.10)') i,j,k,gm(i,j,k)
           write(*,'(3i, 1f15.10)') i,j,k,lm(i,j,k)
      enddo
    enddo
  enddo


end subroutine read_ljmt_param
!===============================================================
!
! This subroutine reads the classical potential parameters from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine read_sw_param(mol_dynamic,nty,sw2,sw3,swa,rc,iunit)
  use constants  
  use molecular_dynamic_module

  type (molecular_dynamic), intent(in) :: mol_dynamic

  integer , intent(in) :: nty
  ! num type
  ! Morse parameters
  real(dp), intent(out) :: sw2(nty,nty,6), &
                           sw3(nty,nty,2), &
                           swa(nty,nty,nty), &
                           rc
  integer, intent(in) :: iunit
  ! counter
  integer :: i, j, k
  ! dummy
  character :: dummy
  ! n_interaction
  integer :: nsum
  
  hg = 0

  ! read paramters from pot file
  !skip header
  do i = 1, 6
    read(iunit,*)
  enddo
  !read global r_morse parameters
  !read(iunit,*) !global N
  ! number of interaction potential
  !do i = 1, nty
  !   do j = 1, nty
  !      if (nty == 1) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !      else if(i .lt. j) then
  !         read(iunit,*) dummy, rm(i,j)
  !         read(iunit,*) dummy, hg ! global smooth cutoff param
  !         write(*,*) rm(i,j), hg 
  !      endif
  !   end do
  !end do
  !read other parameters
  write(*,*) "Reading stiweb potential params" 
  do i = 1, nty
     do j = 1, nty
        if(nty == 1) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,sw2(i,j,1)
           read(iunit,*) dummy,sw2(i,j,2)
           read(iunit,*) dummy,sw2(i,j,3)
           read(iunit,*) dummy,sw2(i,j,4)
           read(iunit,*) dummy,sw2(i,j,5)
           read(iunit,*) dummy,sw2(i,j,6)
           if(mol_dynamic%sc) read(iunit,*) dummy
        else if (i .le. j ) then
           read(iunit,*)
           read(iunit,*) !type
           read(iunit,*) dummy, rc!cutoff
           read(iunit,*) !#rmin
           read(iunit,*) dummy,sw2(i,j,1)
           read(iunit,*) dummy,sw2(i,j,2)
           read(iunit,*) dummy,sw2(i,j,3)
           read(iunit,*) dummy,sw2(i,j,4)
           read(iunit,*) dummy,sw2(i,j,5)
           read(iunit,*) dummy,sw2(i,j,6)
           if(mol_dynamic%sc) read(iunit,*) dummy
        endif
     end do
  end do

  do i = 1, nty
    do j = 1, nty
        if(nty == 1) then
          read(iunit,*)
          read(iunit,*) !type
          read(iunit,*) dummy, rc!cutoff
          read(iunit,*) !#rmin
          read(iunit,*) dummy,sw3(i,j,1)
          read(iunit,*) dummy,sw3(i,j,2)
          write(*,*) sw3(i,j,:)
        else if (i .le. j ) then
          read(iunit,*)
          read(iunit,*) !type
          read(iunit,*) dummy, rc!cutoff
          read(iunit,*) !#rmin
          read(iunit,*) dummy,sw3(i,j,1)
          read(iunit,*) dummy,sw3(i,j,2)
          write(*,*) sw3(i,j,:)
        endif
    enddo
  enddo

  write(*,*) "Read Lambdas" 
  !Read Lambda
  read(iunit,*)
  read(iunit,*) !type
  read(iunit,*) dummy, rc!cutoff
  read(iunit,*) !#rmin
  do i = 1, nty
    do j = 1, nty
      do k = j, nty
           read(iunit,*) dummy,swa(i,j,k)
           write(*,*) "Lambdas" , i, j, k
      enddo
    enddo
  enddo

  ! set parameters to opposite combinations of ia and ja
  do ia = 1, nty
     do ja = ia, nty
        sw2(ja, ia, :) = sw2(ia, ja, :)
        sw3(ja, ia, :) = sw3(ia, ja, :)
     end do
  end do
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           swa(i,k,j) = swa(i,j,k)
      enddo
    enddo
  enddo
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           write(*,*) "ijk and gamma"
           write(*,'(3i, 1f15.10)') i,j,k,swa(i,j,k)
      enddo
    enddo
  enddo


end subroutine read_sw_param
