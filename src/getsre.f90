!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates and reports the regular and
! charge-weighted self-consistent residual error (SRE) following
! each iteration. The charge weighted SRE is the one used to
! determine convergence.
!
!---------------------------------------------------------------
subroutine getsre(elec_st,pot,parallel,hcub,imove,iter)

  use constants
  use electronic_struct_module
  use potential_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! potential related data
  type (potential), intent(in) :: pot
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! (grid spacing)**3
  real(dp), intent(in) :: hcub
  ! movement number, self-consistent iteration number
  integer, intent(in) :: imove, iter
  !
  ! Work variables:
  !
  ! "regular" SRE
  real(dp) :: sre(elec_st%nspin)
  ! temporary data holder
  real(dp) :: dvsum, dv
  ! counters
  integer i,isp,jj
  ! spin identifier for printing
  character(len=2) idsp

  !---------------------------------------------------------------

  ! Sum potential differences squared - as is, and weighted by
  ! the charge density.
  sre(:) = zero
  elec_st%sre(:) = zero

  ! Compute SRE.
  do isp = 1, elec_st%nspin
     jj=isp-1+elec_st%nspin
     do i = 1, parallel%mydim
        dv = pot%vnew(i,isp) - pot%vold(i,isp)
        dvsum = dv*dv
        sre(isp) = sre(isp) + dvsum
    !OS + AJB:
    !if rho is zero and you have a lot of vacuum, 
    !it makes the edges not converge at all?
        elec_st%sre(isp) = elec_st%sre(isp) + &
             elec_st%rho(i,jj)*dvsum*real(elec_st%nrep,dp)
     enddo
  enddo
  call psum(sre,elec_st%nspin,parallel%group_size, &
       parallel%group_comm)
  call psum(elec_st%sre,elec_st%nspin, &
       parallel%group_size,parallel%group_comm)

  ! Normalize by h^3 for the "regular" SRE, by h^3/(# of electrons)
  ! for the charge.weighted SRE. If there are no electrons in this
  ! spin channel, set charge-weighted SRE to zero.
  do isp = 1, elec_st%nspin
     sre(isp) = sqrt(sre(isp)*hcub)
     if (elec_st%totel(isp) == zero) then
        elec_st%sre(isp) = zero
     else
        elec_st%sre(isp) = sqrt(hcub*elec_st%sre(isp)/elec_st%totel(isp))
     endif
     ! Set spin identifier.
     jj=isp-1+elec_st%nspin
     if (jj == 1) idsp = '  '
     if (jj == 2) idsp = 'up'
     if (jj == 3) idsp = 'dn'

     if (parallel%iammaster) then
        write(7,*)
        write(7,10) imove, iter, idsp, sre(isp), elec_st%sre(isp)
        call myflush(7)
     endif
  enddo
10 format(i3,'-',i3,1x,a2,' SRE of pot. & charge weighted pot = ', &
        f14.10,2x, f14.10,/)

   !also keep the non weighted sre
    elec_st%plain_sre=sre 

end subroutine getsre
!===============================================================
