!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Save DFT eigenvalues to unit # funit
!
! Author: O. Guliamov, L. Kronik (2004)
!
!---------------------------------------------------------------
subroutine eigensave(elec_st,imove,funit,outevflag)

  use constants
  use electronic_struct_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(in) :: elec_st

  ! current itteration number
  integer, intent(in) :: imove
  ! unit of output file
  integer, intent(in) :: funit
  ! output flag
  integer, intent(in) :: outevflag
  !
  ! Work variables:
  !
  ! counters
  integer i, isp, jj
  ! number of eigenvalues to be written
  integer nstate
  ! representation index
  integer irp
  ! jrep keeps track of how many eigenstates are already in each
  ! representation
  integer :: jrep(elec_st%nrep)
  ! spin identifier for printing
  character (len=2) :: idsp
  integer nspin
  integer kpnum, kplp
  !---------------------------------------------------------------

  kpnum = max(elec_st%nkpt,1)

  nspin = elec_st%nspin/elec_st%mxwd
  do isp = 1, nspin
     jj = isp - 1 + nspin
     if (jj == 1) idsp = '  '
     if (jj == 2) idsp = 'up'
     if (jj == 3) idsp = 'dn'
     if (mod(outevflag,2) == 0) then
        nstate = elec_st%nstate
     else
        nstate = elec_st%ifmax(isp)
     endif
     do kplp = 1, kpnum
        jrep = 0
        do i = 1,nstate
           irp = elec_st%irep(i,kplp,isp)
           jrep(irp) = jrep(irp) + 1
           if (elec_st%nkpt == 0) then
              write(funit,11) i,elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   elec_st%eig(irp,kplp,isp)%en(jrep(irp))*rydberg,  &
                   elec_st%eig(irp,kplp,isp)%occ(jrep(irp)), &
                   irp,imove, idsp 
           else
              write(funit,12) i,elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   elec_st%eig(irp,kplp,isp)%en(jrep(irp))*rydberg,  &
                   elec_st%eig(irp,kplp,isp)%occ(jrep(irp)), &
                   irp,kplp, imove, idsp 
           endif
        enddo
     enddo
  enddo
  call myflush(funit)
11 format(i5,2x,f18.10,2x,f12.4,2x,f12.4,2x,i2,2x,i5,2x,a2)
12 format(i5,2x,f18.10,2x,f12.4,2x,f12.4,2x,i2,2x,i5,2x,i5,2x,a2)

end subroutine eigensave
!===============================================================
