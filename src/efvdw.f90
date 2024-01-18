subroutine efvdw(clust, elec_st, grid, parallel, pbc, vdw_forces)

use cluster_module
use electronic_struct_module
use grid_module
use parallel_data_module
use pbc_module 
implicit none

type(cluster), intent(in) :: clust
type(electronic_struct), intent(inout) :: elec_st
type(grid_data), intent(in) :: grid
type(parallel_data), intent(in) :: parallel
type(pbc_data), intent(in) :: pbc
real (dp), intent(inout) :: vdw_forces(clust%atom_num,3)


integer :: atmn,ja,itype,iat,ieff,jeff,n_cell,i1,i2,i3,atmn1,atmn2,i_per
integer :: num_cells_super 

! denominator of eq. (8) 
real (dp) :: srhof(parallel%mydim,1)

! hirshfeld atomic partitioning weight for the atoms
real (dp) :: watm(parallel%mydim,clust%atom_num)

! the integrand of the numertor and denominator of eq. (7) 
real (dp) :: v_d(parallel%mydim,clust%atom_num)
real (dp) :: v_u(parallel%mydim,clust%atom_num)
real (dp) :: v_d_pbc(parallel%mydim,clust%atom_num)
real (dp) :: v_u_pbc(parallel%mydim,clust%atom_num)

! free and effective volume
real (dp) :: v_free(clust%atom_num),v_eff(clust%atom_num)
real (dp) :: v_free_pbc(clust%atom_num),v_eff_pbc(clust%atom_num)

! v_eff/v_free
real (dp) :: v_eff_free(clust%atom_num)

! c6, alpha and R0 free
real (dp) :: c6ii_free(clust%atom_num),alpha_free(clust%atom_num),R0(clust%atom_num)

! c6, alpha and R eff
real (dp) :: c6ii_eff(clust%atom_num),alpha_eff(clust%atom_num),R_eff(clust%atom_num)

! c6ij eff
real (dp) :: c6ij_eff

! distance between pairs of atoms
real (dp) :: R_cor

! Reff_i+Reff_j
real (dp) :: R_eff_ij

! f_damp
real (dp) :: f_damp

! the van der waals energy from each pair of atoms
real (dp) :: evdw_add

! forces between pairs of atoms in x y and z directions
real (dp) :: Fx_ij(clust%atom_num,1),Fy_ij(clust%atom_num,1),Fz_ij(clust%atom_num,1)

! total force on each atom in x y and z directions
real (dp) :: Ft_xyz(clust%atom_num,3)

! storge variable
real (dp) :: rr0,rexp,aexp,rfac,ffac,rphi,maxchange
real(dp) :: evdw_pbc,evdw_pbc_add,Ft_xyz_pbc(clust%atom_num,3),F_xyz_pbc(clust%atom_num,3)
real(dp) :: C6_pbc
real(dp) :: R_pbc
real(dp) :: coord_rep_x(clust%atom_num),coord_rep_y(clust%atom_num),coord_rep_z(clust%atom_num),coord_diff(3,1)
real(dp) :: d2_pbc,d6_pbc,d_pbc,exp_arg,exp_pbc,f_dump_pbc
logical pbc_conv

! --------------------------------------------------------------

! what is this	if (parallel%iammaster) open(1,file='efvdw_data.dat')
v_eff(:)=zero


v_free(:)=elec_st%v_free_periodic(1,:)



num_cells_super = (elec_st%hirsh_3d_cell(1,1)*2+1)* &
                  (elec_st%hirsh_3d_cell(1,2)*2+1)* &
                  (elec_st%hirsh_3d_cell(1,3)*2+1)

do ja = 1, clust%atom_num
do i_per = 1, num_cells_super
        v_eff(ja) = v_eff(ja) + &
                sum(elec_st%dist_periodic(:,i_per,ja)**3D0* &
                elec_st%hirshfeld_periodic(:,i_per,ja)* &
                elec_st%rho(:,1)*grid%hcub)
enddo
enddo



call psum(v_eff,clust%atom_num,parallel%group_size,parallel%group_comm)


v_eff_free=v_eff/v_free


if (.not. parallel%iammaster) return
! okay?
ja=0
do itype=1,clust%type_num
  do iat=1,clust%natmi(itype)
    ja=ja+1
!        select case(clust%name(itype)) !<-TODO
    if(clust%name(itype)=='H') then
        alpha_free(ja)=4.5D0
        c6ii_free(ja)=6.5D0
        R0(ja)=3.1D0
    else if(clust%name(itype)=='He') then 
        alpha_free(ja)=1.38D0
        c6ii_free(ja)=1.46D0
        R0(ja)=2.65D0
    else if(clust%name(itype)=='Li') then 
        alpha_free(ja)=164.2D0
        c6ii_free(ja)=1387.D0
        R0(ja)=4.16D0
    else if(clust%name(itype)=='Be') then 
        alpha_free(ja)=38.D0
        c6ii_free(ja)=214.D0
        R0(ja)=4.17D0
    else if(clust%name(itype)=='B') then 
        alpha_free(ja)=21.D0
        c6ii_free(ja)=99.5D0
        R0(ja)=3.89D0
    else if(clust%name(itype)=='C') then 
        alpha_free(ja)=12.D0
        c6ii_free(ja)=46.6D0
        R0(ja)=3.59D0
    else if(clust%name(itype)=='N') then 
        alpha_free(ja)=7.4D0
        c6ii_free(ja)=24.2D0
        R0(ja)=3.34D0
    else if(clust%name(itype)=='O') then 
        alpha_free(ja)=5.4D0
        c6ii_free(ja)=15.6D0
        R0(ja)=3.19D0
    else if(clust%name(itype)=='F') then 
        alpha_free(ja)=3.8D0
        c6ii_free(ja)=9.52D0
        R0(ja)=3.04D0
    else if(clust%name(itype)=='Ne') then 
        alpha_free(ja)=2.67D0
        c6ii_free(ja)=6.38D0
        R0(ja)=2.91D0
    else if(clust%name(itype)=='Na') then 
        alpha_free(ja)=162.7D0
        c6ii_free(ja)=1556.D0
        R0(ja)=3.73D0
    else if(clust%name(itype)=='Mg') then 
        alpha_free(ja)=71.D0
        c6ii_free(ja)=627.D0
        R0(ja)=4.27D0
    else if(clust%name(itype)=='Al') then 
        alpha_free(ja)=60.D0
        c6ii_free(ja)=528.D0
        R0(ja)=4.33D0
    else if(clust%name(itype)=='Si') then 
        alpha_free(ja)=37.D0
        c6ii_free(ja)=305.D0
        R0(ja)=4.2D0
    else if(clust%name(itype)=='P') then 
        alpha_free(ja)=25.D0
        c6ii_free(ja)=185.D0
        R0(ja)=4.01D0
    else if(clust%name(itype)=='S') then 
        alpha_free(ja)=19.6D0
        c6ii_free(ja)=134.D0
        R0(ja)=3.86D0
    else if(clust%name(itype)=='Cl') then 
        alpha_free(ja)=15.D0
        c6ii_free(ja)=94.6D0
        R0(ja)=3.71D0
    else if(clust%name(itype)=='Ar') then 
        alpha_free(ja)=11.1D0
        c6ii_free(ja)=64.3D0
        R0(ja)=3.55D0
    else if(clust%name(itype)=='K') then 
        alpha_free(ja)=292.9D0
        c6ii_free(ja)=3897.D0
        R0(ja)=3.71D0
    else if(clust%name(itype)=='Ca') then 
        alpha_free(ja)=160.D0
        c6ii_free(ja)=2221.D0
        R0(ja)=4.65D0
    else if(clust%name(itype)=='Sc') then 
        alpha_free(ja)=120.D0
        c6ii_free(ja)=1383.D0
        R0(ja)=4.59D0
    else if(clust%name(itype)=='Ti') then 
        alpha_free(ja)=98.D0
        c6ii_free(ja)=1044.D0
        R0(ja)=4.51D0
    else if(clust%name(itype)=='V') then 
        alpha_free(ja)=84.D0
        c6ii_free(ja)=832.D0
        R0(ja)=4.44D0
    else if(clust%name(itype)=='Cr') then 
        alpha_free(ja)=78.D0
        c6ii_free(ja)=602.D0
        R0(ja)=3.99D0
    else if(clust%name(itype)=='Mn') then 
        alpha_free(ja)=63.D0
        c6ii_free(ja)=552.D0
        R0(ja)=3.97D0
    else if(clust%name(itype)=='Fe') then 
        alpha_free(ja)=56.D0
        c6ii_free(ja)=482.D0
        R0(ja)=4.23D0
    else if(clust%name(itype)=='Co') then 
        alpha_free(ja)=50.D0
        c6ii_free(ja)=408.D0
        R0(ja)=4.18D0
    else if(clust%name(itype)=='Ni') then 
        alpha_free(ja)=48.D0
        c6ii_free(ja)=373.D0
        R0(ja)=3.82D0
    else if(clust%name(itype)=='Cu') then 
        alpha_free(ja)=42.D0
        c6ii_free(ja)=253.D0
        R0(ja)=3.76D0
    else if(clust%name(itype)=='Zn') then 
        alpha_free(ja)=40.D0
        c6ii_free(ja)=284.D0
        R0(ja)=4.02D0
    else if(clust%name(itype)=='Ga') then 
        alpha_free(ja)=60.D0
        c6ii_free(ja)=498.D0
        R0(ja)=4.19D0
    else if(clust%name(itype)=='Ge') then 
        alpha_free(ja)=41.D0
        c6ii_free(ja)=354.D0
        R0(ja)=4.2D0
    else if(clust%name(itype)=='As') then 
        alpha_free(ja)=29.D0
        c6ii_free(ja)=246.D0
        R0(ja)=4.11D0
    else if(clust%name(itype)=='Se') then 
        alpha_free(ja)=25.D0
        c6ii_free(ja)=210.D0
        R0(ja)=4.04D0
    else if(clust%name(itype)=='Br') then 
        alpha_free(ja)=20.D0
        c6ii_free(ja)=162.D0
        R0(ja)=3.93D0
    else if(clust%name(itype)=='Kr') then 
        alpha_free(ja)=16.8D0
        c6ii_free(ja)=129.6D0
        R0(ja)=3.82D0
    else if(clust%name(itype)=='Rb') then 
        alpha_free(ja)=319.2D0
        c6ii_free(ja)=4691.D0
        R0(ja)=3.72D0
    else if(clust%name(itype)=='Sr') then 
        alpha_free(ja)=199.D0
        c6ii_free(ja)=3170.D0
        R0(ja)=4.54D0
    else if(clust%name(itype)=='Pd') then 
        alpha_free(ja)=23.68D0
        c6ii_free(ja)=157.5D0
        R0(ja)=3.66D0
    else if(clust%name(itype)=='Ag') then 
        alpha_free(ja)=50.6D0
        c6ii_free(ja)=339.D0
        R0(ja)=3.82D0
    else if(clust%name(itype)=='I') then 
        alpha_free(ja)=35.D0
        c6ii_free(ja)=385.D0
        R0(ja)=4.17D0
    else if(clust%name(itype)=='Xe') then 
        alpha_free(ja)=27.3D0
        c6ii_free(ja)=285.9D0
        R0(ja)=4.08D0
    else if(clust%name(itype)=='Ba') then 
        alpha_free(ja)=275.D0
        c6ii_free(ja)=5727.D0
        R0(ja)=4.77D0
    else if(clust%name(itype)=='Ir') then 
        alpha_free(ja)=42.51D0
        c6ii_free(ja)=359.1D0
        R0(ja)=4.D0
    else if(clust%name(itype)=='Pt') then 
        alpha_free(ja)=39.68D0
        c6ii_free(ja)=347.1D0
        R0(ja)=3.92D0
    else if(clust%name(itype)=='Au') then 
        alpha_free(ja)=36.5D0
        c6ii_free(ja)=298.D0
        R0(ja)=3.86D0
    else if(clust%name(itype)=='Hg') then 
        alpha_free(ja)=33.9D0
        c6ii_free(ja)=392.D0
        R0(ja)=3.98D0
    else if(clust%name(itype)=='Pb') then 
        alpha_free(ja)=61.8D0
        c6ii_free(ja)=697.D0
        R0(ja)=4.31D0
    else if(clust%name(itype)=='Bi') then 
        alpha_free(ja)=49.02D0
        c6ii_free(ja)=571.D0
        R0(ja)=4.32D0
    else
        !write a failing case! case default
    end if
  end do
end do


c6ii_eff=(v_eff_free**two)*c6ii_free


alpha_eff=v_eff_free*alpha_free


R_eff=(v_eff_free**(one/three))*R0


elec_st%evdw=zero

Fx_ij(:,:)=zero
Fy_ij(:,:)=zero
Fz_ij(:,:)=zero

Ft_xyz(:,:)=zero


do ieff=1,clust%atom_num
do jeff=1,clust%atom_num
        if (jeff==ieff) cycle

        c6ij_eff = (two*c6ii_eff(ieff)*c6ii_eff(jeff))/ &
                (alpha_eff(jeff)/alpha_eff(ieff)*c6ii_eff(ieff) + &
                 alpha_eff(ieff)/alpha_eff(jeff)*c6ii_eff(jeff))

        R_cor=dsqrt((clust%xatm(ieff)-clust%xatm(jeff))**two&
                +(clust%yatm(ieff)-clust%yatm(jeff))**two&
                +(clust%zatm(ieff)-clust%zatm(jeff))**two)

        R_eff_ij=R_eff(ieff)+R_eff(jeff)

        f_damp=one/(one+dexp(-20.D0*(R_cor/(0.94D0*R_eff_ij)-one)))

        evdw_add=(f_damp*c6ij_eff)/(R_cor**six)

        elec_st%evdw=elec_st%evdw+evdw_add


        ! forces
        ! Note!! the expresions for the forces were taken from fhi-aims code
        rr0=R_eff_ij*0.94D0

        rexp=dexp(-20.D0*(R_cor/(0.94D0*R_eff_ij)-one))
        aexp=one+rexp
        rphi=-evdw_add

        rfac=-six/R_cor + 20D0*rexp/(aexp*rr0)
        ffac=rphi*rfac/R_cor



        Fx_ij(ieff,1)=Fx_ij(ieff,1)-(clust%xatm(ieff)-clust%xatm(jeff))*ffac
        Fy_ij(ieff,1)=Fy_ij(ieff,1)-(clust%yatm(ieff)-clust%yatm(jeff))*ffac
        Fz_ij(ieff,1)=Fz_ij(ieff,1)-(clust%zatm(ieff)-clust%zatm(jeff))*ffac


enddo 
enddo


Ft_xyz(:,1)=Fx_ij(:,1)
Ft_xyz(:,2)=Fy_ij(:,1)
Ft_xyz(:,3)=Fz_ij(:,1)

! factor 0.5 due to complete sumations
! factor -1 because.
elec_st%evdw=-0.5D0*elec_st%evdw
!factor 2 for converting into Rydberg
elec_st%evdw=two*elec_st%evdw


Ft_xyz_pbc(:,:)=zero        
if (pbc%is_on) then
        pbc_conv=.false.
        n_cell=0

        do while (.not.pbc_conv) 
                n_cell=n_cell+1
                evdw_pbc=0D0
                Ft_xyz_pbc(:,:)=0D0

                do i1=-n_cell,n_cell
                do i2=-n_cell,n_cell
                do i3=-n_cell,n_cell
                        if ( (abs(i1).eq.n_cell) .or. &
                             (abs(i2).eq.n_cell) .or. &
                             (abs(i3).eq.n_cell) ) then

                                do atmn1=1,clust%atom_num,1 
                                do atmn2=1,clust%atom_num,1

C6_pbc=two*c6ii_eff(atmn1)*c6ii_eff(atmn2)/ &
        (alpha_eff(atmn2)/alpha_eff(atmn1)*c6ii_eff(atmn1) + &
         alpha_eff(atmn1)/alpha_eff(atmn2)*c6ii_eff(atmn2))

coord_rep_x(atmn2) = clust%xatm(atmn2) + &
        i1*pbc%latt_vec(1,1)+i2*pbc%latt_vec(1,2)+i3*pbc%latt_vec(1,3)
coord_rep_y(atmn2) = clust%yatm(atmn2) + &
        i1*pbc%latt_vec(2,1)+i2*pbc%latt_vec(2,2)+i3*pbc%latt_vec(2,3)
coord_rep_z(atmn2) = clust%zatm(atmn2) + &
        i1*pbc%latt_vec(3,1)+i2*pbc%latt_vec(3,2)+i3*pbc%latt_vec(3,3)

coord_diff(1,1)=clust%xatm(atmn1)-coord_rep_x(atmn2)
coord_diff(2,1)=clust%yatm(atmn1)-coord_rep_y(atmn2)
coord_diff(3,1)=clust%zatm(atmn1)-coord_rep_z(atmn2)

d2_pbc=coord_diff(1,1)**two+coord_diff(2,1)**two+coord_diff(3,1)**two

d6_pbc=d2_pbc**three

d_pbc=dsqrt(d2_pbc)

exp_arg=-20D0*(d_pbc/(0.94D0*(R_eff(atmn1)+R_eff(atmn2)))-1)
exp_pbc=dexp(exp_arg)
f_dump_pbc=1/(1+exp_pbc)
evdw_pbc_add=-C6_pbc*f_dump_pbc/d6_pbc
evdw_pbc=evdw_pbc+evdw_pbc_add


! forces
! Note!! the expresions for the forces were taken from fhi-aims code
! so someone should double check to see if there were bugfixes there!
rr0=0.94D0*(R_eff(atmn1)+R_eff(atmn2))

rexp=exp_pbc
aexp=one+rexp
rphi=evdw_pbc_add

rfac=-six/d_pbc + 20D0*rexp/(aexp*rr0)
ffac=rphi*rfac/d_pbc


F_xyz_pbc(atmn1,1)=-coord_diff(1,1)*ffac
Ft_xyz_pbc(atmn1,1)=Ft_xyz_pbc(atmn1,1)+F_xyz_pbc(atmn1,1)

F_xyz_pbc(atmn1,2)=-coord_diff(2,1)*ffac
Ft_xyz_pbc(atmn1,2)=Ft_xyz_pbc(atmn1,2)+F_xyz_pbc(atmn1,2)

F_xyz_pbc(atmn1,3)=-coord_diff(3,1)*ffac
Ft_xyz_pbc(atmn1,3)=Ft_xyz_pbc(atmn1,3)+F_xyz_pbc(atmn1,3)

                                enddo
                                enddo
                        endif
                enddo
                enddo
                enddo

                ! factor of 0.5 due to complete sumation        
                evdw_pbc=0.5D0*evdw_pbc
                ! factor of 2 to convert for Rydberg
                evdw_pbc=2D0*evdw_pbc
                pbc_conv=(dabs(evdw_pbc).lt.7.3D-8)
                elec_st%evdw=elec_st%evdw+evdw_pbc
                ! the forces here are in units of Ha/bohr
                Ft_xyz(:,:)=Ft_xyz(:,:)+Ft_xyz_pbc(:,:)
        enddo
endif

vdw_forces(:,:)=zero  
vdw_forces(:,:)=Ft_xyz(:,:)

write(7,*)"NOTE: Writing VDW force info to STDOUT"

write(*,*)"vdW force components in units of eV/A"

do ieff=1,clust%atom_num
        write(*,*) &
                Ft_xyz(ieff,1)*51.42208245D0, &
                Ft_xyz(ieff,2)*51.42208245D0, &
                Ft_xyz(ieff,3)*51.42208245D0
enddo


write(*,*)
write(*,*)"vdW force components in units of Ha/bohr"

do ieff=1,clust%atom_num
        write(*,*)vdw_forces(ieff,1),vdw_forces(ieff,2),vdw_forces(ieff,3)
enddo



write(*,*)
write(*,*)
write(*,*)"Free and effective volumes of each atom"

do ieff=1,clust%atom_num
        write(*,*)v_free(ieff),v_eff(ieff)
enddo

write(*,*)
write(*,*)"Ratio of effective and free volumes"

do ieff=1,clust%atom_num
        write(*,*)v_eff(ieff)/v_free(ieff)
enddo

write(*,*)
write(*,*)



return
end subroutine efvdw
