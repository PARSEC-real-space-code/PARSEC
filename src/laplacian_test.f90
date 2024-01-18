! this subrouting is written especially to test
! the implementation of the laplacian on the 
! non-orthogonal grid.
!
! it is intended just for debugging on a single 
! processor only!!!!!! AMIR

      subroutine laplacian_test(grid,pbc,parallel)

      use constants
      use grid_module
      use pbc_module
      use parallel_data_module
!#ifdef MPI
!     include mpi definitions
!      use mpi
!#endif

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel

!     work variables

      real(dp) :: func(parallel%mydim), lapl(parallel%mydim), &
                  xx(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)
      real(dp) :: u(3),x(3),gv(3),tmpvec(3)
 
      integer i,j
      integer flag,icellx,icelly,icellz

      real(dp) :: tmpmat(3,3),tmpdet,tmptr,tmpf,rx
      real(dp) :: ainv(3,3),ainvdot(3,3)
      character(len=10) ::  f1="sin function", f2="gaussian", f3

!
!     External functions:
!
      real(dp), external :: dnrm2


      func=zero
      xx=zero

      tmpmat=pbc%avec_norm

      write(9,*) "entered laplacian_test function"

      ainv=tmpmat

      call mtrxin(ainv,tmpdet,tmptr)

      ainvdot=matmul(transpose(ainv),ainv)

      flag=0
!     Producing the test function data!

      do i=1,parallel%mydim
        u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
        u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
        u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
        call matvec3('N',tmpmat,u,x)
!        x=matmul(tmpmat,u)

!        func(i)=(x(1))**3
!        xx(i)=6*x(1)
!        func(i)=x(1)**3+x(2)**4
!        xx(i)=6*x(1)+12*x(2)**2

!        preparing a test vector with 
!        arbitrary integers for the sine function.

         if(flag>0) then
            f3=f1
            gv(1)=4
            gv(2)=8
            gv(3)=5
            tmpvec(:)=gv(:)*u(:)/pbc%box_size(:)
            func(i)=sin(2*pi*sum(tmpvec))         

        
!        The analytical form of the laplacian of the sine function.

            tmpvec(:)=(2*pi*gv(:)/pbc%box_size(:))
            call matvec3('N',ainvdot,tmpvec,gv)
            xx(i)=-dot_product(tmpvec,gv)*func(i)
         else
          f3=f2
          do icellx = -1,1
             do icelly = -1,1
                do icellz = -1,1

                   u(1)=(grid%shift(1) + grid%kx(i))* &
                         grid%step(1)+                &
                         real(icellx,dp)*pbc%box_size(1)
                   u(2)=(grid%shift(2) + grid%ky(i))* &
                         grid%step(2)+                &
                         real(icelly,dp)*pbc%box_size(2)
                   u(3)=(grid%shift(3) + grid%kz(i))* &
                         grid%step(3)+                &
                         real(icellz,dp)*pbc%box_size(3)
                   call matvec3('N',tmpmat,u,x)
!                   x=matmul(tmpmat,u)
                   rx=sum(x(:)**2)
                   tmpf=exp(-4*rx)
                   func(i)=func(i)+tmpf
                   xx(i)=xx(i)+(64*rx-24)*tmpf
                 enddo
              enddo
          enddo

         endif

!         func(i)=sin(2*pi*5*u(1)/pbc%box_size(1))
!         xx(i)=-((2*pi*5*ainv(1,1)/pbc%box_size(1))**2+
!     1    (2*pi*5*ainv(1,2)/pbc%box_size(1))**2 +
!     2    (2*pi*5*ainv(1,3)/pbc%box_size(1))**2 ) * func(i)
         
!         func(i)=exp(-(x(1)-3)**2-(x(2)-3)**2)
!         xx(i)=(4*(x(1)-3)**2-2+4*(x(2)-3)**2-2)*func(i)

      end do

      call lapmvs(parallel, func, lapl, grid%coe2, grid%norder, &
                       grid%lap_dir_num, vec)

      j=0
      write(9,*) 'laplacian exception for ',f3,'at:'
      do i=1,parallel%mydim
!         if((abs(lapl(i))-abs(xx(i)))<0.00000000000000000000001) then 
         if(i<0) then 
           j=j+1
         else
           u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
           u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
           u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,103) u, lapl(i), xx(i)
         endif
103   format(5(f15.9,1x))
      end do

      write(9,*) 'number of diffs in laplacian', parallel%mydim-j
      write(9,*) 'total error is:', &
                  real(dnrm2(parallel%mydim,abs(lapl)-abs(xx),1),dp)

      end subroutine laplacian_test

      subroutine lap_symm_check(grid,pbc,parallel)

!     This subroutine checks that the laplacian is symmetric
!     It does so by producing 2 random vectors x and y and 
!     checking that (x,Ay)=(y,Ax).
!
      use constants
      use grid_module
      use pbc_module
      use parallel_data_module

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel

      real(dp) :: x(parallel%mydim), laplx(parallel%mydim), &
                  y(parallel%mydim), laply(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)

      real(dp) :: sumx,sumy
      
      integer i

!
!     External functions:
!
      real(dp), external :: dnrm2


      call random_seed()

      do i=1,parallel%mydim
         call random_number(x(i))
         call random_number(y(i))
      enddo

      call lapmvs(parallel, x, laplx, grid%coe2, grid%norder, &
                        grid%lap_dir_num, vec)

      call lapmvs(parallel, y, laply, grid%coe2, grid%norder, &
                        grid%lap_dir_num, vec)

      sumx=dot_product(x,laply)
      sumy=dot_product(y,laplx)

      write(9,*) 'xAy and yAx are:', sumx, sumy
      write(9,*) 'total diff x-y:', &
                  real(dnrm2(parallel%mydim,x-y,1),dp)

      
      end subroutine lap_symm_check




      subroutine gradient_test(grid,pbc,parallel)

      use constants
      use grid_module
      use pbc_module
      use parallel_data_module

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel

!     work variables

      real(dp) :: func(parallel%mydim), gradl(3,parallel%mydim), &
                  xx(3,parallel%mydim),uu(3,parallel%mydim),     &
                  funcd(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)
      real(dp) :: u(3),x(3),gv(3),tmpvec(3)
 
      integer i,j

      real(dp) :: tmpmat(3,3),tmpdet,tmptr
      real(dp) :: ainv(3,3)

      tmpmat=pbc%avec_norm

      ainv=tmpmat

      call mtrxin(ainv,tmpdet,tmptr)

!     Producing the test function data!

      do i=1,parallel%mydim
        u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
        u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
        u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
        call matvec3('N',tmpmat,u,x)
!        x=matmul(tmpmat,u)

!        func(i)=(x(1))**3
!        xx(i)=6*x(1)
!        func(i)=x(1)**3+x(2)**4
!        xx(i)=6*x(1)+12*x(2)**2

!        preparing a test vector with 
!        arbitrary integers for the sine function.

         gv(1)=4
         gv(2)=8
         gv(3)=5
         tmpvec(:)=gv(:)*u(:)/pbc%box_size(:)
         func(i)=sin(2*pi*sum(tmpvec))         

         funcd(i)=cos(2*pi*sum(tmpvec))

         uu(1:3,i)=2*pi*(gv(1:3)/pbc%box_size(1:3))*funcd(i)

!        The analytical form of the laplacian of the sine function.

!         xx(1,i)=ainv(1,1)*uu(1,i)+ainv(2,1)*uu(2,i)+ainv(3,1)*uu(3,i)
!         xx(2,i)=ainv(1,2)*uu(1,i)+ainv(2,2)*uu(2,i)+ainv(3,2)*uu(3,i)
!         xx(3,i)=ainv(1,3)*uu(1,i)+ainv(2,3)*uu(2,i)+ainv(3,3)*uu(3,i)
         call matvec3('N',grid%grad_bvec_norm,uu(1,i),xx(1,i))
!         xx(1:3,i)=matmul(grid%grad_bvec_norm,uu(1:3,i))
!         xx(1,i)=2*pi*(gv(1)*ainv(1,1)/pbc%box_size(1)+
!     1                 gv(2)*ainv(2,1)/pbc%box_size(2)+
!     2                 gv(3)*ainv(3,1)/pbc%box_size(3))*
!     3             cos(2*pi*(sum(tmpvec)))
!
!         xx(2,i)=2*pi*(gv(1)*ainv(1,2)/pbc%box_size(1)+
!     1                 gv(2)*ainv(2,2)/pbc%box_size(2)+
!     2                 gv(3)*ainv(3,2)/pbc%box_size(3))*
!     3             cos(2*pi*(sum(tmpvec)))
!
!         xx(3,i)=2*pi*(gv(1)*ainv(1,3)/pbc%box_size(1)+
!     1                 gv(2)*ainv(2,3)/pbc%box_size(2)+
!     2                 gv(3)*ainv(3,3)/pbc%box_size(3))*
!     3             cos(2*pi*(sum(tmpvec)))

         

      end do

      call gradmvs(grid,parallel, func, gradl, grid%coe1, grid%norder, &
                   vec)

      j=0
      write(9,*) 'gradient exception at:'
      do i=1,parallel%mydim
         if((abs(gradl(2,i))-abs(xx(2,i)))<0.001) then 
           j=j+1
         else
           u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
           u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
           u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,103) u, gradl(1,i), xx(1,i)
         endif
103   format(5(f15.9,1x))
      end do

      end subroutine gradient_test

      subroutine hpotcg_test(grid,pbc,parallel,elec_st)

      use constants
      use grid_module
      use pbc_module
      use parallel_data_module
      use electronic_struct_module

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel
      type (electronic_struct), intent(in) :: elec_st

!     work variables

      real(dp) :: pot_func(parallel%mydim), lapl(parallel%mydim), &
                  rho(parallel%mydim), pot_orig(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)
      real(dp) :: u(3),x(3),gv(3),tmpvec(3),rx

      integer i,j,flag
      integer icellx,icelly,icellz

      real(dp) :: tmpmat(3,3),tmpdet,tmptr,rhoav
      real(dp) :: ainv(3,3),ainvdot(3,3)

      write(9,*) "entered hpotcg_test function" 

      tmpmat=pbc%avec_norm

      ainv=tmpmat

      call mtrxin(ainv,tmpdet,tmptr)

      ainvdot=matmul(transpose(ainv),ainv)

!     Producing the test function data!

      pot_func=0
      rho=0    
      flag=0

      do i=1,parallel%mydim

        if(flag>0) then
          u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
          u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
          u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
          call matvec3('N',tmpmat,u,x)
!          x=matmul(tmpmat,u)
          gv(1)=4
          gv(2)=8
          gv(3)=5
          tmpvec(:)=gv(:)*u(:)/pbc%box_size(:)

          pot_func(i)=sin(2*pi*sum(tmpvec))


!        The analytical form of the laplacian of the sine function.

          tmpvec(:)=2*pi*gv(:)/pbc%box_size(:)
          call matvec3('N',ainvdot,tmpvec,gv)
          rho(i)=-dot_product(tmpvec,gv)*pot_func(i)
        else
          do icellx = -1,1
             do icelly = -1,1
                do icellz = -1,1

                   u(1)=(grid%shift(1) + grid%kx(i))* &
                         grid%step(1)+ &
                         real(icellx,dp)*pbc%box_size(1)
                   u(2)=(grid%shift(2) + grid%ky(i))* &
                         grid%step(2)+ &
                         real(icelly,dp)*pbc%box_size(2) 
                   u(3)=(grid%shift(3) + grid%kz(i))* &
                         grid%step(3)+ &
                         real(icellz,dp)*pbc%box_size(3)
                   call matvec3('N',tmpmat,u,x)
!                   x=matmul(tmpmat,u)
                   rx=sum(x(:)**2)

                   rho(i)=rho(i)+exp(-0.3*rx)
                 enddo
              enddo
          enddo
       endif    


      end do

      rho=100000*elec_st%rho(:,1)
      
      rhoav=sum(rho(:))/parallel%mydim

      rho=rho-rhoav

      pot_orig=pot_func
      rho=-rho ! doing the minus for rho so it will be a solution.

      pot_func=0 ! zeroing before call to hpotcg to make sure that
                 ! it is the calculated value that is being 
                 ! returned and not something that is left from
                 ! initialization.

      call hpotcg(parallel, rho, pot_func, grid%norder, grid%coe2, &
                        grid%lap_dir_num, vec)

      call lapmvs(parallel, pot_func, lapl, grid%coe2, grid%norder, &
                        grid%lap_dir_num, vec)
   

     

      call write_rho_1(parallel,pbc,grid,pot_func,rho,lapl)

      if(flag>0) then
       j=0
       write(9,*) 'hpotcg exception at:'
       do i=1,parallel%mydim
          if((abs(pot_orig(i))-abs(pot_func(i)))<0.0001) then
            j=j+1
          else
            u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
            u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
            u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,105) u, pot_orig(i), pot_func(i)
          endif
105    format(5(f15.9,1x))
       end do
       write(9,*) "number of diffs in pot: ", parallel%mydim-j
      endif


      j=0
      do i=1,parallel%mydim
         if(abs(lapl(i)-rho(i))<0.01) then
           j=j+1
         else
           u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
           u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
           u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,105) u, lapl(i), rho(i)
         endif
107   format(5(f15.9,1x))
      end do

      write(9,*) "number of diffs in rho: ",parallel%mydim-j

      end subroutine hpotcg_test

      subroutine hpotcg_slab_test(grid,pbc,parallel,elec_st)

      use constants
      use grid_module
      use pbc_module
      use parallel_data_module
      use electronic_struct_module

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel
      type (electronic_struct), intent(in) :: elec_st

!     work variables

      real(dp) :: pot_func(parallel%mydim), lapl(parallel%mydim), &
                  rho(parallel%mydim), pot_orig(parallel%mydim), &
                  brho(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)
      real(dp) :: u(3),x(3),gv(3),tmpvec(3),rx, dslab
      real(dp) :: min_pot, max_pot, pr

      integer i,j,flag,lpole
      integer icellx,icelly,icellz,nz1,nz2

      real(dp) :: tmpmat(3,3),tmpdet,tmptr,rhoav
      real(dp) :: ainv(3,3),ainvdot(3,3)

      write(9,*) "entered hpotcg_test function" 

      tmpmat=pbc%avec_norm

      ainv=tmpmat

      call mtrxin(ainv,tmpdet,tmptr)

      ainvdot=matmul(transpose(ainv),ainv)

!     Producing the test function data!

      pot_func=0
      rho=0    
      flag=1

      lpole=3
      nz1=-10
      nz2=10
      pr=1.0
      dslab=(nz2-nz1)*grid%step(3)
      min_pot=-pr*pi*dslab*grid%step(3)*4
      max_pot=pr*pi*dslab*grid%step(3)*4

!     First we produce a DIPOLE sheet test function. 

      do i=1,parallel%mydim

          u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
          u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
          u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))

          if(grid%kz(i)==nz1) then
              rho(i)=pr
          endif

          if(grid%kz(i)==nz2) then
              rho(i)=-pr
          endif

          if(grid%kz(i)<nz1) then
            pot_orig(i)=min_pot
          else
            if(grid%kz(i)>nz2) then
              pot_orig(i)=max_pot
            else
              pot_orig(i)= &
                min_pot+ &
                pr*4*2*pi*dslab*grid%step(3)*(grid%kz(i)-nz1)/(nz2-nz1)
            endif
          endif

      enddo

      call hartset_slab(grid,pbc,parallel,grid%norder,lpole, &
                        grid%coe2,rho,brho)

      brho=-brho ! doing the minus for rho so it will be a solution.

      call write_rho_1(parallel,pbc,grid,pot_func,rho,brho)

!      brho=-8*pi*rho ! this is just to check what happens if we bypass hartset

      pot_func=0 ! zeroing before call to hpotcg to make sure that
                 ! it is the calculated value that is being 
                 ! returned and not something that is left from
                 ! initialization.

      call hpotcg(parallel, brho, pot_func, grid%norder, grid%coe2, &
                        grid%lap_dir_num, vec)

      call lapmvs(parallel, pot_func, lapl, grid%coe2, grid%norder, &
                        grid%lap_dir_num, vec)
   

     

!      call write_rho_1(parallel,pbc,grid,pot_func,rho,brho)

      if(flag>0) then
       j=0
       write(9,*) 'hpotcg exception at:'
       do i=1,parallel%mydim
          if(abs((abs(pot_orig(i))-abs(pot_func(i))))<-0.0001) then
            j=j+1
          else
            u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
            u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
            u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,121) u, pot_orig(i), pot_func(i)
          endif
121    format(5(f15.9,1x))
       end do
       write(9,*) "number of diffs in pot: ", parallel%mydim-j
      endif


      j=0
      do i=1,parallel%mydim
         if(abs(lapl(i)-rho(i))<0.01) then
           j=j+1
         else
           u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
           u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
           u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,122) u, lapl(i), rho(i)
         endif
122   format(5(f15.9,1x))
      end do

      write(9,*) "number of diffs in rho: ",parallel%mydim-j

      end subroutine hpotcg_slab_test



      subroutine hpotcg_slab_test_2(grid,pbc,parallel,elec_st)

      use constants
      use grid_module
      use pbc_module
      use parallel_data_module
      use electronic_struct_module

      implicit none
!
!     Input/Output variables:
!
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
      type (parallel_data), intent(in) :: parallel
      type (electronic_struct), intent(in) :: elec_st

!     work variables

      real(dp) :: pot_func(parallel%mydim), lapl(parallel%mydim), &
                  rho(parallel%mydim), pot_orig(parallel%mydim), &
                  brho(parallel%mydim)
      real(dp) :: vec(parallel%nwedge+1)
      real(dp) :: u(3),x(3),gv(3),tmpvec(3),rx, dslab
      real(dp) :: min_pot, max_pot, pr, freq, freq_f, b_size
      real(dp) :: zdiff1, zdiff2

      integer i,j,flag,lpole
      integer icellx,icelly,icellz,nz1,nz2

      real(dp) :: tmpmat(3,3),tmpdet,tmptr,rhoav
      real(dp) :: ainv(3,3),ainvdot(3,3)

      write(9,*) "entered hpotcg_test function" 

      tmpmat=pbc%avec_norm

      ainv=tmpmat

      call mtrxin(ainv,tmpdet,tmptr)

      ainvdot=matmul(transpose(ainv),ainv)

!     Producing the test function data!

      pot_func=0
      rho=0    
      flag=1
      freq=2
      freq_f=freq*2*pi/pbc%box_size(2)

      lpole=5
      nz1=-40
      nz2=40
      pr=2.0
      dslab=(nz2-nz1)*grid%step(3)
      b_size=1/grid%step(3)
      min_pot=-pr*pi*dslab*grid%step(3)*4
      max_pot=pr*pi*dslab*grid%step(3)*4

!     First we produce a DIPOLE sheet test function. 

      do i=1,parallel%mydim

          u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
          u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
          u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))

          if(grid%kz(i)==nz1) then
              rho(i)=pr*cos(freq_f*u(2))
          endif

          if(grid%kz(i)==nz2) then
              rho(i)=-pr*cos(freq_f*u(2))
          endif


          zdiff1=abs(grid%kz(i)-nz1)*grid%step(3)
          zdiff2=abs(grid%kz(i)-nz2)*grid%step(3)

      

!          pot_orig(i)=(4*pi*pr*sin(freq_f*u(1))/b_size)*
          pot_orig(i)=((4*pi*pr)/b_size/freq_f)* &
                      (-exp(-freq_f*zdiff1)+exp(-freq_f*zdiff2))* &
                      cos(freq_f*u(2))

      enddo

      call hartset_slab(grid,pbc,parallel,grid%norder,lpole, &
                        grid%coe2,rho,brho)

      brho=-brho ! doing the minus for rho so it will be a solution.

      call write_rho_1(parallel,pbc,grid,pot_func,rho,brho)

!      brho=-8*pi*rho ! this is just to check what happens if we bypass hartset

      pot_func=0 ! zeroing before call to hpotcg to make sure that
                 ! it is the calculated value that is being 
                 ! returned and not something that is left from
                 ! initialization.

      call hpotcg(parallel, brho, pot_func, grid%norder, grid%coe2, &
                        grid%lap_dir_num, vec)

      call lapmvs(parallel, pot_func, lapl, grid%coe2, grid%norder, &
                        grid%lap_dir_num, vec)
   

     

!      call write_rho_1(parallel,pbc,grid,pot_func,rho,brho)

      if(flag>0) then
       j=0
       write(9,*) 'hpotcg exception at:'
       do i=1,parallel%mydim
          if(abs((abs(pot_orig(i))-abs(pot_func(i))))<-0.0001) then
            j=j+1
          else
            u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
            u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
            u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,131) u, pot_orig(i), pot_func(i), &
                        pot_orig(i)-pot_func(i)
          endif
131    format(3(f5.2,1x),3(f15.9,1x))
       end do
       write(9,*) "number of diffs in pot: ", parallel%mydim-j
      endif


      j=0
      do i=1,parallel%mydim
         if(abs(lapl(i)-rho(i))<0.01) then
           j=j+1
         else
           u(1)=((grid%shift(1)+grid%kx(i))*grid%step(1))
           u(2)=((grid%shift(2)+grid%ky(i))*grid%step(2))
           u(3)=((grid%shift(3)+grid%kz(i))*grid%step(3))
           write(9,132) u, lapl(i), rho(i)
         endif
132   format(5(f15.9,1x))
      end do

      write(9,*) "number of diffs in rho: ",parallel%mydim-j

      end subroutine hpotcg_slab_test_2

      subroutine write_rho(elec_st,parallel,pbc,grid,pot)

      use constants
      use cluster_module
      use electronic_struct_module
      use grid_module
      use potential_module
      use pseudo_potential_module
      use parallel_data_module
      use pbc_module
      implicit none

      type (electronic_struct), intent(in) :: elec_st
      type (parallel_data), intent(in) :: parallel
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
!     potential related data
      type (potential), intent(in) :: pot



      integer i,j,ioffset
      real(dp) :: u(3),x(3)

      real(dp) :: tmpmat(3,3)

      character(len=4) :: idstring

      write(idstring,'(I4.4)') parallel%iam

      tmpmat=pbc%avec_norm
      ioffset = parallel%irows(parallel%group_iam) - 1


      open(unit=13,file='chargedata.dat.'//idstring,form='formatted')

      do i = parallel%irows(parallel%group_iam) &
          ,parallel%irows(parallel%group_iam+1)-1
        
         j=i-ioffset
 
         u(1)=(grid%shift(1)+grid%kx(i))*grid%step(1)
         u(2)=(grid%shift(2)+grid%ky(i))*grid%step(2)
         u(3)=(grid%shift(3)+grid%kz(i))*grid%step(3)
         call matvec3('N',tmpmat,u,x)
!         x=matmul(tmpmat,u)

         write(13,706) x(1),x(2),x(3),elec_st%rho(j,1), &
                       pot%vion(j),pot%vhart(j),pot%vxc(j,1)

      enddo

706   format(3(f7.2,','),e13.6,',',e13.6,',',e13.6,',',e13.6)
      close(13)


      end subroutine write_rho

      subroutine write_rho_1(parallel,pbc,grid,vhart,rho,vhartcopy)

      use constants
      use cluster_module
      use grid_module
      use potential_module
      use pseudo_potential_module
      use parallel_data_module
      use pbc_module
      implicit none

      type (parallel_data), intent(in) :: parallel
      type (grid_data), intent(in) :: grid
      type (pbc_data), intent(in) :: pbc
!     potential related data
      real(dp), intent(in) :: vhart(parallel%mydim)
      real(dp), intent(in) :: rho(parallel%mydim)
      real(dp), intent(in) :: vhartcopy(parallel%mydim)


      integer i
      real(dp) :: u(3),x(3)
      real(dp) :: tmpmat(3,3)

      write(9,*) "entered write_rho_1 function"

      tmpmat=pbc%avec_norm

      open(unit=13,file='cgdatatest.dat',form='formatted')

      do i=1,parallel%mydim

         u(1)=(grid%shift(1)+grid%kx(i))*grid%step(1)
         u(2)=(grid%shift(2)+grid%ky(i))*grid%step(2)
         u(3)=(grid%shift(3)+grid%kz(i))*grid%step(3)
         call matvec3('N',tmpmat,u,x)
!         x=matmul(tmpmat,u)

         write(13,707) x(1),x(2),x(3),rho(i),vhart(i),vhartcopy(i)

      enddo

707   format(3(f15.9,','),e13.6,',',e13.6,',',e13.6)
      close(13)

      write(9,*) "end of write_rho_1 function"

      end subroutine write_rho_1

