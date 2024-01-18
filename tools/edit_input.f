c
c     Copyright (C) 2005 Finite Difference Research Group
c     This file is part of parsec, http://www.ices.utexas.edu/parsec/
c
c     Converts old-format parsec.in files into new format
c
c     author: Murilo Tiago, Univ. of Minnesota, Oct 2004
c
c     usage: ./edit_input odlfilename newfilename
c
c     compilation (IBM): xlf -o edit_input edit_input.f
c
      program edit_input

      implicit none
      integer :: ibound,ii,jj,nn,nel,nat,ntotal
      integer, dimension(100) :: ivec

      integer iargc,nargs
      external iargc
      external getarg

      double precision, dimension (100) :: avec
      character (len=255) :: infile, outfile,string
      character (len=2) :: name
      character (len=1) :: loc

      nargs = iargc()
      if (nargs .ne. 2) then
         write(0,*) 'Usage: edit_input oldfilename newfilename'
         stop 
      endif

      call getarg(1,infile)
      call getarg(2,outfile)

      open (10,file=infile,form='formatted',status='old')
      open (20,file=outfile,form='formatted')

      write(20,*) '#'
      write(20,*) '# This is an input file automatically generated from'
      write(20,*) '# ',trim(infile)
      write(20,*) '# The user is strongly encouraged to check sanity of'
      write(20,*) '# this file before using it in a parsec90 run!!!'
      write(20,*) '#'
      write(20,*)

      read(10,*) string
      write(20,*) '# ',trim(string)

      read(10,*) ii
      if (ii == 0) then
         write(20,*) 'Restart_run: .false. '
      else
         write(20,*) 'Restart_run: .true. '
      endif
      write(20,*)

      write(20,*) 'Old_Pseudopotential_Format: .true.'
      write(20,*) 'Binary_KB_Pseudopotential_File: .true.'
      write(20,*)

      read(10,*) ibound
      if (ibound == 0) then
         read(10,*) (avec(ii),ii=1,2)
         write(20,*) 'Periodic_System: .false.'
         write(20,*)
         write(20,*) 'Boundary_Sphere_Radius: ',avec(1)
         write(20,*) 'Grid_Spacing: ',avec(2)
      else
         read(10,*) (avec(ii),ii=1,4)
         write(20,*) 'Periodic_System: .true.'
         write(20,*)
         write(20,*) 'begin Cell_shape'
         write(20,*) (avec(ii),ii=1,3)
         write(20,*) 'end cell_shape'
         write(20,*)
         write(20,*) 'Grid_Spacing: ',avec(4)
      endif
      write(20,*)
      read(10,*) ii
      write(20,*) 'Expansion_Order: ',ii
      write(20,*)

      read(10,*) ivec(1),ivec(2)
      write(20,*) 'States_Num: ',ivec(1)
      write(20,*) 'Net_charges: ',ivec(2)
      read(10,*) avec(1),ivec(1)
      write(20,*) 'Fermi_Temp: ',avec(1)
      write(20,*) 'Max_Iter: ',ivec(1)
      read(10,*) avec(1),avec(2)
      write(20,*) 'Convergence_Criterion: ',avec(1)
      write(20,*) 'Diag_Tolerance: ',avec(2)
      read(10,*) string
      write(20,*) 'Eigensolver: ',trim(string)
      write(20,*)

      read(10,*) ii
      if (ii .eq. 0) then
         write(20,*) 'Mixing_Method: Anderson'
      else
         write(20,*) 'Mixing_Method: Broyden'
      endif
      read(10,*) avec(1),ivec(1)
      write(20,*) 'Mixing_Param: ',avec(1)
      if (ii == 1) write(20,*) 'Memory_Param: ',ivec(1)

      write(20,*)
      read(10,*) nel
      write(20,*) 'Atom_Types_Num: ',nel
      ntotal = 0
      do nn=1,nel
         write(20,*)
         write(20,*) '#------------ new atom type ---------------- '
         read(10,*) name,avec(1)
         write(20,*) 'Atom_Type: ',trim(name)
         write(20,*) 'Core_Cutoff_Radius: ',avec(1)
         read(10,*) ii,loc,(ivec(jj),jj=1,ii)
         write(20,*) 'Local_Component: ',loc
         write(20,*) 'Potential_Num: ',ii
         write(20,*)
         write(20,*) 'begin Electron_Per_Orbital'
         write(20,*) '# S P D F'
         write(20,*) (ivec(jj),jj=1,ii)
         write(20,*) 'end Electron_Per_Orbital'
         write(20,*)
         read(10,*) nat,jj
         ntotal = ntotal + nat
         write(20,*) 'Move_Flag: ',jj
         write(20,*) 'begin Atom_Coord'
         do jj=1,nat
            read(10,*) (avec(ii),ii=1,3)
            write(20,*) (avec(ii),ii=1,3)
         enddo
         write(20,*) 'end Atom_Coord'
         write(20,*)
         write(20,*) '#------------ end of atom type ------------- '
         write(20,*)
      enddo
      write(20,*) 'Total_Atom_Num: ',ntotal

      write(20,*)
      read(10,*) name
      write(20,*) 'Correlation_Type: ',trim(name)
      read(10,*) avec(1)
      write(20,*) 'Ion_Energy_Diff: ',avec(1)
      read(10,*) ii
      if (ii /= 1)  write(20,*) 'Skip_force: .true.'
      read(10,*) ii
      if (ii .eq. 0) then
         write(20,*) 'Minimization: none'
      elseif (ii .eq. 1) then
         write(20,*) 'Minimization: simple'
      elseif (ii .eq. 2) then
         write(20,*) 'Minimization: BFGS'
      elseif (ii .eq. 3) then
         write(20,*) 'Minimization: manual'
      endif

      read(10,*) ivec(1),avec(1),avec(2),avec(3)
      write(20,*) 'Movement_Num: ',ivec(1)
      write(20,*) 'Force_Min: ',avec(1)
      write(20,*) 'Max_Step: ',avec(2)
      write(20,*) 'Min_Step: ',avec(3)
      write(20,*)
      read(10,*) avec(1)
      read(10,*) ii,jj
      if (ii == 1) then
         write(20,*) 'Molecular_Dynamics: .true.'
         if (jj == 1) write(20,*) 'Restart_Mode: .true.'
      endif
      read(10,*) ii

      if (ii .eq. 0) then
         write(20,*) 'Cooling_Method: log'
      elseif (ii .eq. 1) then
         write(20,*) 'Cooling_Method: linear'
      elseif (ii .eq. 2) then
         write(20,*) 'Cooling_Method: stair'
      endif

      read(10,*) (avec(ii),ii=1,3)
      write(20,*) 'Tinit: ',avec(1)
      write(20,*) 'T_Final: ',avec(2)
      write(20,*) 'T_Step: ',avec(3)
      write(20,*)
      read(10,*) avec(1),avec(2),ivec(1)
      write(20,*) 'Time_Step: ',avec(1)
      write(20,*) 'Friction_Coefficient: ',avec(2)
      write(20,*) 'Step_num: ',ivec(1)
      write(20,*)
      read(10,*) ii
      if (ii == 1) write(20,*) 'Polarizability: .true.'
      write(20,*)
      read(10,*) avec(1)
      write(20,*) 'Electric_Field: ',avec(1)
      write(20,*)
      read(10,*) ii
      if (ii == 1) then
         write(20,*) 'Spin_Polarization: .true.'
         read(10,*) ii
         write(20,*) 'Occupied_States: ',ii
         read(10,*) avec(1),avec(2)
         write(20,*) 'Occup_Num_Up: ',avec(1)
         write(20,*) 'Occut_Num_Down: ',avec(2)
      endif
      read(10,*,iostat=ii) nn
      if (nn == 1 .and. ii == 0) write(20,*)
     >     'output_all_states .true.'
      write(20,*)

      end program edit_input
