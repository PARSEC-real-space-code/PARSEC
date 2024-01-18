#!/usr/bin/perl -w
#
# This script takes the information from the parsec.out file and the output from 
# the plotdx software to generate an file called struct.xsf which can be 
# viewed with the XCrySDen visualization package
#
# The file necessary for atomic visualization is just the parsec.out file and 
# this script.  Running ./xcdstruc.pl produces a file struct.xsf which 
# can be opened by XCrySDen to view both the atomic positions (in Angstrom)
# and the forces on the atoms in Hartree/Angstrom
#
# The files necessary for atomic visualization and wavefunction visualization 
# are parsec.out, parsec.dat, plotdx (executable), and this script.  Run the 
# plotdx analysis and specify that the data be written to the file plotdx.dat.  
# Running this script reads the parsec.out, plotdx.dat, and LATTICE_VEC.dx and 
# produces the file struct.xsf which contains the data for xcrysden to plot 
# the atomic positions, forces, and wavefunction grid.
# 
# Ackowledgements:
#
# The part of this script which reads the atomic positions and forces is based 
# upon the script ballnstick.pl written by Murilo L. Tiago August 2005.  The 
# analysis software plotdx was written by Murilo L. Tiago August 2005.
#
# This script is written by Scott P. Beckman August 2006.  Copyright (C) for 
# the contents of this script is held by Scott P. Beckman.  Send comments and 
# bug reports to sbeckman@ices.utexas.edu . 
#
# This work is distributed under the Gnu Public License.
# 

##############
# Front Matter
##############

use POSIX qw(floor);

# List of chemical elements
@elist=("H", "He", 
"Li", "Be", "B", "C", "N", "O", "F", "Ne",
"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
"Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg");

# Scaling factors -- all unite in XCrySDen are Angstrom and Hartree
$scalecart=(1/0.529);
$scaleenergy=0.5;

# Source and target files
$source_file = "parsec.out";
$target_file = "struct.xsf";
$plotdx_file = "plotdx.dat";
$latdx_file = "LATTICE_VEC.dx";

# Open files for read and write
open(SOURCE_FILE,"<$source_file") || die "Can't open $source_file ";
open(TARGET_FILE,">$target_file") || die "Can't open $target_file ";
$nogridplot = 0;
open(PLOTDX,"<$plotdx_file") || ($nogridplot = 1);
#open(PLOTDX,"<$plotdx_file") ;
if($nogridplot == 0)
{
  open(LATDX,"<$latdx_file") || die "Can't open $latdx_file ";
}

###############################
# Reading data from parsec.out
###############################

# Read data from parsec.out
$jj = 0;
while (<SOURCE_FILE>) 
{
# Check if confined system (zero boundaries)
    if (/Confined system with zero boundary condition/)
    {
       $pbc=0;
       $_ = <SOURCE_FILE>;
       @line = split(' ',$_);
#       $radius=$line[4];
    }
# Check if periodic bounaries apply
    if (/Periodic boundary conditions/)
    {
       $pbc=1;
       <SOURCE_FILE>;
       $_ = <SOURCE_FILE>;
       @line = split(' ',$_);
       @a1=($line[0],$line[1],$line[2]);
       $_ = <SOURCE_FILE>;
       @line = split(' ',$_);
       @a2=($line[0],$line[1],$line[2]);
       $_ = <SOURCE_FILE>;
       @line = split(' ',$_);
       @a3=($line[0],$line[1],$line[2]);
    }
# Reading number of atoms, number of each type and checking if type is recognized by XCrySDen
    if (/Total number of atoms =\s*\d+/) 
    {
        @line = split(' ',$_);
        $natom = $line[5];
    }
    if (/There are\s*\d+ \w+ \s*atoms/) 
    {
        $jj = $jj + 1;
        @line = split(' ',$_);
        push(@atomnum,$line[2]);
        push(@spec,$line[3]);
        $warnflag=1;
        foreach(@elist)
        {
           if ($_ eq $line[3])
           {
              $warnflag=0;
           }
        }
        if ($warnflag==1)
        {
           print "WARNING element $line[3] is not recognized by XCrySDen!!!!!\n";
        }
    }
# Reading positions and forces (note that only the last set will be recorded)
    if (/coordinates and total forces \(after setting net force to zero\)/) 
    {
        @xx = () ;
        @yy = () ;
        @zz = () ;
        @fx = () ;
        @fy = () ;
        @fz = () ;
        $_ = <SOURCE_FILE>;
        $_ = <SOURCE_FILE>;
        $_ = <SOURCE_FILE>;
        foreach $ii (0..$natom-1) 
        {
            $_ = <SOURCE_FILE>;
            $_ =~/\s*\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/;
            ($xx[$ii],$yy[$ii],$zz[$ii],$fx[$ii],$fy[$ii],$fz[$ii]) = ($1,$2,$3,$4,$5,$6);
        }
    }
}

###############################
## Reading LATTICE_VEC.dx
################################

while (<LATDX>)
{
   if (/^object  1 class/)
   {
      $_ = <LATDX>;
      @line = split(' ',$_);
      @a1=($line[0],$line[1],$line[2]);
      $_ = <LATDX>;
      @line = split(' ',$_);
      @a2=($line[0],$line[1],$line[2]);
      $_ = <LATDX>;
      @line = split(' ',$_);
      @a3=($line[0],$line[1],$line[2]); 
   }
}		    

###################################
# Write Atomic Positions and Forces
###################################

# WARNING! PLOTTING PERIODIC STRUCTURE WITH GRID IS STILL BETA.

# If boundaries are periodic then print data in "CRYSTAL" format (careful to scale units correctly)
if($pbc == 1)
{
  if($nogridplot == 0)
  {
    print "WARNING! PLOTTING PERIODIC STRUCTURE WITH GRID IS STILL BETA IN XCRYSDEN v. 1.4.1!\n";
  }
	
   print TARGET_FILE "CRYSTAL\n";
   print TARGET_FILE "PRIMVEC\n";
   print TARGET_FILE $a1[0]/$scalecart,"   ",$a1[1]/$scalecart,"   ",$a1[2]/$scalecart,"\n";
   print TARGET_FILE $a2[0]/$scalecart,"   ",$a2[1]/$scalecart,"   ",$a2[2]/$scalecart,"\n";
   print TARGET_FILE $a3[0]/$scalecart,"   ",$a3[1]/$scalecart,"   ",$a3[2]/$scalecart,"\n";
   print TARGET_FILE "CONVVEC\n";
   print TARGET_FILE $a1[0]/$scalecart,"   ",$a1[1]/$scalecart,"   ",$a1[2]/$scalecart,"\n";
   print TARGET_FILE $a2[0]/$scalecart,"   ",$a2[1]/$scalecart,"   ",$a2[2]/$scalecart,"\n";
   print TARGET_FILE $a3[0]/$scalecart,"   ",$a3[1]/$scalecart,"   ",$a3[2]/$scalecart,"\n";
   print TARGET_FILE "PRIMCOORD\n";
   print TARGET_FILE $natom,"  1\n";
   $cumulative=0;
   $cumulativeold=0;
   for($ii=0;$ii<scalar(@atomnum);++$ii)
   {
      if($ii==0)
      {
         $cumulativeold=$cumulative;
         $cumulative=$cumulative+$atomnum[$ii];
         for($kk=0;$kk<$cumulative;++$kk)
         {
            print TARGET_FILE $spec[$ii],"   ",$xx[$kk]/$scalecart,"   ",$yy[$kk]/$scalecart,"   ",$zz[$kk]/$scalecart;
            print TARGET_FILE "   ",$fx[$kk]*$scalecart*$scaleenergy,"   ",$fy[$kk]*$scalecart*$scaleenergy,"   ",$fz[$kk]*$scalecart*$scaleenergy,"\n";
         }
      }
      else
      {
         $cumulativeold=$cumulative;
         $cumulative=$cumulative+$atomnum[$ii];
         for($kk=$cumulativeold;$kk<$cumulative;++$kk)
         {
            print TARGET_FILE $spec[$ii],"   ",$xx[$kk]/$scalecart,"   ",$yy[$kk]/$scalecart,"   ",$zz[$kk]/$scalecart;
            print TARGET_FILE "   ",$fx[$kk]*$scalecart*$scaleenergy,"   ",$fy[$kk]*$scalecart*$scaleenergy,"   ",$fz[$kk]*$scalecart*$scaleenergy,"\n";
         }
      }
   }
}
# If boundaries are aperiodic then print data in "MOLECULE" format (careful to scale units correctly)
else
{
   print TARGET_FILE "MOLECULE\n";
   print TARGET_FILE "ATOMS\n";
   $cumulative=0;
   $cumulativeold=0;
   for($ii=0;$ii<scalar(@atomnum);++$ii)
   {
      if($ii==0)
      {
         $cumulativeold=$cumulative;
         $cumulative=$cumulative+$atomnum[$ii];
         for($kk=0;$kk<$cumulative;++$kk)
         {
            print TARGET_FILE $spec[$ii],"   ",$xx[$kk]/$scalecart,"   ",$yy[$kk]/$scalecart,"   ",$zz[$kk]/$scalecart;
            print TARGET_FILE "   ",$fx[$kk]*$scalecart*$scaleenergy,"   ",$fy[$kk]*$scalecart*$scaleenergy,"   ",$fz[$kk]*$scalecart*$scaleenergy,"\n";
         }
      } 
      else
      {
         $cumulativeold=$cumulative;
         $cumulative=$cumulative+$atomnum[$ii];
         for($kk=$cumulativeold;$kk<$cumulative;++$kk)
         {
            print TARGET_FILE $spec[$ii],"   ",$xx[$kk]/$scalecart,"   ",$yy[$kk]/$scalecart,"   ",$zz[$kk]/$scalecart;
            print TARGET_FILE "   ",$fx[$kk]*$scalecart*$scaleenergy,"   ",$fy[$kk]*$scalecart*$scaleenergy,"   ",$fz[$kk]*$scalecart*$scaleenergy,"\n";
         }
      }
   }
}

###################################################
# Read data from plotdx.dat and write to struct.xsf
###################################################

if($nogridplot == 0)
{

  print TARGET_FILE "\nBEGIN_BLOCK_DATAGRID_3D\n";
  print TARGET_FILE "   written-by-xcdstruc.pl\n";
  print TARGET_FILE "   BEGIN_DATAGRID_3D_this_is_3Dgrid#1\n";

  while (<PLOTDX>) 
  {
    if (/^object 2 class gridconnections counts/)
    {
       @line = split(' ',$_);
       @n=($line[5],$line[6],$line[7]);
    }
  }

  @origin=((-$a1[0]-$a2[0]-$a3[0])/2,(-$a1[1]-$a2[1]-$a3[1])/2,(-$a1[2]-$a2[2]-$a3[2])/2);
  $ntotal=$n[0]*$n[1]*$n[2];

  if (($ntotal % 5) != 0)
  {
     $nlines=floor($ntotal/5)+1;
  }
  else
  {
     $nlines=floor($ntotal/5);
  }

  print TARGET_FILE "      ",$n[0],"  ",$n[1],"  ",$n[2],"\n";
  print TARGET_FILE "      ",$origin[0]/$scalecart,"  ",$origin[1]/$scalecart,"  ",$origin[2]/$scalecart,"\n";
  print TARGET_FILE "      ",$a1[0]/$scalecart,"  ",$a1[1]/$scalecart,"  ",$a1[2]/$scalecart,"\n";
  print TARGET_FILE "      ",$a2[0]/$scalecart,"  ",$a2[1]/$scalecart,"  ",$a2[2]/$scalecart,"\n";
  print TARGET_FILE "      ",$a3[0]/$scalecart,"  ",$a3[1]/$scalecart,"  ",$a3[2]/$scalecart,"\n";


# I know that this is ugly, but there is no rewind command in perl 
# so I'm stuck reopening the data file...
  open(PLOTDX,"<$plotdx_file"); 

  while (<PLOTDX>)
  {
    if (/^object 1 class array/)
    {
       for($i=0;$i<$nlines;++$i)
       {
	  $line=<PLOTDX>;
          print TARGET_FILE $line;
       }
    }
  }

  print TARGET_FILE "   END_DATAGRID_3D\n";                      
  print TARGET_FILE "END_BLOCK_DATAGRID_3D\n";
}


close(TARGET_FILE);



