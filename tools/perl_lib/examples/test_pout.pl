# This PERL script parses parsec out files and prints
# some of the parsed variables.
#
# usage : perl test_pout.pl file1
# to make it usable - change the $PARSECDIR to the actual
# parsec installation directory.
#
use lib '$PARSECROOT/tools/perl_lib';

use parsec::pout;

$i=$#ARGV+1;

if($#ARGV<0) { die "ERROR: must supply parsec outfile name, got only $i params\n"; }

print "got the name $ARGV[0]\n";

open(FILE1,$ARGV[0]) or die "could not open file $ARGV[0]\n";

print "opened successfuly $ARGV[0] \n";

close(FILE1);


$i=0;
  print $i
  $outp = new parsec::pout($ARGV[0]);
  print "Total Energy is: ", $outp->energy, " \n\n";
  @eigen_psec=$outp->eigenvalues;
  %energy_vals=$outp->energy_vals;
  print "$#eigen_psec\n";
  for $i  (0 .. $#eigen_psec) {
   for $j  (0 .. 6) {
    print "$eigen_psec[$i][$j], ";
   }
  print "\n";
  }
  print "\n";

  print "hartree energy is      :  $energy_vals{'hartree'}\n";
  print "eigenvalue energy is   :  $energy_vals{'eigenvalues_energy'}\n";
  print "vxc_rho is             :  $energy_vals{'vxc_rho'}\n";
  print "Exc is                 :  $energy_vals{'e_xc'}\n";
  print "alpha term is          :  $energy_vals{'alpha'}\n";
  print "ion ion energy is      :  $energy_vals{'ion_ion_energy'}\n";
  print "electron_ion energy is :  $energy_vals{'electron_ion_energy'}\n";
  print "total energy is        :  $energy_vals{'total_energy'}\n";

  print "\n\n";

  $j=0;
  @forces_psec=$outp->forcevectors;

  $flen=$#forces_psec+1;
 
  print "forces array length is: $flen\n";
  
  while($j<$#forces_psec+1)
  {
  print "forces[$j] $forces_psec[$j]\n";
  $fnorm=abs($forces_psec[$j]);
  print "norm[$j] $fnorm\n";
  $j=$j+1;
  }

