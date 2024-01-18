# This PERL script compares to parsec output files
# usage : perl test_diff_3.pl file1 file2
# to make it usable - change the $PARSECDIR to the actual
# parsec installation directory.
#
use lib '$PARSECDIR/tools/perl_lib';
use parsec::pout;


$en_tol=0.001;

$force_tol=0.01;

$long_flag=1;

$i=$#ARGV+1;


if($#ARGV<1) { die "ERROR: must supply 2 parsec outfile names, got only $i params\n"; }


open(FILE1,$ARGV[0]) or die "could not open file $ARGV[0]\n";
close(FILE1);

open(FILE1,$ARGV[1]) or die "could not open file $ARGV[1]\n";
close(FILE2);

$i=0;
  
  $outp1 = new parsec::pout($ARGV[0]);
  $outp2 = new parsec::pout($ARGV[1]);

  print "\n\n";

  if(abs($outp1->energy-$outp2->energy)<$en_tol) {
     print "total energy : PASSED\n\n"; }
  else {
     print "total energy : FAILED\n\n"; }

  %en_vals1=$outp1->energy_vals;
  %en_vals2=$outp2->energy_vals;

  foreach $key (keys %en_vals1)
  {
    $val_diff=$en_vals1{$key}-$en_vals2{$key};
    if(abs($val_diff)<$en_tol) {
      print "$key\t\t  : PASSED\n"; }
    else {
      print "$key\t\t  : FAILED\n";
      if($long_flag) { print "diff is: $val_diff\n"; } }
  }
#  print "hartree energy is      :  $energy_vals{'hartree'}\n";
#  print "eigenvalue energy is   :  $energy_vals{'eigenvalues_energy'}\n";
#  print "vxc_rho is             :  $energy_vals{'vxc_rho'}\n";
#  print "Exc is                 :  $energy_vals{'e_xc'}\n";
#  print "alpha term is          :  $energy_vals{'alpha'}\n";
#  print "ion ion energy is      :  $energy_vals{'ion_ion_energy'}\n";
#  print "electron_ion energy is :  $energy_vals{'electron_ion_energy'}\n";
#  print "total energy is        :  $energy_vals{'total_energy'}\n";


  $j=0;
  @forces1_psec=$outp1->forcevectors;
  @forces2_psec=$outp2->forcevectors;

  $flen1=$#forces1_psec+1;
  $flen2=$#forces2_psec+1;
 
  if($flen1-$flen2!=0) { print "forces  : FAILED\n"; }
  
  while($j<$#forces1_psec+1)
  {
  if(abs($forces1_psec[$j]-$forces2_psec[$j])<$force_tol) {
     print "forces[$j]    : PASSED\n";}
  else
  {  print "forces[$j]    : FAILED\n";}

  $j=$j+1;
  }

