# This PERL script compares to parsec output files
# usage : perl test_diff_4.pl file1 file2
# to make it usable - change the $PARSECDIR to the actual
# parsec installation directory.
#

use lib '$PARSECROOT/tools/perl_lib';
use parsec::pout;

$en_tol=0.001;

$force_tol=0.01;

$long_flag=1;

$i=$#ARGV+1;

my %test_data=();


if($#ARGV<1) { die "ERROR: must supply 2 parsec outfile names, got only $i params\n"; }


open(FILE1,$ARGV[0]) or die "could not open file $ARGV[0]\n";
close(FILE1);

open(FILE1,$ARGV[1]) or die "could not open file $ARGV[1]\n";
close(FILE2);

$i=0;
  
  $outp1 = new parsec::pout($ARGV[0]);
  $outp2 = new parsec::pout($ARGV[1]);

#  print "\n\n";

  if(($dd=abs($outp1->energy-$outp2->energy))<$en_tol) {
     #print "total energy : PASSED\n\n"; 
     $test_data{"tot_energy_test"}=0; }
  else {
     #print "total energy : FAILED\n\n";
     $test_data{"tot_energy_test"}=1;
     $test_diff{"tot_energy_test"}=$dd; }

  @eigen_psec1=$outp1->eigenvalues;
  @eigen_psec2=$outp2->eigenvalues;
  for $i  (0 .. $#eigen_psec1) {
   $test_data{"eigen_val[$i]"}=0;
   for $j  (0 .. 6) {
    if(abs($eigen_psec1[$i][$j]-$eigen_psec2[$i][$j])>$en_tol)
    { $test_data{"eigen_val[$i]"}=1; 
      print "eigen val diff: $eigen_psec1[$i][$j] $eigen_psec2[$i][$j]\n"; }
   }
  }



  %en_vals1=$outp1->energy_vals;
  %en_vals2=$outp2->energy_vals;

  foreach $key (keys %en_vals1)
  {
    $val_diff=$en_vals1{$key}-$en_vals2{$key};
    if(abs($val_diff)<$en_tol) {
    #  print "$key\t\t  : PASSED\n"; 
    $test_data{$key}=0; }
    else {
    #  print "$key\t\t  : FAILED\n";
    #  if($long_flag) { print "diff is: $val_diff\n";i
    $test_data{$key}=1;
    $test_diff{$key}=$val_diff;  }
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
 
  if($flen1-$flen2!=0) 
  { #print "forces  : FAILED\n"; 
    $test_data{"force_num"}=1;
    $test_diff{"force_num"}=$flen1-$flen2; }
  
  while($j<$#forces1_psec+1)
  {
  if(($val_diff=abs($forces1_psec[$j]-$forces2_psec[$j]))<$force_tol) {
     #print "forces[$j]    : PASSED\n";
     $test_data{"forces[$j]"}=0; }
  else
  {  #print "forces[$j]    : FAILED\n";
     $test_data{"forces[$j]"}=1;
     $test_diff{"forces[$j]"}=$val_diff; }

  $j=$j+1;
  }

  $test_total=0;

  foreach $key (keys %test_data) {
  $k=$k+$test_data{$key}; }

  if($k>0) 
  { print "FAILED\n"; }
  else
  { print "PASSED\n"; }
