#!/usr/bin/perl

#
# Script to generate spherical cluster of elemental clusters
#
# Written by S.P.Beckman September 2005
# Modified by M.L. Tiago, May 2006
#
# Copyright (C) 2005 S.P.Beckman
#    This is free software; you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#    You can receive a copy of the GNU General Public License
#    by writing to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#

use POSIX;

##############
# define parameters
##############

# These are the unit lattice vectors. We also need to specify the
# shortest size of the cell (for the cell replicas).
#
#           CUBIC CRYSTALS
#
# lattice parameter in Angstrom
#$amin=6.652;

#@a1_vec=($amin,0,0);
#@a2_vec=(0,$amin,0);
#@a3_vec=(0,0,$amin);
#
##   end of CUBIC CRYSTALS

#           HEXAGONAL CRYSTALS
#
# parameters "a" and "c" for HCP structure, in Angstrom
$amin=4.7432;
$c_latt=sqrt(8/3)*$amin;

@a1_vec=($amin/2,-$amin*sqrt(3)/2,0);
@a2_vec=($amin/2, $amin*sqrt(3)/2,0);
@a3_vec=(0,0,$c_latt);
#
##   end of HEXAGONAL CRYSTALS

# Here, specify the center of the cluster, in lattice vector units.
# It is not necessarily on an atom site.

# bridge center, cubic crystals
#$shift1=1/4;
#$shift2=1/4;
#$shift3=1/4;

# cage center, cubic crystals
#$shift1=-1/2;
#$shift2=0.0;
#$shift3=0.0;

# cage center, hexagonal crystals
#$shift1=1/3;
#$shift2=2/3;
#$shift3=1/4;

# atom center
$shift1=0.0;
$shift2=0.0;
$shift3=0.0;

# This is the list of atoms in the cell, in lattice vector units.

# diamond crystal
#@initcoord1=(0+$shift1,1/2+$shift1,1/2+$shift1,0+$shift1,1/4+$shift1,3/4+$shift1,3/4+$shift1,1/4+$shift1);
#@initcoord2=(0+$shift2,1/2+$shift2,0+$shift2,1/2+$shift2,1/4+$shift2,3/4+$shift2,1/4+$shift2,3/4+$shift2);
#@initcoord3=(0+$shift3,0+$shift3,1/2+$shift3,1/2+$shift3,1/4+$shift3,1/4+$shift3,3/4+$shift3,3/4+$shift3);

# BCC crystal
#@initcoord1=(0+$shift1,1/2+$shift1);
#@initcoord2=(0+$shift2,1/2+$shift2);
#@initcoord3=(0+$shift3,1/2+$shift3);

# FCC crystal
#@initcoord1=(0+$shift1,1/2+$shift1,1/2+$shift1,0+$shift1);
#@initcoord2=(0+$shift2,1/2+$shift2,0+$shift2,1/2+$shift2);
#@initcoord3=(0+$shift3,0+$shift3,1/2+$shift3,1/2+$shift3);

# HCP crystal
@initcoord1=($shift1, 1/3+$shift1);
@initcoord2=($shift2,-1/3+$shift2);
@initcoord3=($shift3, 1/2+$shift3);

# radius of sphere in Angstrom
$radius=15;

# Atomic species:
$species=("Fe");

###############################################################################
# end of input parameters; no need to change anything below this point.
###############################################################################

##############
# create atomic block
##############

# because the actual radius will slide some, use 
# twice as large of cell (8x number of atoms)

$radeff=2*$radius;

#figure out how many cells need shifting 

$numcells=ceil($radeff/$amin);
$numcells=$numcells+1 if ($numcells%2 !=0);

# create cell in I quadrant  
# define init* variables as those in a simple cube

# span the cube through-out the rest of quadrant I

for($i=0;$i<$numcells;++$i)
{
    for($j=0;$j<$numcells;++$j)
    {
	for($k=0;$k<$numcells;++$k)
	{
	    foreach(@initcoord1)
	    {
		push(@spec,$species);
	    }
	    foreach(@initcoord1)
	    {
		push(@coord1,$_+$i);
	    }
	    foreach(@initcoord2)
	    {
		push(@coord2,$_+$j);
	    }
	    foreach(@initcoord3)
	    {
		push(@coord3,$_+$k);
	    }
	}
    }    
}

# span the remaining 7 quadrants with those in I
# numquadi is number of atoms in quad I
# tau the translation into the other 7 quads.
# totnum is total number of atoms in all 8 quad

$numquadi=scalar(@spec);
$totnum=8*$numquadi;
$tau=$numcells;

for($i=0;$i<$numquadi;++$i)
{
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]-$tau);
    push(@coord2,$coord2[$i]);
    push(@coord3,$coord3[$i]);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]);
    push(@coord2,$coord2[$i]-$tau);
    push(@coord3,$coord3[$i]);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]);
    push(@coord2,$coord2[$i]);
    push(@coord3,$coord3[$i]-$tau);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]-$tau);
    push(@coord2,$coord2[$i]-$tau);
    push(@coord3,$coord3[$i]);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]);
    push(@coord2,$coord2[$i]-$tau);
    push(@coord3,$coord3[$i]-$tau);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]-$tau);
    push(@coord2,$coord2[$i]);
    push(@coord3,$coord3[$i]-$tau);
	 
    push(@spec,$spec[$i]);
    push(@coord1,$coord1[$i]-$tau);
    push(@coord2,$coord2[$i]-$tau);
    push(@coord3,$coord3[$i]-$tau);
	 
}


# identify the atoms inside and outside the sphere
# change coordinates from lattice vector to cartesian
    
@inspec=();
@inxcord=();
@inycord=();
@inzcord=();
$atomcount=0;
    
for($i=0;$i<$totnum;++$i)
{
	$tmp_spec=@spec[$i];
	$tmp_xcord=$coord1[$i]*$a1_vec[0]+$coord2[$i]*$a2_vec[0]+$coord3[$i]*$a3_vec[0];
	$tmp_ycord=$coord1[$i]*$a1_vec[1]+$coord2[$i]*$a2_vec[1]+$coord3[$i]*$a3_vec[1];
	$tmp_zcord=$coord1[$i]*$a1_vec[2]+$coord2[$i]*$a2_vec[2]+$coord3[$i]*$a3_vec[2];
	
	$mag=sqrt($tmp_xcord*$tmp_xcord+$tmp_ycord*$tmp_ycord+$tmp_zcord*$tmp_zcord);
	
	if($mag<$radius)
	{
	    push(@inspec,$tmp_spec);
	    push(@inxcord,$tmp_xcord);
	    push(@inycord,$tmp_ycord);
	    push(@inzcord,$tmp_zcord);
	    ++$atomcount;	    
	}
}
    
print " total = ",$atomcount,"\n";

#print results to file in xbsa and xyz formats	

print "OPENED FILE and writing!\n";
open OUTFILE,">clusterout.bs";
open OUTFILE_xyz,">clusterout.xyz";

print OUTFILE_xyz $atomcount,"\n";
print OUTFILE_xyz "\n";
	
while(@inspec)
{
	$xdis = shift(@inxcord);
	$ydis = shift(@inycord);
	$zdis = shift(@inzcord);
	$mag=sqrt($xdis*$xdis+$ydis*$ydis+$zdis*$zdis);
	print OUTFILE "atom  ",shift(@inspec),"\t",$xdis,"\t",$ydis,"\t",$zdis,"\n";
	printf(OUTFILE_xyz "  %s %14.9f %14.9f %14.9f %14.9f  \n",$species,$xdis,$ydis,$zdis,$mag);
}
    
print OUTFILE "spec   ",$species,"   0.50   DarkOliveGreen\n";
print OUTFILE "bonds  ",$species," ",$species,"   1.4   2.8   0.2   0.50\n";
	
close(OUTFILE);

close(OUTFILE_xyz);
