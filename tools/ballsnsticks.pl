#!/usr/bin/perl -w
#
# Copyright (C) 2005 Finite Difference Research Group
# This file is part of parsec, http://www.ices.utexas.edu/parsec/
#
# Reads a parsec.out file and creates file BALLS.dx with atomic coordinates,
# to be used with DataExplorer. File parsec.out must have at least one table
# of final forces. Only coordinates from the last table are retained.
#
# Murilo L Tiago, UTexas-Austin, August 2005, mtiago@ices.utexas.edu
#
$source_file = "parsec.out";
$target_file = "BALLS.dx";

open(SOURCE_FILE,"<$source_file") || die "Can't open $source_file: ";
open(TARGET_FILE,">$target_file") || die "Can't open $target_file: ";

print TARGET_FILE "#\n";
print TARGET_FILE "#\n";
print TARGET_FILE "#   BALL AND STICK INFO:\n";

$jj = 0;
while (<SOURCE_FILE>) {
    if (/Total number of atoms =\s*\d+/) {
        @line = split(' ',$_);
        $natom = $line[5];
    }
    if (/There are\s*\d+ \w+ \s*atoms/) {
        $jj = $jj + 1;
        @line = split(' ',$_);
        push(@atomnum,$line[2]);
           foreach $ii (1..$atomnum[-1]) {
              push(@type,$jj);
        }
    }
    if (/coordinates and total forces \(after setting net force to zero\)/) {
        @xx = () ;
        @yy = () ;
        @zz = () ;
        $_ = <SOURCE_FILE>;
        $_ = <SOURCE_FILE>;
        $_ = <SOURCE_FILE>;
        foreach $ii (1..$natom) {
            $_ = <SOURCE_FILE>;
            $_ = /^\s*\d+\s+/;
            ($xx[$ii],$yy[$ii],$zz[$ii]) = split(' ',$');
        }
    }
}

print TARGET_FILE "object \"ballcoord\" array type float rank 1 shape 3 items      $natom data follows\n";
print TARGET_FILE "# the positions of the atoms\n";
foreach $ii (1..$natom) {
    print TARGET_FILE "   $xx[$ii]    $yy[$ii]    $zz[$ii] \n";
}
print TARGET_FILE "object \"data\" array type float rank 0 items     $natom data follows \n";
foreach $ii (1..$natom) {
    print TARGET_FILE "      ","$type[$ii-1]","\.0000  \n";
}
print TARGET_FILE " object \"molecule\" field \n";
print TARGET_FILE " component \"positions\" value \"ballcoord\" \n";
print TARGET_FILE " component \"data\" value \"data\" \n";
print TARGET_FILE "end \n";
