#!/bin/bash
$1 --version | head -n 1 | sed -n "s/GNU Fortran ([^)]*) \([0-9\.]*\)[ ]*\([0-9\.]*\)[ ]*(*[^)]*)*/\1/p" 
