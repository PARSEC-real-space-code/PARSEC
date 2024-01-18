#!/bin/sh
# Author: Charles Lena
# Description: Makes and deploys several compilations of PARSEC.
#              Not intended for use by non-developers. This is
#              handy for testing new additions of code in multiple
#              compilers. This is something that should be done before
#              any major (1.4)/ minor(1.4.3) release. Really it's just
#              good practice
#######################################################################
# Optional clean
make -f config/make.gfortran_debug -f Makefile_Advanced -j12 cleanall

# Compile serial versions with full error checks
make -f config/make.intel12_debug -f Makefile_Advanced -j12
make -f config/make.gfortran_debug -f Makefile_Advanced -j12
#make -f config/make.g95_debug -f Makefile_Advanced -j12

# Compile serial versions with no error checks (fast versions)
make -f config/make.intel12       -f Makefile_Advanced -j12
make -f config/make.gfortran      -f Makefile_Advanced -j12

# Compile openmp versions (coming soon)
echo No OpenMP versions yet

# Compile mpi versions
#make -f config/make.intel12_mpi -f Makefile_Advanced -j12
#make -f config/make.gfortran_mpi_debug -f Makefile_Advanced -j12

# Deploys in a fashion that the test bed can make use of
for file in `ls parsec-*.ser`; do
	nfile=`echo $file | sed -e "s/parsec-\(.*\)\.\([a-z]*\)/\1.\2/g"`
	echo $nfile
	echo ../testbed/parsecs/1.4_$nfile/parsec/src
	mkdir -p ../testbed/parsecs/1.4_$nfile/parsec/src
	mkdir -p ../testbed/parsecs/1.4_$nfile/parsec/benchmarks
	cp ./$file ../testbed/parsecs/1.4_$nfile/parsec/src/parsec.ser
done
for file in `ls parsec-*.mpi`; do
	nfile=`echo $file | sed -e "s/parsec-\(.*\)\.\([a-z]*\)/\1.\2/g"`
	echo $nfile
	echo ../testbed/parsecs/1.4_$nfile/parsec/src
	mkdir -p ../testbed/parsecs/1.4_$nfile/parsec/src
	mkdir -p ../testbed/parsecs/1.4_$nfile/parsec/benchmarks
	cp ./$file ../testbed/parsecs/1.4_$nfile/parsec/src/parsec.mpi
done

# Use scripts to setup benchmarks
cd ../testbed
python ../tools/setup_benchmarks.py benchmarksets/fast --parsecs parsecs
python ../tools/run_benchmarks.py --parsecs parsecs

# Run everything
python -c "print 'No Running Scripts Yet'"
cd ../src


