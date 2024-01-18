

all:
	mpiexec -n 1 ../../src/parsec-ccm_intel-ifort-19.0.5.281.mpi
	tail parsec.out


.PHONY: clean

clean:
	rm -f out.* parsec.out parsec.dat fourier_V.dat
