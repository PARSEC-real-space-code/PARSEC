

all:
	mpiexec -n 8 $(PARSEC_EXEC)
	tail parsec.out


.PHONY: clean

clean:
	rm -f out.* parsec.out parsec.dat eigen.dat fourier_O.dat fourier_Zn.dat
