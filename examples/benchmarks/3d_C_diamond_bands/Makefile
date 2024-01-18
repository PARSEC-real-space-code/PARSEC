

all:
	mpiexec -n 8 $(PARSEC_EXEC)
	tail parsec.out


.PHONY: clean

clean:
	rm -f out.* parsec.out parsec.dat bands.dat eigen.dat bands_plot.dat myjob.e* myjob.o*
