

all:
	mpiexec -n 8 $(PARSEC_EXEC)
	tail parsec.out


.PHONY: clean

clean:
	rm -f out.* parsec.out parsec.dat eigen.dat atom.cor bfgs.dat relax_restart.dat iterate.dat
