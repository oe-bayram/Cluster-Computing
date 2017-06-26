all:
	mpicc -o nbodydsm nbodydsm.c -I/opt/DIS/include/
run:
	mpirun -np 8 -nodes tusci00,tusci08 nbodydsm 100 points10k out
