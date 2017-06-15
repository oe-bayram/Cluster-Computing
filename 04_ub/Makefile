all:
	mpicc -o nbodydsm nbodydsm.c -I/opt/DIS/include/
run:
	mpirun -nodes tusci00,tusci08 nbodydsm 100 points out
