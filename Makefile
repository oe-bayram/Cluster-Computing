all:
	mpicc -o main main.c -I/opt/DIS/include/
run:
	mpirun -nodes tusci00,tusci04,tusci08,tusci12 main