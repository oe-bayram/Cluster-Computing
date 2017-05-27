all:
	mpicc -o main main.c -I/opt/DIS/include/
	mpirun -nodes tusci00,tusci04,tusci08,tusci12 main matrix1 matrix2 matrix3
