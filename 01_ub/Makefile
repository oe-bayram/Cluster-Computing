all:
	mpicc -o main main.c
	mpirun -nodes tusci08,tusci12,tusci04,tusci00 main matrix1 matrix2 result_matrix
