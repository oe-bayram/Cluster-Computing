all:
	mpicc -o main main.c
	mpirun -nodes tusci08,tusci12,tusci04,tusci00 main matrix1 matrix2 result_matrix

big:
	mpicc -o main main.c -I/opt/DIS/include/
	mpirun -nodes tusci00,tusci04,tusci08,tusci12 main big_matrix1 big_matrix2 result_matrix