all:
	git add -A
	git commit -m "added changes"
	git push origin master

run:
	mpicc -o main main.c -I/opt/DIS/include/
	mpirun -nodes tusci12 main matrix1 matrix2 matrix3