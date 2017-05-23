all:
	git add -A
	git commit -m "added changes"
	git push origin master

run:
	mpicc -o main main.c -I/opt/DIS/include/
	mpirun -nodes tusci00,tusci04,tusci08,tusci12 main