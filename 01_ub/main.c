#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define MAX(a,b) (((a)>(b))?(a):(b))
#define INFO_BUFFER 4
#define MATRIX_INFO 0                   //Tag for MPI_Send
#define MATRIX_A    1                   //Tag for MPI_Send
#define MATRIX_B    2                   //Tag for MPI_Send
#define MATRIX_C    3                   //Tag for MPI_Send



//basic matrix struct, used to save the size (rows and columns) and values (*matrix)

struct matrix_s
{
    int *matrix;
    int rows;
    int columns;
} nmatrix = {NULL, 0, 0};

typedef struct matrix_s matrix;


/*
    read_matrix(filename, *A)
    reads a matrix from file <filename> and fills the content into the A matrix struct:
    1- Open File
    2- Get File size for the buffer
    3- Count rows and columns
    4- Populate matrix A
*/

int read_matrix(char *filename, matrix *A)
{
   printf("Read matrix...");
   MPI_File fh;
   MPI_Status status;
   MPI_Offset size;

   MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
   MPI_File_get_size(fh, &size);
 
   printf("file size: %d\n", size);

   char *buffer= (char *) malloc(size * sizeof(char));;
   A->rows = 0;
   A->columns = 0;

   MPI_File_seek(fh, 0, MPI_SEEK_SET);
   MPI_File_read(fh, buffer, size, MPI_CHAR, &status); 
       
   int i;
   for(i = 0; i < size; i++)
   {
       if(A->rows == 0 && buffer[i] == ' ')
       {
   A->columns++;
       }
       if(buffer[i] == '\n')
       {
           A->rows++;
       }
   }

   A->columns++;
   int matrix_size = A->rows * A->columns;
   
   A->matrix = (int *) malloc(matrix_size * sizeof(int));
   //printf(", max byte length: %d\n", max_byte_length);

   int index = 0;
   char *pos = buffer;
       
   while(pos < buffer + size)
   {
       char *end;
       long int value = strtol(pos, &end, 10);
       if( value == 0L && end == pos)
   break;

       A->matrix[index] = (int) value; 
       //printf("%d:%ld ", index, value);
       pos = end;
       index++;
   }

   //printf("Matrix size: %d x %d\n", A->rows, A->columns);
   free(buffer);
   MPI_File_close(&fh);
}

/*
    print_matrix(A)
    For debuging, print matrix A to console
*/

void print_matrix(matrix A)
{
   int i,j;
   printf("\nPrint matrix %d x %d\n", A.rows, A.columns);
   for(i = 0; i < A.rows; i++)
   {
       for(j = 0; j < A.columns; j++)
       {
   printf("%d ", A.matrix[i * A.columns + j]);
       }
       printf("\n");
   }
   printf("\n");
}

/*
    free_matrix(A)
    free the array in a matrix struct
*/

void free_matrix(matrix *A)
{
    if(A->matrix == NULL)
return;

    free(A->matrix);
}

/*
    send_matrix_partts(A, B, *A_rest, comm_size)
    Distribute matrix parts over the nodes
    Matrix B is send to every Node.
    Each Node receives a different part of A
    *A_rest gives the rest back to be processed by the root node
    1- define chunk size
    2- send matrix information to worker (row/column size)
    3- partition A and send to worker
    4- send B
    5- populate *A_rest
*/
void send_matrix_parts(matrix A, matrix B, matrix *A_rest, int comm_size)
{
    int chunk_size = ceil(A.rows / comm_size);

    MPI_Request request;
    int buffer[INFO_BUFFER];
    buffer[0] = chunk_size;
    buffer[1] = A.columns;
    buffer[2] = B.rows;
    buffer[3] = B.columns;

    int worker;
    int *A_pos = A.matrix;
    for(worker = 1; worker < comm_size; worker++)
    {
MPI_Isend(buffer, INFO_BUFFER, MPI_INT, worker, MATRIX_INFO, MPI_COMM_WORLD, 
&request);

        MPI_Isend(A_pos, chunk_size * A.columns, MPI_INT, worker, MATRIX_A, 
MPI_COMM_WORLD, &request);

MPI_Isend(B.matrix, B.rows * B.columns, MPI_INT, worker, MATRIX_B,
MPI_COMM_WORLD, &request);

A_pos += chunk_size * A.columns;
    }

    int A_rest_size = (A.matrix + (A.rows * A.columns)) - A_pos;
    //printf("rest size : %d", A_rest_size);
    A_rest->columns = A.columns;
    A_rest->rows = A_rest_size / A_rest->columns;
    A_rest->matrix = (int *) malloc(A_rest_size * sizeof(int));
    memcpy(A_rest->matrix, A_pos, A_rest_size * sizeof(int));
    
}

/*
    recv_matrix_parts(A, B)
    receive the chunk of matrix A, and the full matrix B from the master node
    1-receive matrix info for A and B
    2-initialize A and B
    3-receive A
    4-receive B
*/

void recv_matrix_parts(matrix *A, matrix *B)
{
   int buffer[INFO_BUFFER];
   MPI_Status status;
   MPI_Recv(buffer, INFO_BUFFER, MPI_INT, 0, MATRIX_INFO, MPI_COMM_WORLD, &status);
       
   //printf("A: %d x %d, B: %d x %d\n", buffer[0], buffer[1], buffer[2], buffer[3]);

   A->rows = buffer[0];
   A->columns = buffer[1];
   int A_size= A->rows * A->columns;
   A->matrix= (int *) malloc(A_size * sizeof(int));
   B->rows = buffer[2];
   B->columns = buffer[3];
   int B_size= B->rows * B->columns;
   B->matrix= (int *) malloc(B_size * sizeof(int));

   MPI_Recv(A->matrix, A_size , MPI_INT, 0, MATRIX_A, MPI_COMM_WORLD, &status);
   MPI_Recv(B->matrix, B_size , MPI_INT, 0, MATRIX_B, MPI_COMM_WORLD, &status);

}

/*
    send_matrix_result(C_part)
    send the result to the root node
*/

void send_matrix_result(matrix C_part)
{
    MPI_Send(C_part.matrix, C_part.rows * C_part.columns, MPI_INT, 0, MATRIX_C, MPI_COMM_WORLD);
}

/*
    recv_build_matrix(C_part, *C, comm_size)
    receive parts of the result matrix C from non root-nodes
    1-foreach worker
        1.1-receive C_part
        1.2-insert C_part in to C
    2-insert own C_part in to C
*/
void recv_build_matrix(matrix C_part, matrix *C, int comm_size)
{
    MPI_Status status;
    int worker;
    int *C_pos;
    int part_size = ceil(C->rows / comm_size) * C->columns;
    int C_size = C->rows * C->columns;
    C->matrix = (int *) malloc(C_size * sizeof(int));

    for(worker = 1, C_pos = C->matrix; worker < comm_size; worker++, C_pos += part_size)
    {
MPI_Recv(C_pos, part_size, MPI_INT, worker, MATRIX_C, MPI_COMM_WORLD, &status);
    }

    memcpy(C_pos, C_part.matrix, C_part.rows * C_part.columns * sizeof(int));
    
}


/*
    multiply_matrix(A, B, *C)
    basic matrix multiplication of A and B, result in C
*/

void multiply_matrix(matrix A, matrix B, matrix *C)
{
    C->rows = A.rows;
    C->columns = B.columns;
    int C_size = C->rows * C->columns;
    C->matrix= (int *) malloc(C_size * sizeof(int));

    //print_matrix(A);
    //print_matrix(B);

    int C_index;
    for(C_index = 0; C_index < C_size; C_index++)
    {
int i;
int value = 0;
int row = C_index / C->columns;
int column = C_index - (row * C->columns);

//printf("row: %d, column: %d\n", row, column);
for(i = 0; i < A.columns; i++)
{
    
    int A_index = row * A.columns + i;
    int B_index = i * B.columns + column;
    //printf("%d: %d,%d ", i,  A_index, B_index);

    value += A.matrix[A_index] * B.matrix[B_index];
}

C->matrix[C_index] = value;
//printf("\n");
    }
}

/*
    write_matrix(filenmae, C)
    write matrix C to File <filename>
*/
void write_matrix(char * filename, matrix C, double time)
{
    MPI_File_delete(filename, MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    int C_size = C.rows * C.columns;
    char *buffer = malloc((C_size * 10 +
      C_size +
  C.rows) *
  sizeof(char));

    MPI_Status status;
    int i;
    char *p;
    for(i = 0, p = buffer; i < C_size; i++)
    {
p += sprintf(p, "%d", C.matrix[i]);
if(i % C.columns == C.columns -1)
    p += sprintf(p, "\n");
else
    p += sprintf(p, " ");
    }

	// write calculation time
	p += sprintf(p, "\n\nCalculation time: %lf", time);
	
    MPI_File_write(fh, buffer, p - buffer, MPI_CHAR, &status);

    free(buffer);
    MPI_File_close(&fh);

}

/*
    main
    Init MPI
    if is root node (node == 0)
        read matrix A and B
        send matrix parts to worker
        calculate own part
        receive result
        build result matrix
        write matrix to File
    else
        receive matrix parts
        calculate part
        send result back to root
    End MPI
*/
int main(int argc, char **argv)
{
   int node;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &node);
   MPI_Status status;

   matrix A = nmatrix;
   matrix B = nmatrix;
   
   if(node == 0)
   {

       read_matrix(argv[1], &A);
       read_matrix(argv[2], &B);
       
       matrix A_rest = nmatrix;
       matrix C_part = nmatrix;
      
       int comm_size;
       MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
       
       printf("start calculation\n");
       double time = MPI_Wtime();

       send_matrix_parts(A, B, &A_rest, comm_size);
       multiply_matrix(A_rest, B, &C_part);

       matrix C = nmatrix;
       C.rows = A.rows;
       C.columns= B.columns;

       recv_build_matrix(C_part, &C, comm_size);
       
       time = MPI_Wtime() - time;
       printf("calculation on %d nodes: %.2f seconds\n", comm_size, time); 
       //print_matrix(C);
       write_matrix(argv[3], C, time);


       free_matrix(&A_rest);
       free_matrix(&C_part);
       free_matrix(&C);
   }    
   else
   {
       matrix C_part = nmatrix;
       recv_matrix_parts(&A, &B);
       multiply_matrix(A, B, &C_part);
       send_matrix_result(C_part);

       free_matrix(&C_part);
   }

   //printf("node: %d finalize\n", node);

   free_matrix(&A);
   free_matrix(&B);
   
   MPI_Finalize();
   return 0;
}
