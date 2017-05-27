#include <sisci_api.h>
#include <sisci_error.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NO_FLAGS 0
#define NO_CALLBACK 0
#define NO_ARG 0

#define SEGMENT_ID 126
#define ADAPTER_NO 0

#define MAX(a,b) (((a)>(b))?(a):(b))
#define INFO_BUFFER 4
#define MATRIX_INFO 0                   //Tag for MPI_Send
#define MATRIX_A    1                   //Tag for MPI_Send
#define MATRIX_B    2                   //Tag for MPI_Send
#define MATRIX_C    3                   //Tag for MPI_Send
#define MASTER_NODE 5                   //Tag for MPI_Send


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

   char *buffer		= (char *) malloc(size * sizeof(char));;
   A->rows 		= 0;
   A->columns 		= 0;

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
    A_rest->columns 	= A.columns;
    A_rest->rows 	= A_rest_size / A_rest->columns;
    A_rest->matrix 	= (int *) malloc(A_rest_size * sizeof(int));
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

   A->rows 	= buffer[0];
   A->columns 	= buffer[1];
   int A_size	= A->rows * A->columns;
   A->matrix	= (int *) malloc(A_size * sizeof(int));
   B->rows 	= buffer[2];
   B->columns 	= buffer[3];
   int B_size	= B->rows * B->columns;
   B->matrix	= (int *) malloc(B_size * sizeof(int));

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
    int part_size 	= ceil(C->rows / comm_size) * C->columns;
    int C_size 		= C->rows * C->columns;
    C->matrix 		= (int *) malloc(C_size * sizeof(int));

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
    C->rows 	= A.rows;
    C->columns 	= B.columns;
    int C_size 	= C->rows * C->columns;
    C->matrix	= (int *) malloc(C_size * sizeof(int));

    //print_matrix(A);
    //print_matrix(B);

    int C_index;
    for(C_index = 0; C_index < C_size; C_index++)
    {
	int i;
	int value 	= 0;
	int row 	= C_index / C->columns;
	int column 	= C_index - (row * C->columns);

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
void write_matrix(char * filename, matrix C)
{
    MPI_File_delete(filename, MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    int C_size 		= C.rows * C.columns;
    char *buffer 	= malloc((C_size * 10 +
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

    MPI_File_write(fh, buffer, p - buffer, MPI_CHAR, &status);

    free(buffer);
    MPI_File_close(&fh);
}

/*
    write_raw_matrix_segment(local_address, A, B)
    write matrixes A and B and their sizes to segment
*/
void write_raw_matrix_segment(int *segment, matrix A, matrix B)
{
	int *pos = segment;
	segment[0] = A.rows;
	segment[1] = A.columns;
	segment[2] = B.rows;
	segment[3] = B.columns;
	pos += 4;
	memcpy(pos, A.matrix, A.rows * A.columns * sizeof(int));
	pos += A.rows * A.columns;
	memcpy(pos, B.matrix, B.rows * B.columns * sizeof(int));
}

/*
    write_result_segment(segment, C, comm_size, node_id)
    write result (matrix C) to segment depending on node id
*/
void write_result_segment(int *segment, matrix C, int comm_size, int node_id)
{
	int chunk_size = ceil(segment[0] / comm_size);
	int *C_pos = segment;
	C_pos += 4 + segment[0] * segment[1] + segment[2] * segment[3]; // set position to begin of segment part for C
	C_pos += (node_id-1) * chunk_size * segment[3]; // set position depending on calculated part
	memcpy(C_pos, C.matrix, C.rows * C.columns * sizeof(int));
}

/*
    read_matrix_segment(segment, A, position)
    read matrix from segment (matrix A)
*/
void read_matrix_segment(int *segment, matrix *A, int position)
{
	int *A_pos = segment;
	A_pos += position;
	A.matrix = (int *) malloc(A.rows * A.columns * sizeof(int));
	memcpy(A.matrix, A_pos, A.rows * A.columns * sizeof(int));
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
   printf("main started\n");
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &node);
   MPI_Status status;
   
   int comm_size;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
   sci_desc_t	v_dev;
   sci_error_t error;

   SCIInitialize(NO_FLAGS, &error);
   printf("SCI initialized\n");
   if(error != SCI_ERR_OK)
   printf("error Init\n");

   SCIOpen(&v_dev, NO_FLAGS, &error);
   printf("SCI opened\n");
   if(error != SCI_ERR_OK) {
   printf("Error Open\n");
      return 1;
   }

   unsigned int local_node_id;
   sci_query_adapter_t query;
   query.subcommand = SCI_Q_ADAPTER_NODEID;
   query.localAdapterNo = ADAPTER_NO;
   query.data = &local_node_id;
   SCIQuery(SCI_Q_ADAPTER, &query, NO_FLAGS, &error);

   matrix A = nmatrix;
   matrix B = nmatrix;
   matrix C = nmatrix;
   
   if(node == 0) {
      sci_local_segment_t local_segment;
      int *local_address;
      sci_map_t local_map;
      printf("Reading matrix A\n");
      read_matrix(argv[1], &A);
	  print_matrix(A);
	  printf("Reading matrix B\n");
      read_matrix(argv[2], &B);
	  print_matrix(B);
      matrix A_rest 	= nmatrix;
      matrix C_part 	= nmatrix;

      printf("start calculation\n");
      double time = MPI_Wtime();

      // send_matrix_parts(A, B, &A_rest, comm_size);
      int A_size = A.columns * A.rows;
      int B_size = B.columns * B.rows;
      int C_size = B.columns * A.rows;
      C.rows = A.rows;
      C.columns = B.columns;
	   
      unsigned int SEGMENT_SIZE = A_size+B_size+C_size+4;
	   
      printf("prepare segment: %d \n", SEGMENT_SIZE);
	   
      SCICreateSegment(v_dev, &local_segment, SEGMENT_ID, SEGMENT_SIZE, NO_CALLBACK,
				NO_ARG, NO_FLAGS, &error);
      if(error != SCI_ERR_OK)
         printf("Master error! %d\n", error);

      SCIPrepareSegment(local_segment, ADAPTER_NO, NO_FLAGS, &error);

      local_address = (int *) SCIMapLocalSegment(local_segment, &local_map, 0, 
      SEGMENT_SIZE, 0, NO_FLAGS, &error);
	  
	  write_raw_matrix_segment(local_address, A, B);
		
      SCISetSegmentAvailable(local_segment, ADAPTER_NO, NO_FLAGS, &error);

      // send segment information to other nodes
      MPI_Bcast(&local_node_id, 1, MPI_INT, node, MPI_COMM_WORLD);

	  int chunk_size = ceil(A.rows / comm_size);
	  A_rest.rows = A.rows - (comm_size-1) * chunk_size;
	  A_rest.columns = A.columns;
	  
	  int *A_pos = A.matrix;
	  A_pos += (comm_size-1) * chunk_size * A.columns;
	  A_rest.matrix = (int *) malloc(A_rest.rows * A.columns * sizeof(int));
	  memcpy(A_rest.matrix, A_pos, A_rest.rows * A.columns * sizeof(int));
	  
      multiply_matrix(A_rest, B, &C_part);
	  print_matrix(A_rest);
	  print_matrix(C_part);

	  write_result_segment(local_address, C_part, comm_size, comm_size);
	  
      MPI_Barrier(MPI_COMM_WORLD);
	  
	  matrix C = nmatrix;
	  C.rows = A.rows;
      C.columns	= B.columns;
	  int *C_end_pos = local_address;
	  C_end_pos += 4 + A_size + B_size;
	  C.matrix = (int *) malloc(C_size * sizeof(int));
	  memcpy(C.matrix, C_end_pos, C_size * sizeof(int));
      time = MPI_Wtime() - time;
      printf("calculation on %d nodes: %.2f seconds\n", comm_size, time); 
      print_matrix(C);
	      
      write_matrix(argv[3], C);

      free_matrix(&A_rest);
      free_matrix(&C_part);
      free_matrix(&C);
   } else {
      matrix C_part = nmatrix;
      //recv_matrix_parts(&A, &B);
      unsigned int master_node_id;
      sci_remote_segment_t remote_segment;
      sci_map_t remote_map;
      unsigned int segment_size;
      volatile int *remote_address;
      MPI_Status status;
      MPI_Bcast(&master_node_id, 1, MPI_INT, 0, MPI_COMM_WORLD);

      printf("received master_node_id: %d\n", master_node_id);

      SCIConnectSegment(v_dev, &remote_segment, master_node_id, SEGMENT_ID, ADAPTER_NO,
		    NO_CALLBACK, NO_ARG, SCI_INFINITE_TIMEOUT, NO_FLAGS, &error);
			
      if(error != SCI_ERR_OK)
         printf("error: %d !!! \n", error);

      segment_size = SCIGetRemoteSegmentSize(remote_segment);
	
      remote_address = (volatile int *) SCIMapRemoteSegment(remote_segment, 
		&remote_map, 0, segment_size, 0, NO_FLAGS, &error);

      printf("Node: %d, first two array value: %d, %d\n", local_node_id,
      remote_address[0], remote_address[1]);
	  
	  A.rows = ceil(remote_address[0] / comm_size);
	  A.columns = remote_address[1];
	  B.rows = remote_address[2];
	  B.columns = remote_address[3];

	  int position = (node-1) * A.rows * A.columns + 4;
	  read_matrix_segment(segment, &A, position);
	  
	  int *B_pos = remote_address;
	  B_pos += 4 + remote_address[0] * remote_address[1];
	  B.matrix = (int *) malloc(B.rows * B.columns * sizeof(int));
	  memcpy(B.matrix, B_pos, B.rows * B.columns * sizeof(int));
	  
	  print_matrix(A);
      multiply_matrix(A, B, &C_part);
	  write_result_segment(remote_address, C_part, comm_size, node);
	  
      MPI_Barrier(MPI_COMM_WORLD);
      free_matrix(&C_part);
   }

   //printf("node: %d finalize\n", node);

   free_matrix(&A);
   free_matrix(&B);
   
   SCIClose(v_dev, NO_FLAGS, &error);
   SCITerminate();
   MPI_Finalize();
   return 0;
}