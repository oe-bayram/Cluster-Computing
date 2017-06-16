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

#define MIN(a,b) (((a)<(b))?(a):(b))

#define MASTER_ID 0

#define DIST_THRESHOLD 5

// point struct
typedef struct
{
    float x;
    float y;
    float weight;
} point;

// vector struct, used for movement direction
typedef struct
{
    float x;
    float y;
} vector;


//  compute the acceleration induced by b to a
//  result vector = direction * magnitude
vector acceleration(point *a, point *b)
{
    float dx    = b->x - a->x;
    float dy    = b->y - a->y;
    float dist  = sqrt(dx*dx + dy*dy);

    //  If two points are close the bigger weight absorbs the other.
    //  The pointer address is used in case of a tie to decide which point absorbs the other.
    //  This way we are sure that every node will chose the same point.
    if(dist < DIST_THRESHOLD)
    {
         //printf("THRESHOLD\n");
         if(a->weight > b->weight || (a->weight == b->weight && a < b))
         {
              a->weight *= b->weight;
              b->weight = 0;
         }
         else
         {
              b->weight *= a->weight;
              a->weight = 0;
         }

         vector direction = {0, 0};
         return direction;
    }

    vector direction;
    direction.x = dx / dist;
    direction.y = dy / dist;

    float acc = b->weight / (dist * dist);

    direction.x *= acc;
    direction.y *= acc;

    //printf("weight: %.1f, dist: %.1f, computed acc: %.2f\n", b->weight, dist, acc);

    return direction;
}

// Add one vector to another
void actualise_vel(vector *vel, vector acc)
{
    //printf("acc: %.1f:%.1f\n", acc.x, acc.y);
    vel->x += acc.x;
    vel->y += acc.y;

    //printf("vel: %.1f:%.1f\n", vel->x, vel->y);
}

//  Update point position
//  Absorbed points ( weight == 0) are ignored
void compute_movement(  point *points, vector *point_vel, unsigned int offset,
                             unsigned int compute_size, unsigned int point_size, int *segment, int node_id) 
{
    unsigned int i, j;

    //compute new point_vel for each point between offset and offset + compute_size
    for(i = offset; i < offset + compute_size; i++)
    {
         point *p1 = &points[i];
// if weight == 0, point is absorbed and should be ignored
         if(p1->weight == 0) continue;
         
         vector acc = {0, 0};
         for(j = 0; j < point_size; j++)
         {
              point *p2 = &points[j];
    // if absorbed or same point ignore
              if (p2->weight == 0 || p1 == p2) continue;
              // add acceleration to total acceleration
              actualise_vel(&acc, acceleration(p1, p2));
         }
// add actual computed acceleration to velocity
         actualise_vel(&point_vel[i - offset], acc);
    }
    
    // Hier ein Barrier setzen, da Berechnungen der nächsten Iteration die aktuelle manipulieren würden
    MPI_Barrier(MPI_COMM_WORLD);

    // compute new point position
    printf("offset is: %d and offset + compute_size is: %d\n", offset, offset + compute_size);
    for(i = offset; i < offset + compute_size; i++)
    {
        printf("#######  %d: Iteration: %d  #######\n", node_id, i);
        point *p = &points[i];
         
        if(p->weight == 0) continue;
        // apply movement
        p->x += point_vel[i - offset].x;
        p->y += point_vel[i - offset].y;
         
        printf("point values are: %.1f %.1f %.1f\n", p->x, p->y, p->weight);
        // write new position to segment
        printf("offset is here: %u\n", i);
        
        write_point_segment(segment, p, i);
    }
}

// Read point from file
void read_point(char *filename, point **points, int *full_size)
{
    FILE *fp;
    long size;
    fp = fopen(filename, "r");

    fseek(fp, 0, SEEK_END);
    size = ftell(fp); // get current file pointer
    rewind(fp);

    char *buffer = (char *) malloc(size * sizeof(float));
    *points = (point *) malloc(size * sizeof(point));

    int result = fread(buffer, 1, size, fp);
    if(result != size) { perror("read error\n"); }

    char *pos = buffer;
    int index = 0;

    //printf("read: %s\n", buffer);

    while (pos < buffer + size)
    {
         char *end;
         float x = strtod(pos, &end);
         pos = end;
         float y = strtod(pos, &end);
         pos = end;
         float weight = strtod(pos, &end);

         if (end == pos)
              break;
         
         point p = { x, y, weight };
//printf("x: %.1f, y: %.1f, weight: %.1f\n", x, y, weight);
         (*points)[index] = p;

         pos = end;
         index++;
    }

    *full_size = index;

    free(buffer);
    fclose(fp);

    //printf("end read points\n");
}

// Send point from Master to Worker node
// 1 send: number of points
// 2 send: full point array
// Each struct is 3 float, so we send 3 * size * MPI_FLOAT
void send_point(point *points, int size)
{
    //printf("start send point\n");
    MPI_Bcast(&size, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(points, size * 3, MPI_FLOAT, MASTER_ID, MPI_COMM_WORLD);    
}

//  Initialize velocity vector
void init_vel(vector **point_vel, int size)
{
    *point_vel = (vector *) malloc(size * sizeof(vector));

    int i;
    for(i = 0; i < size; i++)
    {
         vector v = {0, 0};        // null velocity at start
         (*point_vel)[i] = v;
    }
}

// Receive initialisation points from Master
void receive_points(point **points, int *full_size)
{
    int size;
    // receive broadcast
    MPI_Bcast(&size, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    // alocate array
    *points = (point *) malloc(size * sizeof(point));
    MPI_Bcast(*points, size * 3, MPI_FLOAT, MASTER_ID, MPI_COMM_WORLD); // Each point consist of 3 float

    *full_size = size;

    //printf("worker received points\n");
}

// Every Node Bcast its part to the other. From 0 to n Nodes
void update_points(int comm_size, point *points, int size)
{

    size *= 3;// Each point consist of 3 float
    int chunk =  (size / comm_size);
    int node_id;
    for(node_id = 0; node_id < comm_size; node_id++)
    {
         int start    = chunk * node_id;
         
        if(node_id == comm_size -1)
            chunk = size - start;

        printf("node_id: %d, start: %d, size: %d, count: %d\n", node_id, start, size, chunk);
        MPI_Bcast(points + (start/3), chunk, MPI_FLOAT, node_id, MPI_COMM_WORLD);
    }
}

// The work done by Workers and Master
// for each iteration a new position is computed for the points, and send between the Nodes
// after all iterations are done exit
void work(int node_id, int comm_size, point *points, int full_size, int iteration, int *segment)
{
    // overhead counter
    double overhead = 0;
    // point velocities, allocated only of size N/n
    vector *point_vel;
    // how much points this Node need to compute
    int compute_size;
    int offset;
    
    compute_size = full_size / comm_size;
    offset = node_id * compute_size;

    // ceil
    if(node_id == comm_size -1)
    compute_size = (full_size + comm_size -1 ) / comm_size ;
    
    //printf("full_size: %d, comm_size: %d, node_id: %d, compute_size: %d\n", full_size, comm_size, node_id, compute_size);

    init_vel(&point_vel, compute_size);
    
    // Main loop
    int i;
    for(i = 0; i < iteration; i++)
    {
         compute_movement(points, point_vel, offset, compute_size, full_size, segment, node_id);
    }
    

    free(point_vel);
    //printf("%d end work\n", node_id);
    if(node_id == 0)
    printf("communication overhead: %.1f sec\n", overhead);
    
}

// write point back to output
void write_point( char *filename, point *points, int size) 
{

    FILE *fp;
    fp = fopen(filename, "w");

    int i;
    for(i = 0; i < size; i++)
    {
point p = points[i];
fprintf(fp, "%.1f %.1f %.1f\n", p.x, p.y, p.weight);
    }

    fclose(fp);
}

/*
    write_point_segment(local_address, A, B)
    write matrixes A and B and their sizes to segment
*/
void write_point_segment(int *segment, point *p, unsigned int offset)
{
    printf("got this point: %.1f %.1f %.1f and offset is: %u\n", p->x, p->y, p->weight, offset);
    int *pos = segment;
    pos += 1;
    pos += offset;
    memcpy(pos, p, 1 * sizeof(point));
}

/*
    write_points_segment(local_address, A, B)
    write matrixes A and B and their sizes to segment
*/
void write_points_segment(int *segment, point *points, int full_size)
{
    int *pos = segment;
    segment[0] = full_size;
    pos += 1;
    memcpy(pos, points, full_size * sizeof(point));
}

/*
    write_points_segment(local_address, A, B)
    write matrixes A and B and their sizes to segment
*/
void read_points_segment(int *segment, point **points, int *full_size)
{
    int *pos = segment;
    int size = segment[0];
    *full_size = size;
    pos += 1;
    *points = (point *) malloc(size * sizeof(point));
    memcpy(*points, pos, size * sizeof(point));
}

/*
Master:
- read points from file
- send point to worker
- Work
- write point back to file
Worker:
- receive points from Master
- Work
Work:
- compute new position
- share result
- repeat
*/

int main(int argc, char **argv)
{
    //printf("START MAIN\n");
    int node_id;              // MPI rank from 0 to n nodes
    int comm_size;            // number of nodes
    point *points;            // array holding all the points
    int full_size;            // number of points


    //  MPI initialisation
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    sci_desc_t v_dev;
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
    
    int iteration = (int) strtol(argv[1], NULL, 10);

    if(node_id == MASTER_ID)
    {
        sci_local_segment_t local_segment;
        int *local_address;
        sci_map_t local_map;
        
        //printf("Master started\n");
         read_point(argv[2], &points, &full_size);
         
        unsigned int SEGMENT_SIZE = full_size;
        printf("prepare segment: %d \n", SEGMENT_SIZE);
        
        SCICreateSegment(v_dev, &local_segment, SEGMENT_ID, SEGMENT_SIZE, NO_CALLBACK,
                   NO_ARG, NO_FLAGS, &error);
        if(error != SCI_ERR_OK)
           printf("Master error! %d\n", error);

        SCIPrepareSegment(local_segment, ADAPTER_NO, NO_FLAGS, &error);

        local_address = (int *) SCIMapLocalSegment(local_segment, &local_map, 0, 
        SEGMENT_SIZE, 0, NO_FLAGS, &error);
        
        write_points_segment(local_address, points, full_size);
        
        // Check if written to segment
        point *p1 = &local_address[1];
        printf("First coordinate is: %.1f %.1f %.1f\n", p1->x, p1->y, p1->weight);
        
        SCISetSegmentAvailable(local_segment, ADAPTER_NO, NO_FLAGS, &error);

        // send segment information to other nodes
        MPI_Bcast(&local_node_id, 1, MPI_INT, node_id, MPI_COMM_WORLD);

        // Take time
        double time = MPI_Wtime();

        // Main work
        // send_point(points, full_size);
        work(node_id, comm_size, points, full_size, iteration, local_address);

        // Final time
        double final_time = MPI_Wtime() - time;
        printf("Simulation took: %.1f sec, for: %d iterations with: %d nodes\n", final_time, iteration, comm_size);
        //read_points_segment(local_address, &points, &full_size);
        write_point(argv[3], points, full_size);
    }
    else
    {
        unsigned int master_node_id;
        sci_remote_segment_t remote_segment;
        sci_map_t remote_map;
        unsigned int segment_size;
        volatile int *remote_address;
        MPI_Status status;
        MPI_Bcast(&master_node_id, 1, MPI_INT, 0, MPI_COMM_WORLD);

        printf("%d: received master_node_id: %d\n", node_id, master_node_id);

        SCIConnectSegment(v_dev, &remote_segment, master_node_id, SEGMENT_ID, ADAPTER_NO,
            NO_CALLBACK, NO_ARG, SCI_INFINITE_TIMEOUT, NO_FLAGS, &error);
        
        if(error != SCI_ERR_OK)
         printf("error: %d !!! \n", error);
        segment_size = SCIGetRemoteSegmentSize(remote_segment);
        remote_address = (volatile int *) SCIMapRemoteSegment(remote_segment, 
            &remote_map, 0, segment_size, 0, NO_FLAGS, &error);
        
        read_points_segment(remote_address, &points, &full_size);
        work(node_id, comm_size, points, full_size, iteration, remote_segment);
    }
    
    //printf("finalize\n");
    free(points);

    MPI_Finalize();
}
