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

#define DIST_THRESHOLD 1

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
vector acceleration(volatile point *a, volatile point *b)
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
void compute_movement(  volatile point *points, vector *point_vel, unsigned int offset,
                             unsigned int compute_size, unsigned int point_size, int node_id, int iteration) 
{
    unsigned int i, j;

    //compute new point_vel for each point between offset and offset + compute_size
    for(i = offset; i < offset + compute_size; i++)
    {
         volatile point *p1 = &points[i];
         // if weight == 0, point is absorbed and should be ignored
         if(p1->weight == 0) continue;
         
         vector acc = {0, 0};
         for(j = 0; j < point_size; j++)
         {
              volatile point *p2 = &points[j];
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
    for(i = offset; i < offset + compute_size; i++)
    {
        volatile point *p = &points[i];
        if(p->weight == 0){
            continue;
        }
        // apply movement
        p->x += point_vel[i - offset].x;
        p->y += point_vel[i - offset].y;
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

// The work done by Workers and Master
// for each iteration a new position is computed for the points, and send between the Nodes
// after all iterations are done exit
void work(int node_id, int comm_size, volatile point *points, int full_size, int iteration)
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
    
    init_vel(&point_vel, compute_size);
    
    // Main loop
    int i;
    for(i = 0; i < iteration; i++)
    {
        compute_movement(points, point_vel, offset, compute_size, full_size, node_id, i);
    }
        
    free(point_vel);
    //printf("%d end work\n", node_id);
    if(node_id == 0)
    printf("communication overhead: %.1f sec\n", overhead);
    
}

// write point back to output
void write_point( char *filename, volatile point *points, int size, double final_time, int iteration, int comm_size) 
{

    FILE *fp;
    fp = fopen(filename, "w");

    int i;
    for(i = 0; i < size; i++)
    {
point p = points[i];
fprintf(fp, "%.1f %.1f %.1f\n", p.x, p.y, p.weight);
    }
    
    fprintf(fp, "\n\n\nSimulation took: %.1f sec, for: %d iterations with: %d nodes\n", final_time, iteration, comm_size);

    fclose(fp);
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
    
    unsigned int sci_id;
    sci_query_adapter_t query;
    query.subcommand = SCI_Q_ADAPTER_NODEID;
    query.localAdapterNo = ADAPTER_NO;
    query.data = &sci_id;
    SCIQuery(SCI_Q_ADAPTER, &query, NO_FLAGS, &error);
    
    int iteration = (int) strtol(argv[1], NULL, 10);

    if(node_id == MASTER_ID)
    {
        sci_local_segment_t local_segment;
        volatile point *local_address;
        sci_map_t local_map;
        
        //printf("Master started\n");
         read_point(argv[2], &points, &full_size);
         
        unsigned int SEGMENT_SIZE = full_size * sizeof(point);
        printf("prepare segment: %d \n", SEGMENT_SIZE);
        
        SCICreateSegment(v_dev, &local_segment, SEGMENT_ID, SEGMENT_SIZE, NO_CALLBACK,
                   NO_ARG, NO_FLAGS, &error);
        if(error != SCI_ERR_OK)
           printf("Master error! %d\n", error);

        SCIPrepareSegment(local_segment, ADAPTER_NO, NO_FLAGS, &error);

        local_address = (volatile point *) SCIMapLocalSegment(local_segment, &local_map, 0, 
        SEGMENT_SIZE, 0, NO_FLAGS, &error);

        memcpy(local_address, points, SEGMENT_SIZE);
        free(points);

        SCISetSegmentAvailable(local_segment, ADAPTER_NO, NO_FLAGS, &error);

        // send metadata to worker SCI root and size of point array
        int metadata[2] = {sci_id, full_size};
        MPI_Bcast(metadata, 2, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

        // Take time
        double time = MPI_Wtime();

        // Main work 
        work(node_id, comm_size, local_address, full_size, iteration);

        // Final time
        double final_time = MPI_Wtime() - time;
        printf("Simulation took: %.1f sec, for: %d iterations with: %d nodes\n", final_time, iteration, comm_size);
        
        write_point(argv[3], local_address, full_size, final_time, iteration, comm_size);
    }
    else
    {
        unsigned int master_sci_id;
        sci_remote_segment_t remote_segment;
        sci_map_t remote_map;
        unsigned int segment_size;
        volatile point *remote_address;
        MPI_Status status;
        
        // receive meta data ( sci master id and size of point array )
        int metadata[2];
        MPI_Bcast(metadata, 2, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
        master_sci_id = metadata[0];
        full_size = metadata[1];

        printf("%d: received master_sci_id: %d\n", node_id, master_sci_id);

        SCIConnectSegment(v_dev, &remote_segment, master_sci_id, SEGMENT_ID, ADAPTER_NO,
            NO_CALLBACK, NO_ARG, SCI_INFINITE_TIMEOUT, NO_FLAGS, &error);
        
        if(error != SCI_ERR_OK)
         printf("error: %d !!! \n", error);

        segment_size = SCIGetRemoteSegmentSize(remote_segment);
        remote_address = (volatile point *) SCIMapRemoteSegment(remote_segment, 
            &remote_map, 0, segment_size, 0, NO_FLAGS, &error);
 
        work(node_id, comm_size, remote_address, full_size, iteration);
        
    }
    
    //printf("finalize\n");

    MPI_Finalize();
}
