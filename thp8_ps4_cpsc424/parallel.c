#include "timing.h"
#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"

#define ROOT 0
#define DEBUG 0

double LENGTH(double x, double y, double z){
  return (sqrt((x * x) + (y * y) + (z * z)));
}
double POSDIST(double x){
  if (x > 0) return 0; else return ((-1) * x);
}
double NEGDIST(double x){
  if (x < 0) return 0; else return x;
}


typedef struct body {
  double mass;

  double position_x;
  double position_y;
  double position_z;
  
  double velocity_x;
  double velocity_y;
  double velocity_z;
} body;

MPI_Datatype body_type;

double distance(body a, body b){
  double dx, dy, dz;

  dx = b.position_x - a.position_x;
  dy = b.position_y - a.position_y;
  dz = b.position_z - a.position_z;

  return sqrt(dx*dx + dy*dy + dz*dz);
}

double force_x(body a, body b){
  double dist = distance(a, b);

  return (a.mass * b.mass * (b.position_x - a.position_x))/(dist * dist * dist);
}

double force_y(body a, body b){
  double dist = distance(a, b);

  return (a.mass * b.mass * (b.position_y - a.position_y))/(dist * dist * dist);
}

double force_z(body a, body b){
  double dist = distance(a, b);

  return (a.mass * b.mass * (b.position_z - a.position_z))/(dist * dist * dist);
}

double stepVelocity_x(body b, double forc_x, double DT){
  return b.velocity_x + (DT * (forc_x / b.mass));
}

double stepVelocity_y(body b, double forc_y, double DT){
  return b.velocity_y + (DT * (forc_y / b.mass));
}

double stepVelocity_z(body b, double forc_z, double DT){
  return b.velocity_z + (DT * (forc_z / b.mass));
}

double stepCoords_x(body b, double forc_x, double DT){
  double avgVelocity_x = b.velocity_x + (0.5 * DT * (forc_x / b.mass));
  return b.position_x + avgVelocity_x * DT;
}

double stepCoords_y(body b, double forc_y, double DT){
  double avgVelocity_y = b.velocity_y + (0.5 * DT * (forc_y / b.mass));
  return b.position_y + avgVelocity_y * DT;
}

double stepCoords_z(body b, double forc_z, double DT){
  double avgVelocity_z = b.velocity_z + (0.5 * DT * (forc_z / b.mass));
  return b.position_z + avgVelocity_z * DT;
}

void defineBodyType(){
  const int nitems = 7;
  int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};

  MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets[7];

  offsets[0] = offsetof(body, mass);
  offsets[1] = offsetof(body, position_x);
  offsets[2] = offsetof(body, position_y);
  offsets[3] = offsetof(body, position_z);
  offsets[4] = offsetof(body, velocity_x);
  offsets[5] = offsetof(body, velocity_y);
  offsets[6] = offsetof(body, velocity_z);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &body_type);
  MPI_Type_commit(&body_type);
}

double distToQuad(body b, int quad){
  if(quad < 0 || quad > 7){
    fprintf(stderr, "Bad quadrant: %d\n", quad);
    exit(-1);
  }
  
  double distance = -1.0;
  
  switch(quad){
    case 0: distance = LENGTH(POSDIST(b.position_x),
                POSDIST(b.position_y),
                POSDIST(b.position_z));
            break;
    case 1: distance = LENGTH(NEGDIST(b.position_x),
                POSDIST(b.position_y),
                POSDIST(b.position_z));
            break;
    case 2: distance = LENGTH(NEGDIST(b.position_x),
                NEGDIST(b.position_y),
                POSDIST(b.position_z));
            break;
    case 3: distance = LENGTH(POSDIST(b.position_x),
                NEGDIST(b.position_y),
                POSDIST(b.position_z));
            break;
    case 4: distance = LENGTH(POSDIST(b.position_x),
                POSDIST(b.position_y),
                NEGDIST(b.position_z));
            break;
    case 5: distance = LENGTH(NEGDIST(b.position_x),
                POSDIST(b.position_y),
                NEGDIST(b.position_z));
            break;
    case 6: distance = LENGTH(NEGDIST(b.position_x),
                NEGDIST(b.position_y),
                NEGDIST(b.position_z));
            break;
    case 7: distance = LENGTH(POSDIST(b.position_x),
                NEGDIST(b.position_y),
                NEGDIST(b.position_z));
            break;
    default:
            fprintf(stderr, "Programming error, shouldn't reach here\n");
            exit(-2);
  }

  return distance;
}

void tests();
void master_setup();
void worker_process(int rank);
void worker_receive(int rank);

body *bodies;
int howManyWeHave, N, K;
double DT;

int main(int argc, char **argv) {
  int size, rank;
  double wctime, wctime_end, cputime;
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&size);

  defineBodyType(); /* define the body MPI type */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == ROOT) {
    tests();
    timing(&wctime, &cputime);
    master_setup();
    worker_process(rank);
    timing(&wctime_end, &cputime);
    printf("Time for %d timesteps with %d bodies: %lf\n", K, N, wctime_end - wctime);
  }else{
    worker_receive(rank);
    worker_process(rank);
  }
  
  MPI_Finalize();
}

/*
 * Read in all data, send it to all the processes, and figure out
 * how to distribute bodies between processes
 */
void master_setup() {
  
  // Part 1 - Read in all the data from stdin!
  
  body *inputBodies, *sortedBodies;
  int i, j, c;
  
  scanf("%d\n%d\n%lf", &N, &K, &DT);

  inputBodies = calloc(N, sizeof(body));
  sortedBodies = calloc(N, sizeof(body));
  bodies = calloc(N, sizeof(body));

  for(i = 0; i < N; i++){
    scanf("%lf\n", &(inputBodies[i].mass));
  }
  for(i = 0; i < N; i++){
    scanf("%lf %lf %lf\n", &(inputBodies[i].position_x), &(inputBodies[i].position_y), &(inputBodies[i].position_z)); 
  }
  for(i = 0; i < N; i++){
    scanf("%lf %lf %lf\n", &(inputBodies[i].velocity_x), &(inputBodies[i].velocity_y), &(inputBodies[i].velocity_z)); 
  }
  if (DEBUG) printf("Master done reading in bodies! N = %d\n", N);

  // Part 2 - sort the array by assignment to processor! (must send basic info first)

  MPI_Bcast(&N, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&K, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&DT, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  int *taken = calloc(N, sizeof(int)); // TODO free this!
  int offsets[8] = {0};
  int counts[8] = {0};

  c = 0;
  for(i = 0; i < 8; i++){
    offsets[i] = c;
    for(j = 0; j < N; j++){
      if(distToQuad(inputBodies[j], i) == 0){ // belongs to processor i
        if(taken[j]){
          fprintf(stderr, "Error, element %d assigned twice!\n", j);
          exit(-3);
        }else{
          taken[j] = 1; // mark as taken!
        }

        sortedBodies[c] = inputBodies[j];
        counts[i]++;
        c++;
      }
    }
    if(DEBUG) printf("Master sending %d bodies to processor %d\n", counts[i], i);
  }

  // quick sanity check!
  for(i = 0; i < N; i++){
    if(taken[i] != 1){
      fprintf(stderr, "Error, element %d not assigned to any processor!\n", i);
      exit(-4);
    }
  }
  if(c != N){
    fprintf(stderr, "Error, not all elements assigned! c = %d, N = %d\n", c, N);
    exit(-5);
  }

  // ... ok we're all good!
  //
  // Part 3 - send the assignments to each processor!

  MPI_Bcast(offsets, 8, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(counts, 8, MPI_INT, ROOT, MPI_COMM_WORLD);

  MPI_Scatterv(sortedBodies, counts, offsets, body_type, bodies, counts[ROOT], body_type, ROOT, MPI_COMM_WORLD);

  howManyWeHave = counts[ROOT]; // number of bodies assigned to master
}

void worker_receive(int rank){
  MPI_Bcast(&N, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&K, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&DT, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  bodies = calloc(N, sizeof(body));
  
  int offsets[8] = {0};
  int counts[8] = {0};
  
  MPI_Bcast(offsets, 8, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(counts, 8, MPI_INT, ROOT, MPI_COMM_WORLD);

  MPI_Scatterv(NULL, counts, offsets, body_type, bodies, counts[rank], body_type, ROOT, MPI_COMM_WORLD);

  if (DEBUG) printf("Process %d received %d elements!\n", rank, counts[rank]);
  for(int i = 0; i < counts[rank]; i++){
    // printf("\tProcess %d received element with coords (%e, %e, %e)\n", rank, bodies[i].position_x, bodies[i].position_y, bodies[i].position_z);
  }
  howManyWeHave = counts[rank];
}

void postUpdate(int rank, int round){
  int i;
  double numerator_x = 0.0;
  double numerator_y = 0.0;
  double numerator_z = 0.0;
  
  double numerator_x_total = 0.0;
  double numerator_y_total = 0.0;
  double numerator_z_total = 0.0;
  
  double numerator_x_velocity = 0.0;
  double numerator_y_velocity = 0.0;
  double numerator_z_velocity = 0.0;
  
  double numerator_x_total_velocity = 0.0;
  double numerator_y_total_velocity = 0.0;
  double numerator_z_total_velocity = 0.0;
  
  double denominator = 0.0;
  double denominator_total = 0.0;

  int quadrantList[8] = {0};

  MPI_Gather(&howManyWeHave, 1, MPI_INT, quadrantList, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  for(i = 0; i < howManyWeHave; i++){
    numerator_x += bodies[i].mass * bodies[i].position_x;
    numerator_y += bodies[i].mass * bodies[i].position_y;
    numerator_z += bodies[i].mass * bodies[i].position_z;

    numerator_x_velocity += bodies[i].velocity_x;
    numerator_y_velocity += bodies[i].velocity_y;
    numerator_z_velocity += bodies[i].velocity_z;

    denominator += bodies[i].mass;
  }

  MPI_Reduce(&numerator_x, &numerator_x_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&numerator_y, &numerator_y_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&numerator_z, &numerator_z_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  
  MPI_Reduce(&numerator_x_velocity, &numerator_x_total_velocity, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&numerator_y_velocity, &numerator_y_total_velocity, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&numerator_z_velocity, &numerator_z_total_velocity, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

  MPI_Reduce(&denominator, &denominator_total, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  
  if(rank == ROOT){
    printf("Conditions after timestep %d (time = %e):\n\n", round, round * DT);
    printf("\tCenter of Mass: (%e, %e, %e)\n", 
        numerator_x_total / denominator_total, numerator_y_total / denominator_total, numerator_z_total / denominator_total);
    printf("\tAverage Velocity: (%e, %e, %e)\n\n\n", 
        numerator_x_total_velocity / N, numerator_y_total_velocity / N, numerator_z_total_velocity / N);
    printf("\tQuadrant breakdown (0-7): ");
    for(i = 0; i < 8; i++){
      printf("%d, ", quadrantList[i]);
    }
    printf("\n\n\n");
  }
}

void worker_process(int rank){
  int round = 0;
  for(round = 0; round < K; round++){

    // Part 0 - Post update if necessary!
    if(round % 128 == 0){
      postUpdate(rank, round);
    }

    // Part 1 - Advertise to each neighbor how many bodies he'll be notified of!
    int i, j;

    int agentsToAdvertiseCount = 0;
    int agentsToAdvertiseCounts[8] = {0};
    int agentsToAdvertiseOffsets[8] = {0};
    body *agentsToAdvertise = calloc(howManyWeHave * 7, sizeof(body)); // TODO make this smaller...

    for(i = 0; i < 8; i++){
      if (i == rank) continue;

      agentsToAdvertiseOffsets[i] = agentsToAdvertiseCount;
      for(j = 0; j < howManyWeHave; j++){
        if(distToQuad(bodies[j], i) < 5.0){ // gotta notify processor i!
          agentsToAdvertise[agentsToAdvertiseCount] = bodies[j];
          agentsToAdvertiseCounts[i]++;
          agentsToAdvertiseCount++;
        }
      }

      if(DEBUG) printf("Processor %d is notifying processor %d of %d bodies!\n", rank, i, agentsToAdvertiseCounts[i]);
    }

    int agentsToBeAdvertisedOf[8] = {0};

    MPI_Alltoall(agentsToAdvertiseCounts, 1, MPI_INT, agentsToBeAdvertisedOf, 1, MPI_INT, MPI_COMM_WORLD);

    // quick sanity check!
    if(agentsToBeAdvertisedOf[rank] != 0){
      fprintf(stderr, "Somehow rank %d is getting notified %d elements...\n", rank, agentsToBeAdvertisedOf[rank]);
      exit(-7);
    }

    for(i = 0; i < 8; i++){
      // printf("Process %d is getting %d notifications from processor %d\n", rank, agentsToBeAdvertisedOf[i], i);
    }
    
    // ... ok, done!
    //
    // Part 3 - Doing the advertising!

    int totalToReceive = 0;
    for(i = 0; i < 8; i++){
      totalToReceive += agentsToBeAdvertisedOf[i];
    }

    body *advertisements = calloc(totalToReceive, sizeof(body));
    int agentsToBeAdvertisedOfOffsets[8] = {0};
    int agentsToBeAdvertisedOfOffsetCounter = 0;
    for(i = 0; i < 8; i++){
      agentsToBeAdvertisedOfOffsets[i] = agentsToBeAdvertisedOfOffsetCounter;
      agentsToBeAdvertisedOfOffsetCounter += agentsToBeAdvertisedOf[i];
    }
    
    MPI_Alltoallv(agentsToAdvertise, agentsToAdvertiseCounts, agentsToAdvertiseOffsets, body_type,
        advertisements, agentsToBeAdvertisedOf, agentsToBeAdvertisedOfOffsets, body_type, MPI_COMM_WORLD);

    // sanity check!
    for(i = 0; i < totalToReceive; i++){
      if(distToQuad(advertisements[i], rank) > 5.0){
        fprintf(stderr, "Error, processor %d was notified of a body more than 5 DUs away\n", rank);
        exit(-7);
      }else if(distToQuad(advertisements[i], rank) == 0.0){
        fprintf(stderr, "Error, processor %d was notified of a body in its own quadrant\n", rank);
        exit(-8);
      }
    }

    // Part 4 - doing the math stepping!!
    
    double *new_x = calloc(howManyWeHave, sizeof(double));
    double *new_y = calloc(howManyWeHave, sizeof(double));
    double *new_z = calloc(howManyWeHave, sizeof(double));

    double *new_velocity_x = calloc(howManyWeHave, sizeof(double));
    double *new_velocity_y = calloc(howManyWeHave, sizeof(double));
    double *new_velocity_z = calloc(howManyWeHave, sizeof(double));
    
    for(i = 0; i < howManyWeHave; i++){
      double forc_x, forc_y, forc_z;
      forc_x = 0; forc_y = 0; forc_z = 0;

      for(j = 0; j < howManyWeHave; j++){
        if(i == j) continue;

        double dist = distance(bodies[i], bodies[j]);

        if(dist > 5.0) continue;

        forc_x += force_x(bodies[i], bodies[j]);
        forc_y += force_y(bodies[i], bodies[j]);
        forc_z += force_z(bodies[i], bodies[j]);
      }

      for(j = 0; j < totalToReceive; j++){
        double dist = distance(bodies[i], advertisements[j]);
        
        if(dist > 5.0) continue;
        
        forc_x += force_x(bodies[i], advertisements[j]);
        forc_y += force_y(bodies[i], advertisements[j]);
        forc_z += force_z(bodies[i], advertisements[j]);
      }

      new_velocity_x[i] = stepVelocity_x(bodies[i], forc_x, DT);
      new_velocity_y[i] = stepVelocity_y(bodies[i], forc_y, DT);
      new_velocity_z[i] = stepVelocity_z(bodies[i], forc_z, DT);

      new_x[i] = stepCoords_x(bodies[i], forc_x, DT);
      new_y[i] = stepCoords_y(bodies[i], forc_y, DT);
      new_z[i] = stepCoords_z(bodies[i], forc_z, DT);
    }

    for(i = 0; i < howManyWeHave; i++){
      if(DEBUG) printf("\tProcessor %d is moving a body from x = %e to %e, velocity x = %e to %e\n", 
          rank, bodies[i].position_x, new_x[i], bodies[i].velocity_x, new_velocity_x[i]);
      bodies[i].position_x = new_x[i];
      bodies[i].position_y = new_y[i];
      bodies[i].position_z = new_z[i];

      bodies[i].velocity_x = new_velocity_x[i];
      bodies[i].velocity_y = new_velocity_y[i];
      bodies[i].velocity_z = new_velocity_z[i];
    }

    // Part 5 - the math stepping is done! - Now we need to redistribute the bodies

    int numLost = 0;
    body *bodiesToSend = calloc(howManyWeHave, sizeof(body));

    for(i = 0; i < howManyWeHave; i++){
      if(distToQuad(bodies[i], rank) != 0.0){
        if(DEBUG) printf("Body at position (%e, %e, %e) left quadrant %d\n", bodies[i].position_x, bodies[i].position_y, bodies[i].position_z, rank);
        bodiesToSend[numLost] = bodies[i];
        numLost++;
      }
    }

    // we need to remove the old ones we don't need anymore (they're saved in bodiesToSend)
    int index = 0;
    for(i = 0; i < howManyWeHave; i++){
      if(distToQuad(bodies[i], rank) == 0.0){
        bodies[index] = bodies[i];
        index++;
      }else{
        if(DEBUG) printf("Body at position (%e, %e, %e) is being skipped over by processor %d\n",
            bodies[i].position_x, bodies[i].position_y, bodies[i].position_z, rank);
      }
    }

    if(index != howManyWeHave - numLost){
      fprintf(stderr, "Processor %d copied over %d elements, but we originally had %d and lost %d\n", rank, index, howManyWeHave, numLost);
      exit(-9);
    }

    howManyWeHave = howManyWeHave - numLost;

    int numSending[8] = {0};

    MPI_Allgather(&numLost, 1, MPI_INT, numSending, 1, MPI_INT, MPI_COMM_WORLD);

    if(numSending[rank] != numLost){
      fprintf(stderr, "Process %d getting rid of %d bodies, but according to allgather we're getting rid of %d", rank, numLost, numSending[rank]);
    }

    // now send 'em over!
    
    int send_offset_helper = 0; // this ultimately becomes the total number leaving all quadrants!
    int send_offsets[8] = {0};

    for(i = 0; i < 8; i++){
      send_offsets[i] = send_offset_helper;
      send_offset_helper += numSending[i];
    }

    body *recvBuff = calloc(send_offset_helper, sizeof(body));
    MPI_Allgatherv(bodiesToSend, numLost, body_type, recvBuff, numSending, send_offsets, body_type, MPI_COMM_WORLD);

    for(i = 0; i < send_offset_helper; i++){
      if(distToQuad(recvBuff[i], rank) == 0.0){
        if(DEBUG) printf("Processor %d is adopting body with coordinates (%e, %e, %e)\n",
            rank, recvBuff[i].position_x, recvBuff[i].position_y, recvBuff[i].position_z);
        bodies[howManyWeHave] = recvBuff[i];
        howManyWeHave++;
      }
    }
  }
  postUpdate(rank, round);
}


void tests(){
  body b;
  b.position_x = 3.0;
  b.position_y = 3.0;
  b.position_z = 3.0;

  if(distToQuad(b, 0) != 0) fprintf(stderr, "Failed test! DistToQuad failed attempt 1\n");
  if(distToQuad(b, 1) != 3.0) fprintf(stderr, "Failed test! DistToQuad failed attempt 2\n");
  if(distToQuad(b, 2) != sqrt(18.0)) fprintf(stderr, "Failed test! DistToQuad failed attempt 3\n");
  if(distToQuad(b, 3) != 3.0) fprintf(stderr, "Failed test! DistToQuad failed attempt 4\n");
  if(distToQuad(b, 4) != 3.0) fprintf(stderr, "Failed test! DistToQuad failed attempt 5\n");
  if(distToQuad(b, 5) != sqrt(18.0)) fprintf(stderr, "Failed test! DistToQuad failed attempt 6\n");
  if(distToQuad(b, 6) != sqrt(27.0)) fprintf(stderr, "Failed test! DistToQuad failed attempt 7\n");
  if(distToQuad(b, 7) != sqrt(18.0)) fprintf(stderr, "Failed test! DistToQuad failed attempt 8\n");

  if(DEBUG) printf("All tests complete!\n");
}
