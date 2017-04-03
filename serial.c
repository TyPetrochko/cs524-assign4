#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timing.h"

typedef struct vector {
  float x, y, z;
} vector;

struct body {
  float mass;
  vector position;
  vector velocity;
};

float distance(vector a, vector b){
  float dx, dy, dz;

  dx = b.x - a.x;
  dy = b.y - a.y;
  dz = b.z - a.z;

  return sqrt(dx*dx + dy*dy + dz*dz);
}

// force that acts on a (not b)
vector force(struct body a, struct body b){
  vector toReturn;

  float dist = distance(a.position, b.position);

  toReturn.x = (a.mass * b.mass * (b.position.x - a.position.x))/(dist * dist * dist);
  toReturn.y = (a.mass * b.mass * (b.position.y - a.position.y))/(dist * dist * dist);
  toReturn.z = (a.mass * b.mass * (b.position.z - a.position.z))/(dist * dist * dist);

  return toReturn;
}

// sum two vectors
vector vecsum(vector a, vector b){
  vector toReturn;
  toReturn.x = a.x + b.x;
  toReturn.y = a.y + b.y;
  toReturn.z = a.z + b.z;
  return toReturn;
}

// get new velocity
vector stepVelocity(struct body a, vector forc, float DT){
  vector toReturn;

  toReturn.x = a.velocity.x + (DT * (forc.x / a.mass));
  toReturn.y = a.velocity.y + (DT * (forc.y / a.mass));
  toReturn.z = a.velocity.z + (DT * (forc.z / a.mass));

  return toReturn;
}

// get new coords
vector stepCoords(struct body a, vector forc, float DT){
  vector avgVelocity, toReturn;

  avgVelocity.x = a.velocity.x + (0.5 * DT * (forc.x / a.mass));
  avgVelocity.y = a.velocity.y + (0.5 * DT * (forc.y / a.mass));
  avgVelocity.z = a.velocity.z + (0.5 * DT * (forc.z / a.mass));

  toReturn.x = a.position.x + avgVelocity.x * DT;
  toReturn.y = a.position.y + avgVelocity.y * DT;
  toReturn.z = a.position.z + avgVelocity.z * DT;

  return toReturn;
}

void step(int N, float DT, double wctime, struct body *bodies){
  int i, j;
  float dist;
  vector forc, accel, veloc;

  vector *newCoords = calloc(N, sizeof(vector));
  vector *newVelocities = calloc(N, sizeof(vector));

  for(i = 0; i < N; i++){
    forc.x = 0; forc.y = 0; forc.z = 0;
    // sum the force of all bodies on i
    for(j = 0; j < N; j++){
      if(i == j) continue; // skip myself

      dist = distance(bodies[i].position, bodies[j].position);

      if(dist < 5.0) continue;

      forc = vecsum(forc, force(bodies[i], bodies[j]));
    }
    newVelocities[i] = stepVelocity(bodies[i], forc, DT);
    newCoords[i] = stepCoords(bodies[i], forc, DT);
  }

  for(i = 0; i < N; i++){
    bodies[i].position = newCoords[i];
    bodies[i].velocity = newVelocities[i];
  }
}

vector calcCenterOfMass(struct body *bodies, int N){
  vector toReturn;
  double numerator, denominator;
  int i;

  toReturn.x = 0.0;
  toReturn.y = 0.0;
  toReturn.z = 0.0;

  numerator = 0.0; denominator = 0.0;
  for(i = 0; i < N; i++){
    numerator += (bodies[i].mass * bodies[i].position.x);
    denominator += bodies[i].mass;
    printf("X numerator %lf, denominator %lf\n", numerator, denominator);
  }
  toReturn.x = numerator / denominator;
  
  numerator = 0; denominator = 0;
  for(i = 0; i < N; i++){
    numerator += (bodies[i].mass * bodies[i].position.y);
    denominator += bodies[i].mass;
  }
  toReturn.y = numerator / denominator;
  
  numerator = 0; denominator = 0;
  for(i = 0; i < N; i++){
    numerator += (bodies[i].mass * bodies[i].position.z);
    denominator += bodies[i].mass;
  }
  toReturn.z = numerator / denominator;
  return toReturn;
}

vector calcAvgVelocity(struct body *bodies, int N){
  vector toReturn;
  int i;
  
  toReturn.x = 0;
  toReturn.y = 0;
  toReturn.z = 0;

  for(i = 0; i < N; i++){
    toReturn.x += bodies[i].velocity.x;
    toReturn.y += bodies[i].velocity.y;
    toReturn.z += bodies[i].velocity.z;
  }

  toReturn.x = toReturn.x / (float) N;
  toReturn.y = toReturn.y / (float) N;
  toReturn.z = toReturn.z / (float) N;

  return toReturn;
}

void postUpdate(int step, int N, double wctime, float DT, struct body *bodies){
  vector centerOfMass, avgVelocity;
  double endtime, cputime;

  timing(&endtime, &cputime);

  centerOfMass = calcCenterOfMass(bodies, N);
  avgVelocity = calcAvgVelocity(bodies, N); 

  printf("Conditions after timestep %d (time = %f):\n\n", step, step * DT);
  printf("\tCenter of Mass: (%f, %f, %f)\n", centerOfMass.x, centerOfMass.y, centerOfMass.z);
  printf("\tAverage Velocity: (%f, %f, %f)\n\n\n", avgVelocity.x, avgVelocity.y, avgVelocity.z);
}

int main (int argc, char **argv){
  int N, K, i;
  float DT;
  struct body *bodies;

  double wctime, cputime;

  scanf("%d\n%d\n%f", &N, &K, &DT);
 
  bodies = calloc(N, sizeof(struct body));

  for(i = 0; i < N; i++){
    scanf("%f\n", &(bodies[i].mass));
  }
  for(i = 0; i < N; i++){
    scanf("%f %f %f\n", &(bodies[i].position.x), &(bodies[i].position.y), &(bodies[i].position.z));
  }
  for(i = 0; i < N; i++){
    scanf("%f %f %f\n", &(bodies[i].velocity.x), &(bodies[i].velocity.y), &(bodies[i].velocity.z));
  }

  timing(&wctime, &cputime);
  for(i = 0; i < K; i++){
    if(i % 128 == 0){
      postUpdate(i, N, wctime, DT, bodies);
    }
    step(N, DT, wctime, bodies);
  }
  postUpdate(i, N, wctime, DT, bodies);
  return 0;
}

