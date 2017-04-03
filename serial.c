#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VECZERO {0}

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
vector stepVelocity(vector a, vector b, vector forc, float DT){
  vector toReturn;

  toReturn.x = a.velocity.x + (DT * (forc.x / a.mass));
  toReturn.y = a.velocity.y + (DT * (forc.y / a.mass));
  toReturn.z = a.velocity.z + (DT * (forc.z / a.mass));

  return toReturn;
}

// get new coords
vector stepCoords(vector a, vector b, vector forc, float DT){
  vector avgVelocity, toReturn;

  avgVelocity.x = a.velocity.x + (0.5 * DT * (forc.x / a.mass));
  avgVelocity.y = a.velocity.y + (0.5 * DT * (forc.y / a.mass));
  avgVelocity.z = a.velocity.z + (0.5 * DT * (forc.z / a.mass));

  toReturn.x = a.x + avgVelocity.x * DT;
  toReturn.y = a.y + avgVelocity.y * DT;
  toReturn.z = a.z + avgVelocity.z * DT;

  return toReturn;
}

void step(int N, float DT, struct body *bodies){
  int i, j;
  float dist;
  vector forc, accel, veloc;

  vector *newCoords = calloc(N, sizeof(vector));
  vector *newVelocities = calloc(N, sizeof(vector));

  for(i = 0; i < N; i++){
    forc = VECZERO;

    // sum the force of all bodies on i
    for(j = 0; j < N; j++){
      if(i == j) continue; // skip myself

      dist = distance(bodies[i].position, bodies[j].position);

      if(dist < 5.0) continue;

      forc = vecsum(forc, force(bodies[i], bodies[j]));
    }

    newVelocities[i] = stepVelocity(a, b, forc);
    newCoords[i] = stepCoords(a, b, forc);
  }
}

int main (int argc, char **argv){
  int N, K, i;
  float DT;
  struct body *bodies;

  scanf("%d\n%d\n%f", &N, &K, &DT);
 
  bodies = calloc(N, sizeof(struct body));

  for(i = 0; i < N; i++){
    scanf("%f\n", &(bodes[i].mass))
  }
  for(i = 0; i < N; i++){
    scanf("%f %f %f\n", &(bodes[i].position.x), &(bodes[i].position.y), &(bodes[i].position.z));
  }
  for(i = 0; i < N; i++){
    scanf("%f %f %f\n", &(bodes[i].velocity.x), &(bodes[i].velocity.y), &(bodes[i].velocity.z));
  }

  return 0;
}
