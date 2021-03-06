#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
using namespace std;

//
//  tuned constants
//
#define density 0.0005
#define mass 0.01
#define cutoff 0.01
#define min_r (cutoff / 100)
#define dt 0.0005

// calculate particle's bin number
int binNum(particle_t &p, int bpr)
{
  return (floor(p.x / cutoff) + bpr * floor(p.y / cutoff));
}

void clearBins(vector<particle_t *> *bins, int numbins)
{
  for (int m = 0; m < numbins; m++)
    bins[m].clear();
};

void placeParticlesInBins(int n, vector<particle_t *> *bins, particle_t *particles, int bpr)
{
  for (int i = 0; i < n; i++)
    bins[binNum(particles[i], bpr)].push_back(particles + i);
}

void computeForces(int n, vector<particle_t *> *bins, particle_t *particles, int bpr)
{
  for (int p = 0; p < n; p++)
  {
    particles[p].ax = particles[p].ay = 0;

    // find current particle's bin, handle boundaries
    int cbin = binNum(particles[p], bpr);
    int lowi = -1, highi = 1, lowj = -1, highj = 1;
    if (cbin < bpr)
      lowj = 0;
    if (cbin % bpr == 0)
      lowi = 0;
    if (cbin % bpr == (bpr - 1))
      highi = 0;
    if (cbin >= bpr * (bpr - 1))
      highj = 0;

    // apply nearby forces
    for (int i = lowi; i <= highi; i++)
      for (int j = lowj; j <= highj; j++)
      {
        int nbin = cbin + i + bpr * j;
        for (int k = 0; k < bins[nbin].size(); k++)
          apply_force(particles[p], *bins[nbin][k]);
      }
  }
}

void moveParticles(int n, particle_t *particles)
{
  for (int p = 0; p < n; p++)
    move(particles[p]);
}

//
//  benchmarking program
//
int main(int argc, char **argv)
{
  if (find_option(argc, argv, "-h") >= 0)
  {
    printf("Options:\n");
    printf("-h to see this help\n");
    printf("-n <int> to set the number of particles\n");
    printf("-o <filename> to specify the output file name\n");
    return 0;
  }

  int n = read_int(argc, argv, "-n", 50000);

  char *savename = read_string(argc, argv, "-o", NULL);

  FILE *fsave = savename ? fopen(savename, "w") : NULL;
  particle_t *particles = (particle_t *)malloc(n * sizeof(particle_t));
  set_size(n);
  init_particles(n, particles);

  // create spatial bins (of size cutoff by cutoff)
  double size = sqrt(density * n);
  int bpr = ceil(size / cutoff);
  int numbins = bpr * bpr;
  vector<particle_t *> *bins = new vector<particle_t *>[numbins];

  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer();
  for (int step = 0; step < NSTEPS; step++)
  {

    // clear bins at each time step
    clearBins(bins, numbins);

    // place particles in bins
    placeParticlesInBins(n, bins, particles, bpr);

    //
    //  compute forces
    //
    computeForces(n, bins, particles, bpr);

    //
    //  move particles
    //
    moveParticles(n, particles);

    //
    //  save if necessary
    //
    if (fsave && (step % SAVEFREQ) == 0)
      save(fsave, n, particles);
  }
  simulation_time = read_timer() - simulation_time;

  printf("n = %d, simulation time = %g seconds\n", n, simulation_time);

  free(particles);
  delete[] bins;
  if (fsave)
    fclose(fsave);

  return 0;
}
