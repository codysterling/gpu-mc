#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include "cody-mc-functions.c"
#include <omp.h>

struct timeval start,stop;

int main (int argc, char* argv[]) {
  //////////
  //// Defining variables
  gettimeofday(&start, NULL);
  int seed = 1; // seed for random number reproducability
//  srand(seed);

  // For MC run
  int natoms = atoi(argv[1]); // number of atoms in the box
  int nsteps = atoi(argv[2]); // total number of steps to take
  int wsteps = atoi(argv[3]); // step number to sample RDF
  printf("Running total %i steps (sample each %i) with %i atoms\n",nsteps,wsteps,natoms);
  float dist = 0.5; // max size of random move in Å
  int naccept = 0;
  int prodstep = nsteps/2; // step to start production (writing trajectory, RDFS)

  // Argon stuff, LJ and mass
  float sigma = 3.405; // in Å
  float eps = 1.654E-1; // in kgÅ^2/s^2
  float mass = 6.6335209E-26; // in kg
  float kb = 1.3806482E-3; // in Å^2*kg/s^2*K
  float t = 94.4; // in K
  // Box definitions
  float p = pow(natoms,1/3.);
  float xlen = pow(2,1/6.)*sigma*(floor(p)+ceil(p-floor(p))-1);
  float ylen = pow(2,1/6.)*sigma*(floor(p)+ceil(p-floor(p))-1);
  float zlen = pow(2,1/6.)*sigma*(floor(p)+ceil(p-floor(p))-1);

  // RDF variables
  float dr = 0.1; // bin width in Å
  float rdist = tmin(3,xlen/2,ylen/2,zlen/2); // in Å, radius of inscribed circle
  int nblks = rdist/dr;
  float g[nblks];
  for (int i=0;i<sizeof(g)/sizeof(float);i++) {g[i] = 0;}

  //////////
  //// Starting program
  // Initializing atom position arrays
  float atompos[natoms][3];
  float natompos[natoms][3];

  // Initializing random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_env_setup();
  gsl_rng_set(r,seed);

  // Initializing atom positions into the atompos array
  initpos(atompos,natoms,xlen,ylen,zlen,r);

  // Now doing each MC step
  // Production (with printing):
  #pragma omp parallel firstprivate(atompos,natompos,naccept) shared(g)
  {
  // Getting OMP variables
  int nthreads = omp_get_num_threads();
  int pid = omp_get_thread_num();

  // Initializing the random number generator for each thread
  gsl_rng *tr = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(tr,pid*nthreads);

  // Getting total number of steps
  int tsteps = prodstep + (nsteps-prodstep)/nthreads;

  printf("Thread %i doing %i total steps\n",pid,tsteps);
  printf("Thread %i doing %i equilibration steps\n",pid,prodstep);
  for (int i=0;i<tsteps;i++) {
    // Sampling RDF only if in production part
    if (i >= prodstep && i%wsteps == 0) {
			rdf(atompos,natoms,xlen,ylen,zlen,dr,rdist,g);
    }

    // Choosing random atom to move
    int randn = gsl_rng_uniform(tr)*natoms;
    // Getting initial energy of that atom
    float ben = epot(atompos,natoms,randn,sigma,eps,xlen,ylen,zlen);

    // Copying positions to new array and nudging selected atom
    memcpy(natompos,atompos,sizeof(natompos));
    nudge(natompos,natoms,randn,dist,tr);
    // Getting new energy of that atom
    float aen = epot(natompos,natoms,randn,sigma,eps,xlen,ylen,zlen);

    // Checking step acceptance, updating positions if needed
    if (gsl_rng_uniform(tr) < fmin(exp(-(aen-ben)/(kb*t)),1.)) {
      memcpy(atompos,natompos,sizeof(atompos));
      naccept += 1;
    }
    if (i == prodstep-1) {printf("naccept after equilibration: %i, ratio: %f\n",naccept,(float)naccept/(float)prodstep);ptime(start,"after eq:");printf("Thread %i starting %i production steps\n",pid,tsteps-prodstep);}
  }
  printf("naccept after production: %i, ratio: %f\n",naccept,(float)naccept/(float)nsteps);
  }

  // Writing RDF
  char fname[45];
  snprintf(fname,sizeof(fname),"serial-rdf-n%i-s%i-w%i.rdf",natoms,nsteps,wsteps);
  FILE *rdfile = fopen(fname,"w");
  for (int i=0;i<nblks;i++) {
    fprintf(rdfile,"%f %f\n",i*dr,g[i]);
  }
  fclose(rdfile);

  // Timing
  gettimeofday(&stop, NULL);
  double secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
  printf("Total elapsed time: %f s\n",secs);
}
