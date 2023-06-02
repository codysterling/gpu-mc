#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

//////////
// Defining functions
// new power function since exponents are always integers, should be faster
float npow (float base, int power) {
  float ans = base;
  for (int i=1;i<power;i++) {
    ans *= base;
  }
  return ans;
}

void ptime (struct timeval start, char* msg) {
  struct timeval ct;
  gettimeofday(&ct,NULL);
  float secs = (double)(ct.tv_usec - start.tv_usec) / 1000000 + (double)(ct.tv_sec - start.tv_sec);
  printf("%s: %f s\n",msg,secs);
}

// min function retooled from https://www.geeksforgeeks.org/variable-length-argument-c/
float tmin (int n, ...) {
  double min, a;

  va_list ap;
  va_start(ap,n);
  min = va_arg(ap,double);

  for (int i=2;i<=n;i++) {
    if ((a = va_arg(ap,double)) < min) {
      min = a;
    }
  }

  va_end(ap);
  return min;
}

void amult (float* arr, int size, float arg) {
  for (int i=0;i<size;i++) {
    arr[i] = arr[i]*arg;
  }
}

float vlen (float* arr, int size) {
  float len = 0.0;

  for (int i=0;i<size;i++) {
    len += npow(arr[i],2);
  }
  len = sqrt(len);
  return len;
}

float ljpot (float r2, float sigma, float eps) {
  float ljp;

  ljp = 4*eps*(npow(sigma/r2,12) - npow(sigma/r2,6));
  return ljp;
}

void minimg (float r[3], float xlen, float ylen, float zlen) {
  r[0] = r[0] - round(r[0]/xlen)*xlen;
  r[1] = r[1] - round(r[1]/ylen)*ylen;
  r[2] = r[2] - round(r[2]/zlen)*zlen;
}

void initpos (float atompos[][3], int natoms, float xlen, float ylen, float zlen, gsl_rng *r) {
  float p = pow(natoms,1/3.);
  int nside = floor(p)+ceil(p-floor(p));

  for (int i=0;i<nside;i++) {
    for (int j=0;j<nside;j++) {
      for (int k=0;k<nside;k++) {
        int sum = i*npow(nside,2) + j*nside + k;
        // Checking if we've reached the total number of atoms to place
        if (sum >= natoms) {
          return;
        }
        // Updating atompos array
        atompos[sum][0] = i*xlen/(nside-1)+(gsl_rng_uniform(r)-0.5);
        atompos[sum][1] = j*ylen/(nside-1)+(gsl_rng_uniform(r)-0.5);
        atompos[sum][2] = k*zlen/(nside-1)+(gsl_rng_uniform(r)-0.5);
      }
    }
  }
}

float epot (float atompos[][3], int natoms, int n, float sigma, float eps, float xlen, float ylen, float zlen) {
  float en = 0.0;
  float npos[3] = {atompos[n][0], atompos[n][1], atompos[n][2]};

  for (int i=0;i<natoms;i++) {
    if (i != n) {
      float ipos[3] = {atompos[i][0], atompos[i][1], atompos[i][2]};
      float r[3] = {npos[0]-ipos[0],npos[1]-ipos[1],npos[2]-ipos[2]};
      minimg(r,xlen,ylen,zlen);
      float r2 = vlen(r,3);
      en += ljpot(r2,sigma,eps);
    }
  }

  return en;
}

void nudge (float atompos[][3], int natoms, int n, float dist, gsl_rng *r) {
  atompos[n][0] += (gsl_rng_uniform(r) - 0.5)*2*dist;
  atompos[n][1] += (gsl_rng_uniform(r) - 0.5)*2*dist;
  atompos[n][2] += (gsl_rng_uniform(r) - 0.5)*2*dist;
}

void wrappositions(float atompos[][3], int natoms, float xlen, float ylen, float zlen) {
  for (int i=0;i<natoms;i++) {
    atompos[i][0] -= floor(atompos[i][0]/xlen)*xlen;
    atompos[i][1] -= floor(atompos[i][1]/ylen)*ylen;
    atompos[i][2] -= floor(atompos[i][2]/zlen)*zlen;
  }
}

void rdf(float atompos[][3], int natoms, float xlen, float ylen, float zlen, float dr, float dist, float g[]) {
  int nbins = dist/dr;

  for (int i=0;i<natoms-1;i++) {
    float ipos[3] = {atompos[i][0], atompos[i][1], atompos[i][2]};
    for (int j=i+1;j<natoms;j++) {
      float jpos[3] = {atompos[j][0], atompos[j][1], atompos[j][2]};
      float r[3] = {ipos[0]-jpos[0],ipos[1]-jpos[1],ipos[2]-jpos[2]};
      minimg(r,xlen,ylen,zlen);
      float r2 = vlen(r,3);

      if (r2 <= dist) {
        int ig = r2/dr;
        g[ig] += 2;
      }
    }
  }

  for (int i=0;i<nbins;i++) {
    float vb = ((npow(i+1,3) - npow(i,3))*npow(dr,3));
    float nid = (4*M_PI*vb*natoms)/(3*xlen*ylen*zlen);
    g[i] /= natoms*nid;
  }
}
