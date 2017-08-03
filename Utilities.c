#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/time.h>
#include "Constants.h"
#include "Utilities.h"

// Returns time in seconds (double)
double cclock() {
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

// Print the Matrix (for debugging)
void printM(const int L, REAL complex * M) {
  for(int i = 0; i < L; i++) {
    for(int j = 0; j < L; j++)
      if( cimag(M[i*L+j])>=0) printf( " %.2f + %.2fi\t ", creal( M[i*L+j] ),  cimag( M[i*L+j] ) );
      else                    printf( " %.2f - %.2fi\t ", creal( M[i*L+j] ), -cimag( M[i*L+j] ) );
    printf("\n");
  }
}

// To compare two floats for sorting
int compare_floats(const void *a, const void *b) {
  const float *da = (const float *) a;
  const float *db = (const float *) b;
  return (*da > *db) - (*da < *db);
}

/* Delta approximation for EV distribtion */
float delta_approx(float x) {
  float abs = x > 0 ? 1.f-x/EPSILON : 1.f+x/EPSILON;
  return ( abs > 0 ? abs/EPSILON : 0 );
}

// Erase the content of a complex square matrix with side length l
void nullify(const int l, float complex * m) {
  for(int i=0;i<l*l;++i) {
    m[i] = 0.0 + 0.0 * I;
  }
}
