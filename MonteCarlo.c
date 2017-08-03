#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/time.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"



// Returns time in seconds (double)
double cclock() {
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

// Print the Matrix (for debugging)
void printM(const int L, float complex * M) {
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

// Initialise all Matrices, their Eigenvalues and the action
void Matrices_Initialisation(
    struct pcg32_random_t *rng,
    float complex *Matrices,
    float *action,
    int NUM_H,
    int NUM_L
    )
{
  // Set high-temperature initial state for all
  // matrices, where random elements are in the
  // range [-1-i, 1+i], and calculate initial action

  int itimesN;
  int Nplus1 = N+1;
  int offset;

  for( int n = 0; n < NUM_M; ++n )
  {
    offset = n*SWEEP;
    for( int i = 0; i < N; ++i )
    {
      Matrices[i*Nplus1+offset] = (pcg32_boundedrand_r(&rng[n],2)?-1:1) * ldexp(pcg32_random_r(&rng[n]),-32)
                            + I * 0.0 ;
      itimesN = i*N;
      for( int j = i + 1; j < N; ++j )
      {
        Matrices[itimesN+j+offset] = (pcg32_boundedrand_r(&rng[n],2)?-1:1) * ldexp(pcg32_random_r(&rng[n]),-32)
                               + I * (pcg32_boundedrand_r(&rng[n],2)?-1:1) * ldexp(pcg32_random_r(&rng[n]),-32);
        Matrices[j*N+i+offset]     =  conj( Matrices[itimesN+j+offset] ); /* This is hermitian! */
      }
    }
  }

  //*action = G2 * traceD2(Matrices, NUM_H, NUM_L) + G4 * traceD4(Matrices, NUM_H, NUM_L);
}


// Creates a new Markov chain element
void Get_Next_MCMC_Element(struct pcg32_random_t *rng, float complex *Matrices, float *action,
                           int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L, int *acc, float step_size)
{
  int pos_x, pos_y;
  int pos_upper, pos_lower;
  float p;                    // Random float for accepting MC elemement
  float complex temp;         // To save random value that is changed
  float delta_action;

  /* For each Matrix change a value in the upper triangle randomly *
   * calculate the the change of the action and finally decide if  *
   * the new matrix should be accpted.                             */
  for(int n=0;n<NUM_M;++n)
  {
    /* Set the offset to write to the right matrix */
    int offset = n * SWEEP;

    /* Calculate random float in [0,1) for Monte Carlo Move Decision */
    p = ldexp(pcg32_random_r(&rng[n]),-32);

    /* Calculate two random integers and generate position in upper and lower half */
    pos_x = pcg32_boundedrand_r(&rng[n],N);
    pos_y = pcg32_boundedrand_r(&rng[n],N);
    pos_upper = pos_x <= pos_y ? pos_x * N+pos_y : pos_y * N+pos_x;
    pos_lower = pos_x >  pos_y ? pos_x * N+pos_y : pos_y * N+pos_x;

    if(pos_x != pos_y) {
      temp  = step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32)
        + I * step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
    } else {
      temp  = step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32) + 0.0 * I;
    }
    temp += Matrices[pos_upper+offset];

    delta_action  = G2 * delta_action_traceD2(Matrices, n, temp, pos_x, pos_y, NUM_H, NUM_L);
    delta_action += G4 * delta_action_traceD4(Matrices, n, temp, pos_x, pos_y, sigmaAB, sigmaABCD, NUM_H, NUM_L);

    /* Finally test if new action is smaller or except randomly if exp supressed, *
     * if yes write new element in upper and lower half and copy new eigenvalues  */
    if( delta_action<=0 || expf(-delta_action) > p ) {
      *acc += 1;
      *action += delta_action;
      if(pos_x != pos_y) {
        Matrices[pos_upper+offset] = temp;
        Matrices[pos_lower+offset] = conj(temp);
      } else {
        Matrices[pos_upper+offset] = temp;
      }
    }
  }

}


void tune_step_size(
    float acceptance_rate, // acceptance rate so far
    float* step_size      // reference to step_size to be tuned
    )
{
  // Assuming there is an ideal aceptance rate for Metropolis-Hastings,
  // we want to decrease the step size if the measured rate is smaller than the ideal one and
  // we want to increase the step size if the measured rate is larger  than the ideal one.
  // We assume this rate to be 23%
  *step_size *= acceptance_rate / 0.23;
}
