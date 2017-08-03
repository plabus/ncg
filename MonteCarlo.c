#include <complex.h>
#include <math.h>
#include <inttypes.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"


// Initialise all Matrices, their Eigenvalues and the action
void Matrices_Initialisation(
    struct pcg32_random_t *rng,
    REAL complex *Matrices,
    double *action,
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
void Get_Next_MCMC_Element(struct pcg32_random_t *rng, REAL complex *Matrices, double *action,
                           int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L, int *acc, double step_size)
{
  int pos_x, pos_y;
  int pos_upper, pos_lower;
  double p;                 // Random double for accepting MC elemement
  REAL complex temp;       // To save random value that is changed
  double delta_action;

  /* For each Matrix change a value in the upper triangle randomly *
   * calculate the the change of the action and finally decide if  *
   * the new matrix should be accpted.                             */
  for(int n=0;n<NUM_M;++n)
  {
    /* Set the offset to write to the right matrix */
    int offset = n * SWEEP;

    /* Calculate random double in [0,1) for Monte Carlo Move Decision */
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
    if( delta_action<=0 || exp(-delta_action) > p ) {
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
    double acceptance_rate, // acceptance rate so far
    double* step_size      // reference to step_size to be tuned
    )
{
  // Assuming there is an ideal aceptance rate for Metropolis-Hastings,
  // we want to decrease the step size if the measured rate is smaller than the ideal one and
  // we want to increase the step size if the measured rate is larger  than the ideal one.
  // We assume this rate to be 23%
  *step_size *= acceptance_rate / 0.23;
}
