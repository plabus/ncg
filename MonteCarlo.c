#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <inttypes.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"

// Return random signed double in range (-1, 1) uniformly distributed
inline double signed_uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    )
{
  return ( pcg32_boundedrand_r(rng, 2) ? -1 : 1 ) * ldexp( pcg32_random_r(rng), -32 );
}

// Initialise all Matrices, their Eigenvalues and the action
void Matrices_Initialisation(
    struct pcg32_random_t *rng, // array of rng's, one for each matrix
    REAL complex *Matrices      // array of matrices
    )
{
  // Set high-temperature initial state for all
  // matrices, where random elements are in the
  // range [-1-i, 1+i]

  uint64_t Nplus1 = N + 1;

  for( uint64_t n = 0; n < NUM_M; ++n )
  {
    uint64_t offset = n*SWEEP;

    for( uint64_t i = 0; i < N; ++i )
    {
      uint64_t itimesN = i * N;
      Matrices[i*Nplus1+offset] = signed_uniform( &rng[n] ) + I * 0.0 ;

      for( uint64_t j = i + 1; j < N; ++j )
      {
        double re = signed_uniform( &rng[n] );
        double im = signed_uniform( &rng[n] );

        Matrices[itimesN+j+offset] = re + I * im;
        Matrices[j*N+i+offset]     = re - I * im;
      }

    }

  }

}


// Creates a new Markov chain element
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng,
    REAL complex *Matrices,
    double *action,
    int *sigmaAB,
    int **sigmaABCD,
    int NUM_H,
    int NUM_L,
    uint64_t* acc,
    double step_size
    )
{
  int pos_x, pos_y;
  int pos_upper, pos_lower;
  double p;                 // Random double for accepting MC elemement
  REAL complex temp;       // To save random value that is changed
  REAL complex old;

  /* For each Matrix change a value in the upper triangle randomly *
   * calculate the the change of the action and finally decide if  *
   * the new matrix should be accepted.                            */
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

    old = Matrices[pos_upper+offset];

    if(pos_x != pos_y) {
      temp  = step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32)
        + I * step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
    } else {
      temp  = step_size * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32) + 0.0 * I;
    }

    if(pos_x != pos_y) {
      Matrices[pos_upper+offset] = temp;
      Matrices[pos_lower+offset] = conj(temp);
    } else {
      Matrices[pos_upper+offset] = temp;
    }

    double action_new = traceD2(Matrices, NUM_H, NUM_L);
    double delta_action = action_new - (*action);

    /* Finally test if new action is smaller or except randomly if exp supressed, *
     * if yes write new element in upper and lower half and copy new eigenvalues  */
    // fprintf( stdout, "  deltaS = %f, p = %f, ", delta_action, p );
    if( delta_action<=0.0 || exp(-delta_action) > p )
    {
      // accepted
      *acc += 1;
      *action = (*action) + delta_action;
    }
    else
    {
      // rejected
      if(pos_x != pos_y) {
        Matrices[pos_upper+offset] = old;
        Matrices[pos_lower+offset] = conj(old);
      } else {
        Matrices[pos_upper+offset] = old;
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
