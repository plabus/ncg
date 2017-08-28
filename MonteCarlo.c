#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <inttypes.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"


// Initialise a random HERMITIAN matrix of type H or type L (that is tracefree or
// not), where the random elements lie in range [ -step_length_diag, step_length_diag ]
// on the diagonal and in [ -step_length_off*(1+i), spep_length_off*(1+i) ] for
// the off-diagonal elements
void random_matrix(
    struct pcg32_random_t* rng,    // pointer to random number generator
    REAL complex *M,               // pointer to first element of N x N matrix
    int const length,              // side length N
    double const step_length_diag, // range of diagonal elements
    double const step_length_off,  // range of off-diagonal elements
    enum Matrix_Type const type    // H or L? Non-traceless or traceless?
    )
{
  uint64_t const Nplus1 = length + 1;

  for( uint64_t i = 0; i < length; ++i )
  {
    uint64_t const itimesN = i * length;
    M[i*Nplus1] = step_length_diag * signed_uniform(rng) + I * 0.0 ; // M[i,i]

    for( uint64_t j = i + 1; j < length; ++j )
    {
      double const re = step_length_off * signed_uniform(rng);
      double const im = step_length_off * signed_uniform(rng);
      M[itimesN+j]  = re + I * im; // M[i,j]
      M[j*length+i] = re - I * im; // M[j,i]
    }
  }

  // If the matrix is of L_TYPE, the trace needs to
  // be zero. We calculate the sum of the first
  // (N-1) matrix elements and set the last diagonal
  // element to the negative of this sum.
  if( type == L_TYPE )
  {
    double complex trace = 0.0;
    for( uint64_t i = 0; i < length - 1; ++i )
    {
      trace += M[i*Nplus1]; // M[i,i]
    }
    M[length*length-1] = -trace;
  }
}

// Initialise (num_h + num_l) HERMITIAN matrices, where num_h are
// of H_TYPE and num_l are of L_TYPE. The range of the diagonal
// and off-diagonal elements can be varied of H_TYPE and L_TYPE
// independently
void random_matrices(
    struct pcg32_random_t *rng,      // array of rng's, one for each matrix
    REAL complex *Ms,                // pointer to first element of (num_h + num_l) N x N matrices
    int const num_h,                 // number of matrices of H_TYPE
    int const num_l,                 // number of matrices of L_TYPE
    int const length,                // side length N (the same for all matrices)
    double const step_length_diag_h, // range of diagonal elements for matrix H_TYPE
    double const step_length_off_h,  // range of off-diagonal elements for matrix H_TYPE
    double const step_length_diag_l, // range of diagonal elements for matrix H_TYPE
    double const step_length_off_l   // range of off-diagonal elements for matrix H_TYPE
    )
{
  for( uint64_t n = 0; n < num_h; ++n )
  {
    uint64_t const offset = n * length * length;
    random_matrix( &rng[n], &Ms[offset], length, step_length_diag_h, step_length_off_h, H_TYPE );
  }

  for( uint64_t n = num_h; n < num_h + num_l; ++n )
  {
    uint64_t const offset = n * length * length;
    random_matrix( &rng[n], &Ms[offset], length, step_length_diag_l, step_length_off_l, L_TYPE );
  }
}

// Initialise all Matrices for the MCMC, and calculate the initial action
void Matrices_Initialisation(
    struct pcg32_random_t *rng, // array of rng's, one for each matrix
    REAL complex *Matrices,     // array of matrices
    int const num_h,            // number of matrices of H_TYPE
    int const num_l,            // number of matrices of L_TYPE
    int const length,           // side length N (the same for all matrices)
    double* action              // pointer to the action variable
    )
{
  // Set high-temperature initial state for all  matrices, where random elements are in the range [-1-i, 1+i]
  random_matrices( rng, Matrices, num_h, num_l, length, 1.0, 1.0, 1.0, 1.0 );

  // Calculate the inital action
  *action = G2 * traceD2( Matrices, num_h, num_l, length ) + G4 * traceD4( Matrices, num_h, num_l, length );
}


// Creates a new Markov chain element
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng, // array of rng's, one for each matrix
    REAL complex *Matrices,     // array of matrices
    double *action,             // pointer to the value of the action of the matrices
    int *sigmaAB,               // precalculated Clifford products of 2 Gamma matrices
    int **sigmaABCD,            // precalculated Clifford products of 4 Gamma matrices
    const int NUM_H,            // number of matrices of H_TYPE
    const int NUM_L,            // number of matrices of L_TYPE
    uint64_t* acc,              // pointer to number of accepted steps
    const double step_size      // length of each Monte Carlo step
    )
{
  REAL complex temp;       // To save random value that is changed

  /* For each Matrix change a value in the upper triangle randomly *
   * calculate the the change of the action and finally decide if  *
   * the new matrix should be accepted.                            */
  for( int n = 0; n < NUM_M; ++n )
  {
    /* Set the offset to write to the right matrix */
    const int offset = n * SWEEP;

    /* Calculate random double in [0,1) for Monte Carlo Move Decision */
    const double p = uniform( &rng[n] );

    /* Calculate two random integers and generate position in upper and lower half */
    const int pos_x = uniform_int( &rng[n], N );
    const int pos_y = uniform_int( &rng[n], N );
    const int pos_upper = pos_x <= pos_y ? pos_x * N+pos_y : pos_y * N+pos_x;
    const int pos_lower = pos_x >  pos_y ? pos_x * N+pos_y : pos_y * N+pos_x;

    const REAL complex old = Matrices[pos_upper+offset];

    if(pos_x != pos_y) {
      temp  = step_size * signed_uniform( &rng[n] ) + I * step_size * signed_uniform( &rng[n] );
    } else {
      temp  = step_size * signed_uniform( &rng[n] );
    }

    if(pos_x != pos_y) {
      Matrices[pos_upper+offset] = temp;
      Matrices[pos_lower+offset] = conj(temp);
    } else {
      Matrices[pos_upper+offset] = temp;
    }

    double action_new = traceD2(Matrices, NUM_H, NUM_L, N);
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
