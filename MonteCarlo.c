#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <inttypes.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"

// Initialise a random HERMITIAN matrix of type H or type L (that is tracefree or not),
// where the random elements lie in range [ -step_length_diag, step_length_diag ]
// on the diagonal and in [ -step_length_off*(1+i), step_length_off*(1+i) ] for
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
    struct pcg32_random_t *rng,          // array of rng's, one for each matrix
    REAL complex *Ms,                    // pointer to first element of (num_h + num_l) N x N matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    double const step_length_diag_h,     // range of diagonal elements for matrix H_TYPE
    double const step_length_off_h,      // range of off-diagonal elements for matrix H_TYPE
    double const step_length_diag_l,     // range of diagonal elements for matrix H_TYPE
    double const step_length_off_l       // range of off-diagonal elements for matrix H_TYPE
    )
{
  // Unpack property parameters
  size_t const num_h = prop.num_h;
  size_t const num_l = prop.num_l;
  size_t const n     = prop.n;

  // All matrices of H_TYPE
  for( uint64_t m = 0; m < num_h; ++m )
  {
    uint64_t const offset = m * n * n;
    random_matrix( &rng[m], &Ms[offset], n, step_length_diag_h, step_length_off_h, H_TYPE );
  }

  // All matrices of L_TYPE
  for( uint64_t m = num_h; m < num_h + num_l; ++m )
  {
    uint64_t const offset = m * n * n;
    random_matrix( &rng[m], &Ms[offset], n, step_length_diag_l, step_length_off_l, L_TYPE );
  }
}

// Initialise all Matrices for the MCMC, and return the initial action
double Matrices_Initialisation(
    struct pcg32_random_t *rng,          // array of rng's, one for each matrix
    REAL complex *Matrices,              // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    )
{
  // Set high-temperature initial state for all  matrices, where random elements are in the range [-1-i, 1+i]
  random_matrices( rng, Matrices, prop, 1.0, 1.0, 1.0, 1.0 );

  // Calculate the initial action and return
  return Calculate_Action( Matrices, prop, sigmaAB, sigmaABCD );
}

// Generate a new Monte Carlo candidate by changing one matrix element in one matrix
struct Matrix_State Generate_Candidate(
    struct pcg32_random_t *rng,                // pointer to one rng
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    const double step_size,                    // length of each Monte Carlo step
    const int matrix                           // in which matrix an element should be changed
    )
{
  // Generate position of change and save old element:
  // -------------------------------------------------

  // Set the offset to write to the right matrix
  const int length = parameters.n;
  const int offset = matrix * length * length;

  // Calculate two random integers and generate position in upper and lower half
  const int pos_x = uniform_int( rng, length );
  const int pos_y = uniform_int( rng, length );
  const int pos_upper = pos_x <= pos_y ? pos_x * length + pos_y : pos_y * length + pos_x;
  const int pos_lower = pos_x >  pos_y ? pos_x * length + pos_y : pos_y * length + pos_x;

  // Save the old value
  const REAL complex old_element = Matrices[ offset + pos_upper ];
  const struct Matrix_State old_state = {
    .matrix         = matrix,
    .pos_x          = pos_x,
    .pos_y          = pos_y,
    .pos_upper      = pos_upper,
    .pos_lower      = pos_lower,
    .matrix_element = old_element
  };

  // Generate new matrix element:
  // ----------------------------

  REAL complex temp;

  // Off-diagonal case
  if(pos_upper != pos_lower)
  {
    temp = step_size * signed_uniform( rng ) + I * step_size * signed_uniform( rng );
    Matrices[ offset + pos_upper ] += temp;
    Matrices[ offset + pos_lower ] += conj(temp);
  }
  // Diagonal case
  else
  {
    temp = step_size * signed_uniform( rng ) + I * 0.0;
    Matrices[ offset + pos_upper ] += temp;

    // FIXME: Ensure traceless-ness, but be careful when restoring result!
    // // For L_TYPE ensure that the trace remains zero
    // enum Matrix_Type const type = matrix < parameters.num_h ? H_TYPE : L_TYPE;
    // if( type == L_TYPE )
    // {
    //   const int pos_diag = uniform_int( rng, length );
    //   Matrices[ offset + pos_diag * length + pos_diag ] -= temp;
    // }
  }

  return old_state;
}

// Restore the Matrices as they have been before using Generate_Candidate
void Restore_Matrices(
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    const struct Matrix_State old              // old state
    )
{
  // Set the offset to write to the right matrix
  const int length = parameters.n;
  const int offset = old.matrix * length * length;

  // Off-diagonal case
  if(old.pos_upper != old.pos_lower)
  {
    Matrices[ old.pos_upper + offset ] = old.matrix_element;
    Matrices[ old.pos_lower + offset ] = conj(old.matrix_element);
  }
  // Diagonal case
  else
  {
    Matrices[ old.pos_upper + offset ] = old.matrix_element;
  }
}

// Creates a new Markov chain element
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng,                // array of rng's, one for each matrix
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    int *sigmaAB,                              // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD,                           // pre-calculated Clifford products of 4 Gamma matrices
    uint64_t* accepted,                        // pointer to number of accepted steps
    const double step_size                     // length of each Monte Carlo step
    )
{
  // For each matrix, separately generate a new candidate,
  // calculate the resulting action and perform a Monte Carlo step
  for( size_t number_matrix = 0; number_matrix < parameters.num_h + parameters.num_l; ++number_matrix )
  {

    // Generate new candidate and calculate new action:
    // ------------------------------------------------

    const double action_old = Calculate_Action( Matrices, parameters, sigmaAB, sigmaABCD );
    const struct Matrix_State old =
      Generate_Candidate( &rng[number_matrix], Matrices, parameters, step_size, number_matrix );
    const double action_new = Calculate_Action( Matrices, parameters, sigmaAB, sigmaABCD );
    const double delta_action = action_new - action_old;

    // const double delta_action_naive = action_new - action_old;
    // const double delta_action = Calculate_Delta_Action( Matrices, parameters, sigmaAB, sigmaABCD, old );
    // if( fabs(delta_action_naive-delta_action) > 1e-08 )
    // {
    //   printf("  !!!!!!!!!!!!!!!!!!!!!! deltadeltaS = %g\n", delta_action_naive-delta_action);
    // }


    // Monte Carlo move decision:
    // --------------------------

    const double p = uniform( &rng[number_matrix] );

    // accept
    if( delta_action <= 0.0 || exp(-delta_action) > p )
    {
      *accepted += 1;
    }
    // reject
    else
    {
      Restore_Matrices( Matrices, parameters, old );
    }

  }
}


void tune_step_size(
    double acceptance_rate, // acceptance rate so far
    double* step_size       // reference to step_size to be tuned
    )
{
  // Assuming there is an ideal acceptance rate for Metropolis-Hastings,
  // we want to decrease the step size if the measured rate is smaller than the ideal one and
  // we want to increase the step size if the measured rate is larger  than the ideal one.
  // We assume this rate to be 23%

  // If the acceptance rate was zero set step size to small value:
  if( acceptance_rate < 1e-6 )
  {
    *step_size *= 0.0001;
  }
  else
  {
    *step_size *= acceptance_rate / 0.23;
  }
}
