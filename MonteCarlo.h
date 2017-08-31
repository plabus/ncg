#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"
#include "Constants.h"
#include "Matrix_Properties.h"

// Types of HERMITIAN matrices
// H_TYPE has no constraints on the trace
// L_TYPE is trace-free
enum Matrix_Type { H_TYPE, L_TYPE };

// This is used to save a state of matrices
// when using an ultra-local update
struct Matrix_State
{
  size_t matrix;
  size_t pos_upper;
  size_t pos_lower;
  REAL complex matrix_element;
};

// Return random signed double in range (-1, 1) uniformly distributed
double signed_uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    );

// Initialise a random HERMITIAN matrix of type H or type L (that is tracefree or
// not), where the random elements lie in range [ -step_length_diag, step_length_diag ]
// on the diagonal and in [ -step_length_off*(1+i), step_length_off*(1+i) ] for
// the off-diagonal elements
void random_matrix(
    struct pcg32_random_t* rng,    // pointer to random number generator
    REAL complex *M,               // pointer to first element of N x N matrix
    int const length,              // side length N
    double const step_length_diag, // range of diagonal elements
    double const step_length_off,  // range of off-diagonal elements
    enum Matrix_Type const type    // H or L? Non-traceless or traceless?
    );

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
    );

// Initialise all Matrices for the MCMC, and return the initial action
double Matrices_Initialisation(
    struct pcg32_random_t *rng,         // array of rng's, one for each matrix
    REAL complex *Matrices,             // array of matrices
    struct Matrix_Properties const prop // includes num_h, num_l, n and k
    );

// Generate a new Monte Carlo candidate by changing one matrix element in one matrix
struct Matrix_State Generate_Candidate(
    struct pcg32_random_t *rng,                // pointer to one rng
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    const double step_size,                    // length of each Monte Carlo step
    const int matrix                           // in which matrix an element should be changed
    );

// Restore the Matrices as they have been before using Generate_Candidate
void Restore_Matrices(
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    const struct Matrix_State old              // old state
    );

// For each of the NUM_H + NUM_L matrices generate a new matrix
// with exactly one entry changed and make a accept/reject step
// individually. The acceptance counter is update for each matrix
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng,                // array of rng's, one for each matrix
    REAL complex *Matrices,                    // array of matrices
    struct Matrix_Properties const parameters, // includes num_h, num_l, n and k
    int *sigmaAB,                              // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD,                           // pre-calculated Clifford products of 4 Gamma matrices
    uint64_t* accepted,                        // pointer to number of accepted steps
    const double step_size                     // length of each Monte Carlo step
    );

// Tune the step size based on some acceptance rate measured
void tune_step_size(
    double accepance_rate, // acceptance rate so far
    double* step_size      // reference to step_size to be tuned
    );

#endif
