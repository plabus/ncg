#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"
#include "Constants.h"

// Types of HERMITIAN matrices
// H_TYPE has no constraints on the trace
// L_TYPE is trace-free
enum Matrix_Type { H_TYPE, L_TYPE };


// Return random signed double in range (-1, 1) uniformly distributed
double signed_uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    );

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
    );

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
    );

// Initialise all Matrices for the MCMC, and calculate the initial action
void Matrices_Initialisation(
    struct pcg32_random_t *rng, // array of rng's, one for each matrix
    REAL complex *Matrices,     // array of matrices
    int const num_h,            // number of matrices of H_TYPE
    int const num_l,            // number of matrices of L_TYPE
    int const length,           // side length N (the same for all matrices)
    double* action              // pointer to the action variable
    );

// For each of the NUM_H + NUM_L matrices generate a new matrix
// with exactly one entry changed and make a accept/reject step
// individually. The acceptance counter is update for each matrix
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng,
    REAL complex *Matrices,
    double *action,
    int *sigmaAB,
    int **sigmaABCD,
    int NUM_H,
    int NUM_L,
    uint64_t *acc,
    double step_size
    );

// Tune the step size based on some acceptance rate measured
void tune_step_size(
    double accepance_rate, // acceptance rate so far
    double* step_size      // reference to step_size to be tuned
    );

#endif
