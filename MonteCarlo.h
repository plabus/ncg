#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"
#include "Constants.h"

// Return random signed double in range (-1, 1) uniformly distributed
double signed_uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    );

// Initialise each of the NUM_H + NUM_L matrices
// randomly with elements in the range [-1-i, 1+i]
void Matrices_Initialisation(
    struct pcg32_random_t *rng, // array of rng's, one for each matrix
    REAL complex *Matrices      // array of matrices
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
