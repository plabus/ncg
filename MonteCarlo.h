#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"

double cclock();
void printM(const int L, float complex * M);
int compare_floats(const void * a, const void * b);
float delta_approx(float);
void nullify(const int l, float complex * m);


// Initialise each of the NUM_H + NUM_L matrices
// randomly with elements in the range [-1-i, 1+i],
// and calculate (and set) initialial action
void Matrices_Initialisation(
    struct pcg32_random_t *rng,
    float complex *Matrices,
    float *action,
    int NUM_H,
    int NUM_L
    );

// For each of the NUM_H + NUM_L matrices generate a new matrix
// with exactly one entry changed and make a accept/reject step
// individually. The acceptance counter is update for each matrix
void Get_Next_MCMC_Element(
    struct pcg32_random_t *rng,
    float complex *Matrices,
    float *action,
    int *sigmaAB,
    int **sigmaABCD,
    int NUM_H,
    int NUM_L,
    int *acc,
    float step_size
    );

// Tune the step size based on some acceptance rate measured
void tune_step_size(
    float accepance_rate, // acceptance rate so far
    float* step_size      // reference to step_size to be tuned
    );

#endif
