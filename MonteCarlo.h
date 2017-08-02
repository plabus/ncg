#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"

double cclock();
void printM(const int L, float complex * M);
int compare_floats(const void * a, const void * b);
float delta_approx(float);
void nullify(const int l, float complex * m);


void Matrices_Initialisation(struct pcg32_random_t *rng, float complex *Matrices, float *action, int NUM_H, int NUM_L);
void Get_Next_MCMC_Element(struct pcg32_random_t *rng, float complex *Matrices, float *action,
                           int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L, int *acc);

void Arrange_Anticommutator(float complex * m, float complex * acomm);
void Arrange_Commutator(float complex * m, float complex * comm);
void Arrange_Dirac_Matrix(float complex *gamma_passed, float complex *Matrices, float complex *Matrix_Operators, float complex *Dirac, int NUM_H, int NUM_L);

void Measure_Eigenvalues_Dirac(float complex *Gamma_Matrices, float complex *Matrices,
		               float complex *Matrix_Operators, float complex *Dirac,
		               float *evs_D, float *evs_D_avrg, float *evs_D_avrg2,
                               int NUM_H, int NUM_L);
void Measure_Eigenvaluedistribution_Dirac(float *support_points, float *evs_D,
                                          float *dist_evs_D_avrg);
void Measure_Orderparameter_Frac(float complex *Matrices, double *frac,
                                 double *frac_squared, int NUM_H);

#endif
