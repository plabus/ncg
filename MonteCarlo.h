#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "Random.h"

double cclock();
void printM(const int L, MKL_Complex8 * M);
int compare_floats(const void * a, const void * b);
float delta_approx(float);
void nullify(const int l, MKL_Complex8 * m);
 

void Matrices_Initialisation(struct pcg32_random_t *rng, MKL_Complex8 *Matrices, float *action, int NUM_H, int NUM_L);
void Get_Next_MCMC_Element(struct pcg32_random_t *rng, MKL_Complex8 *Matrices, float *action,
                           int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L, int *acc);

void Arrange_Anticommutator(MKL_Complex8 * m, MKL_Complex8 * acomm);
void Arrange_Commutator(MKL_Complex8 * m, MKL_Complex8 * comm);
void Arrange_Dirac_Matrix(float complex *gamma_passed, MKL_Complex8 *Matrices, MKL_Complex8 *Matrix_Operators, MKL_Complex8 *Dirac, int NUM_H, int NUM_L);

void Measure_Eigenvalues_Dirac(float complex *Gamma_Matrices, MKL_Complex8 *Matrices,
		               MKL_Complex8 *Matrix_Operators, MKL_Complex8 *Dirac,
		               float *evs_D, float *evs_D_avrg, float *evs_D_avrg2,
                               int NUM_H, int NUM_L);
void Measure_Eigenvaluedistribution_Dirac(float *support_points, float *evs_D,
                                          float *dist_evs_D_avrg);
void Measure_Orderparameter_Frac(MKL_Complex8 *Matrices, double *frac,
                                 double *frac_squared, int NUM_H);
  
#endif
