#ifndef _MEASUREMENTS_H_
#define _MEASUREMENTS_H_


// Given a (NxN) matrix m compute {m,.}
void Arrange_Anticommutator(
    float complex * m,
    float complex * acomm
    );

// Given a (NxN) matrix m compute [m,.]
void Arrange_Commutator(
    float complex * m,
    float complex * comm
    );

// Arrange Matrix D (kn^2 x kn^2) from NUM_H Matrices H and NUM_L Matrices L of type (P,Q)
void Arrange_Dirac_Matrix(
    float complex *gamma_passed,
    float complex *Matrices,
    float complex *Matrix_Operators,
    float complex *Dirac,
    int NUM_H,
    int NUM_L
    );


void Measure_Eigenvalues_Dirac(
    float complex *Gamma_Matrices,
    float complex *Matrices,
    float complex *Matrix_Operators,
    float complex *Dirac,
    float *evs_D,
    float *evs_D_avrg,
    float *evs_D_avrg2,
    int NUM_H,
    int NUM_L
    );

void Measure_Eigenvaluedistribution_Dirac(
    float *support_points,
    float *evs_D,
    float *dist_evs_D_avrg
    );

void Measure_Orderparameter_Frac(
    float complex *Matrices,
    double *frac,
    double *frac_squared,
    int NUM_H
    );

#endif
