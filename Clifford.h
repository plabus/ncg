#ifndef __CLIFF_HEADER__
#define __CLIFF_HEADER__

#include "Matrix_Properties.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// enum to distinguish even, odd and full Clifford group
enum Clifford_Type { CLIFFORD_ODD, CLIFFORD_EVEN, CLIFFORD_FULL };


void printMint(int *matrix, int n);
void printMcomplex(float complex *matrix, int n);

int factorial(int n);
int binomial(int a, int b);
void combination(int* c,int n,int p, int x);

void commutator(float complex *A, float complex *B, float complex *C, int k);
void antisymmetrise(float complex *gammas, int dim, int k, int num_indices, int *sequence, float complex *matrix);
void Count_Hs_and_Ls(
    size_t p,
    size_t q,
    int *seq,
    size_t n,
    size_t *num_h,
    size_t *num_l
    );


void Generate_Gammas(
    float complex *gammas,
    struct Matrix_Properties const prop
    );

void Reshuffle_Clifford_Group(
    float complex *big_gammas,
    struct Matrix_Properties const prop,
    enum Clifford_Type const group_type
    );

void Generate_Clifford_Group(
    float complex *big_gammas,
    struct Matrix_Properties* prop
    );

void Generate_Clifford_Odd_Group(
    float complex *big_gammas,
    struct Matrix_Properties* prop
    );

void Calculate_Trace_Gamma_ABAB(
    float complex const *Gamma_Matrices,
    struct Matrix_Properties const prop,
    int *sigmaAB
    );

int calculate_sigmaABCD(
    float complex const *Gamma_Matrices,
    struct Matrix_Properties const prop,
    size_t const a,
    size_t const b,
    size_t const c,
    size_t const d
    );

void Calculate_Trace_Gamma_ABCD(
    float complex const *Gamma_Matrices,
    struct Matrix_Properties const prop,
    int **sigmaABCD
    );
#endif
