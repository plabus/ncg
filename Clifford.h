#ifndef __CLIFF_HEADER__
#define __CLIFF_HEADER__

#include "Matrix_Properties.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void printMint(int *matrix, int n);
void printMcomplex(float complex *matrix, int n);


int factorial(int n);
int binomial(int a, int b);
void combination(int* c,int n,int p, int x);


void commutator(float complex *A, float complex *B, float complex *C, int k);
void antisymmetrise(float complex *gammas, int dim, int k, int num_indices, int *sequence, float complex *matrix);
void Count_Hs_and_Ls(int p, int q, int *seq, int n, int *num_h, int *num_l);


void Generate_Gammas(int p, int q, float complex *gammas);
void Reshuffle_Clifford_Group(int p, int q, float complex *big_gammas, int num_h, int num_l, int ODD);
void Generate_Clifford_Group(int p, int q, float complex *big_gammas, int *num_h, int *num_l);
void Generate_Clifford_Odd_Group(int p, int q, float complex *big_gammas, int *num_h, int *num_l);

void Calculate_Trace_Gamma_ABAB(
    float complex const *Gamma_Matrices,
    struct Matrix_Properties const prop,
    int *sigmaAB
    );

int calculate_sigmaABCD(
    float complex const *Gamma_Matrices,
    size_t const a,
    size_t const b,
    size_t const c,
    size_t const d,
    struct Matrix_Properties const prop
    );

void Calculate_Trace_Gamma_ABCD(
    float complex const *Gamma_Matrices,
    struct Matrix_Properties const prop,
    int **sigmaABCD
    );
#endif
