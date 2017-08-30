#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include "Constants.h"
#include "MonteCarlo.h"

double tr1(REAL complex *Matrices, int pos1, int NUM_H);
double tr2(REAL complex *Matrices, int pos1, int pos2);
double tr3(REAL complex *Matrices, int pos1, int pos2, int pos3);

double tr2_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
double tr2_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
double tr3_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
double tr3_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
double row_norm_squared(REAL complex *Matrices, int pos1, int pos_row);

// Calculate the trace of a matrix H
double tr_H(
    REAL complex const *Matrix,
    uint64_t const length
    );

// Calculate the trace of a matrix H^2
double tr_H2(
    REAL complex const *Matrix,
    uint64_t const length
    );

// Calculate the action Tr D^2
double traceD2(
    REAL complex const *Matrices, // Array of all matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    uint64_t const length              // side length of all matrices
    );

// Calculate the action Tr D^4
double traceD4(
    REAL complex const *Matrices, // Array of all matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length              // side length of all matrices
    );

// Calculate the change of the action S = Tr D^2,
// if in one matrix one element is changed (t),
// where the matrix is in its already changed state
double deltaS_traceD2(
    REAL complex *Matrices,        // array of matrices
    const int num_h,               // number of matrices of H_TYPE
    const int num_l,               // number of matrices of L_TYPE
    const int length,              // side length N of one matrix in Matrices (the same for all matrices)
    const struct Matrix_State old  // old state
    );

double delta_action_traceD4(REAL complex *Matrices, int positionA, REAL complex temp, int pos_x, int pos_y, int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L);
#endif
