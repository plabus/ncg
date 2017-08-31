#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include "Constants.h"
#include "MonteCarlo.h"

// Calculates the trace of a single matrix in the array Matrices at position pos1
double tr1(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1                // position of the matrix inside the Matrices array
    );

// Calculates the trace of the product of two matrix in the array Matrices
// at positions pos1 and pos2
double tr2(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1,               // position of the 1st matrix inside the Matrices array
    int const pos2                // position of the 2nd matrix inside the Matrices array
    );

// Calculates the trace of the product of three matrix in the array Matrices
// at positions pos1, pos2 and pos3
double tr3(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1,               // position of the 1st matrix inside the Matrices array
    int const pos2,               // position of the 2nd matrix inside the Matrices array
    int const pos3                // position of the 3rd matrix inside the Matrices array
    );

// Calculates the trace of the product of four matrix in the array Matrices
// at positions pos1, pos2, pos3 and pos4
double tr4(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1,               // position of the 1st matrix inside the Matrices array
    int const pos2,               // position of the 2nd matrix inside the Matrices array
    int const pos3,               // position of the 3rd matrix inside the Matrices array
    int const pos4                // position of the 4th matrix inside the Matrices array
    );

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
    int const n,                  // side length of all matrices
    int const k                   // dimension k of gamma matrices
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
