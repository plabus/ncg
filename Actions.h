#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include "MonteCarlo.h"
#include "Precision.h"

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

// Calculates \sum_n Re( A_in B_nj )
double tr2_real_ij(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    size_t const pos1,                   // position of matrix A within Matrices array
    size_t const pos2,                   // position of matrix B within Matrices array
    size_t const pos_x,                  // one of the indices contributing to (i, j)
    size_t const pos_y                   // the other of the indices contributing to (i, j)
    );

// Calculates \sum_n Im( A_in B_nj )
double tr2_imag_ij(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    size_t const pos1,                   // position of matrix A within Matrices array
    size_t const pos2,                   // position of matrix B within Matrices array
    size_t const pos_x,                  // one of the indices contributing to (i, j)
    size_t const pos_y                   // the other of the indices contributing to (i, j)
    );

// Calculates \sum_n \sum_m Re( A_in B_nm C_mj )
double tr3_real_ij(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    size_t const pos1,                   // position of matrix A within Matrices array
    size_t const pos2,                   // position of matrix B within Matrices array
    size_t const pos3,                   // position of matrix C within Matrices array
    size_t const pos_x,                  // one of the indices contributing to (i, j)
    size_t const pos_y                   // the other of the indices contributing to (i, j)
    );

// Calculates \sum_n \sum_m Im( A_in B_nm C_mj )
double tr3_imag_ij(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    size_t const pos1,                   // position of matrix A within Matrices array
    size_t const pos2,                   // position of matrix B within Matrices array
    size_t const pos3,                   // position of matrix C within Matrices array
    size_t const pos_x,                  // one of the indices contributing to (i, j)
    size_t const pos_y                   // the other of the indices contributing to (i, j)
    );

// Calculates \sum_n |A_in|^2
double row_norm_squared(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    size_t const pos,                    // position of matrix A within Matrices array
    size_t const pos_row                 // index i, indicating the row of the matrix to be squared
    );

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
    REAL complex const *Matrices,       // Array of all matrices
    struct Matrix_Properties const prop // includes num_h, num_l, n and k
    );

// Calculate the action Tr D^4
double traceD4(
    REAL complex const *Matrices,        // Array of all matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    );

// Calculate the change of the action S = Tr D^2,
// if in one matrix one element is changed (t),
// where the matrix is in its already changed state
double delta_traceD2(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    struct Matrix_State const old        // old state
    );

// Calculate the change of the action S = Tr D^2,
// if in one matrix one element is changed ( t = a + I * b ),
// where the matrix is in its already changed state
double delta_traceD4(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    const struct Matrix_State old_state, // struct containing information about old state
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    );

// Wrapper function for the full action:
// calculates S = g2 * Tr(D^2) + g4 * Tr(D^4)
double Calculate_Action(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    );

// Wrapper function for the change in the action:
// calculates delta_S = g2 * Tr(delta_D^2) + g4 * Tr(delta_D^4)
double Calculate_Delta_Action(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD,                     // pre-calculated Clifford products of 4 Gamma matrices
    const struct Matrix_State old        // struct containing information about old state
    );
#endif
