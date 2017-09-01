#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <inttypes.h>
#include "Actions.h"
#include "Constants.h"

/**********************************************************************************************/

// Calculates the trace of a single matrix in the array Matrices at position pos1
double tr1(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1                // position of the matrix inside the Matrices array
    )
{
  // Quickly returning for trace-less matrices
  if( pos1 >= num_h )
  {
    return 0.0;
  }

  uint64_t const offset1 = pos1 * length * length;
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    trace += creal( Matrices[ offset1 + i * length + i ] );
  }

  return trace;
}

// Calculates the trace of the product of two matrix in the array Matrices
// at positions pos1 and pos2
double tr2(
    REAL complex const *Matrices, // array of matrices
    int const num_h,              // number of matrices of H_TYPE
    int const num_l,              // number of matrices of L_TYPE
    int const length,             // side length N of one matrix in Matrices (the same for all matrices)
    int const pos1,               // position of the 1st matrix inside the Matrices array
    int const pos2                // position of the 2nd matrix inside the Matrices array
    )
{
  uint64_t const offset1 = pos1 * length * length;
  uint64_t const offset2 = pos2 * length * length;
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    for( uint64_t j = 0; j < length; ++j )
    {
      trace += creal( Matrices[ offset1 + i * length + j ] ) * creal( Matrices[ offset2 + j * length + i ] );
      trace -= cimag( Matrices[ offset1 + i * length + j ] ) * cimag( Matrices[ offset2 + j * length + i ] );
    }
  }

  return trace;
}

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
    )
{
  uint64_t const offset1 = pos1 * length * length;
  uint64_t const offset2 = pos2 * length * length;
  uint64_t const offset3 = pos3 * length * length;
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    for( uint64_t j = 0; j < length; ++j )
    {
      for( uint64_t k = 0; k < length; ++k )
      {
        trace += creal( Matrices[ offset1 + i * length + j ] ) *
                 creal( Matrices[ offset2 + j * length + k ] ) *
                 creal( Matrices[ offset3 + k * length + i ] );
        trace -= cimag( Matrices[ offset1 + i * length + j ] ) *
                 cimag( Matrices[ offset2 + j * length + k ] ) *
                 creal( Matrices[ offset3 + k * length + i ] );
        trace -= creal( Matrices[ offset1 + i * length + j ] ) *
                 cimag( Matrices[ offset2 + j * length + k ] ) *
                 cimag( Matrices[ offset3 + k * length + i ] );
        trace -= cimag( Matrices[ offset1 + i * length + j ] ) *
                 creal( Matrices[ offset2 + j * length + k ] ) *
                 cimag( Matrices[ offset3 + k * length + i ] );
      }
    }
  }

  return trace;
}

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
    )
{
  uint64_t const offset1 = pos1 * length * length;
  uint64_t const offset2 = pos2 * length * length;
  uint64_t const offset3 = pos3 * length * length;
  uint64_t const offset4 = pos4 * length * length;
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    for( uint64_t j = 0; j < length; ++j )
    {
      for( uint64_t k = 0; k < length; ++k )
      {
        for( uint64_t l = 0; l < length; ++l )
        {
          trace += creal( Matrices[ offset1 + i * length + j ] ) *
                   creal( Matrices[ offset2 + j * length + k ] ) *
                   creal( Matrices[ offset3 + k * length + l ] ) *
                   creal( Matrices[ offset4 + l * length + i ] );

          trace += cimag( Matrices[ offset1 + i * length + j ] ) *
                   cimag( Matrices[ offset2 + j * length + k ] ) *
                   cimag( Matrices[ offset3 + k * length + l ] ) *
                   cimag( Matrices[ offset4 + l * length + i ] );

          trace -= creal( Matrices[ offset1 + i * length + j ] ) *
                   creal( Matrices[ offset2 + j * length + k ] ) *
                   cimag( Matrices[ offset3 + k * length + l ] ) *
                   cimag( Matrices[ offset4 + l * length + i ] );

          trace -= creal( Matrices[ offset1 + i * length + j ] ) *
                   cimag( Matrices[ offset2 + j * length + k ] ) *
                   creal( Matrices[ offset3 + k * length + l ] ) *
                   cimag( Matrices[ offset4 + l * length + i ] );

          trace -= cimag( Matrices[ offset1 + i * length + j ] ) *
                   creal( Matrices[ offset2 + j * length + k ] ) *
                   creal( Matrices[ offset3 + k * length + l ] ) *
                   cimag( Matrices[ offset4 + l * length + i ] );

          trace -= creal( Matrices[ offset1 + i * length + j ] ) *
                   cimag( Matrices[ offset2 + j * length + k ] ) *
                   cimag( Matrices[ offset3 + k * length + l ] ) *
                   creal( Matrices[ offset4 + l * length + i ] );

          trace -= cimag( Matrices[ offset1 + i * length + j ] ) *
                   creal( Matrices[ offset2 + j * length + k ] ) *
                   cimag( Matrices[ offset3 + k * length + l ] ) *
                   creal( Matrices[ offset4 + l * length + i ] );

          trace -= cimag( Matrices[ offset1 + i * length + j ] ) *
                   cimag( Matrices[ offset2 + j * length + k ] ) *
                   creal( Matrices[ offset3 + k * length + l ] ) *
                   creal( Matrices[ offset4 + l * length + i ] );
        }
      }
    }
  }

  return trace;
}

/**********************************************************************************************/

double tr2_real_ij(REAL complex const *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
  double sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    sum += creal( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+j+off2] );
    sum -= cimag( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+j+off2] );
  }

  return sum;
}


double tr2_imag_ij(REAL complex const *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
  double sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    sum += creal( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+j+off2] );
    sum += cimag( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+j+off2] );
  }

  return sum;
}


double tr3_real_ij(REAL complex const *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
  double sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;
  int off3 = pos3*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    for(int c=0;c<N;c++) {
      sum += creal( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+c+off2] ) * creal( Matrices[c*N+j+off3] );
      sum -= cimag( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+c+off2] ) * creal( Matrices[c*N+j+off3] );
      sum -= creal( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+c+off2] ) * cimag( Matrices[c*N+j+off3] );
      sum -= cimag( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+c+off2] ) * cimag( Matrices[c*N+j+off3] );
    }
  }

  return sum;
}


double tr3_imag_ij(REAL complex const *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
  double sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;
  int off3 = pos3*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    for(int c=0;c<N;c++) {
      sum -= cimag( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+c+off2] ) * cimag( Matrices[c*N+j+off3] );
      sum += creal( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+c+off2] ) * cimag( Matrices[c*N+j+off3] );
      sum += creal( Matrices[i*N+b+off1] ) * cimag( Matrices[b*N+c+off2] ) * creal( Matrices[c*N+j+off3] );
      sum += cimag( Matrices[i*N+b+off1] ) * creal( Matrices[b*N+c+off2] ) * creal( Matrices[c*N+j+off3] );
    }
  }

  return sum;
}


double matrix_norm_squared(REAL complex const *Matrices, int pos1) {
  double sum_off_diag = 0.f;
  double trace_squared = 0.f;
  int off1 = pos1*SWEEP;

  for(int a=0;a<N;a++) {
    for(int b=a+1;b<N;b++) {
      sum_off_diag += creal( Matrices[a*N+b+off1] ) * creal( Matrices[a*N+b+off1] );
      sum_off_diag += cimag( Matrices[a*N+b+off1] ) * cimag( Matrices[a*N+b+off1] );
    }
  }

  for(int a=0;a<N;a++) trace_squared += creal( Matrices[a*N+a+off1] ) * creal( Matrices[a*N+a+off1] );

  return 2*sum_off_diag+trace_squared;
}


double row_norm_squared(REAL complex const *Matrices, int pos1, int pos_row) {
  double sum = 0.f;
  int off = pos1*SWEEP + pos_row*N;

  for(int a=0;a<N;a++) {
    sum += creal( Matrices[off+a] ) * creal( Matrices[off+a] );
    sum += cimag( Matrices[off+a] ) * cimag( Matrices[off+a] );
  }

  return sum;
}

/**********************************************************************************************/

// Calculate the trace of a matrix H
double tr_H(
    REAL complex const *Matrix,
    uint64_t const length
    )
{
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    trace += creal( Matrix[ i * length + i ] );
  }

  return trace;
}

// Calculate the trace of a matrix H^2
double tr_H2(
    REAL complex const *Matrix,
    uint64_t const length
    )
{
  double trace = 0.0;

  for( uint64_t i = 0; i < length; ++i )
  {
    for( uint64_t j = 0; j < length; ++j )
    {
      trace += creal( Matrix[ i * length + j ] * Matrix[ j * length + i] );
    }
  }

  return trace;
}

/**********************************************************************************************/

/*********************************************************************************/
/********                                                                *********/
/********                      ACTION S = Tr D^2                         *********/
/********                                                                *********/
/*********************************************************************************/

// Calculate the full action S = Tr D^2
double traceD2(
    REAL complex const *Matrices,       // Array of all matrices
    struct Matrix_Properties const prop // includes num_h, num_l, n and k
    )
{
  //------------------------------------------------------//
  // Calculating the action:                              //
  // S = Tr D^2 = 2 k [ N \sum Tr(H^2) + \sum (Tr H)^2 ]  //
  //------------------------------------------------------//

  // Unpack property parameters
  size_t const num_h = prop.num_h;
  size_t const num_l = prop.num_l;
  size_t const n     = prop.n;
  size_t const k     = prop.k;

  // For the Tr H^2 part
  // -------------------
  double sum_trace_H2 = 0.0;

  for( uint64_t m = 0; m < num_h + num_l; ++m )  // This is for all matrices H and L
  {
    uint64_t const offset = m * n * n;
    sum_trace_H2 += tr_H2( &Matrices[offset], n );
  }

  // For the (Tr H)^2 part
  // ---------------------
  double sum_trace_H_sq  = 0.0;

  for( uint64_t m = 0; m < num_h; ++m )  // This is for matrices H only
  {
    uint64_t const offset = m * n * n;
    double const trace_H = tr_H( &Matrices[offset], n );
    sum_trace_H_sq += trace_H * trace_H;
  }

  return 2 * k * ( n * sum_trace_H2 + sum_trace_H_sq );
}

// Calculate the change of the action S = Tr D^2,
// if in one matrix one element is changed (t),
// where the matrix is in its already changed state
double delta_traceD2(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    struct Matrix_State const old        // old state
    )
{
  // Unpack needed information
  const size_t n              = prop.n;
  const size_t k              = prop.k;
  const enum Matrix_Type type = old.matrix < prop.num_h ? H_TYPE : L_TYPE;
  const size_t offset            = old.matrix * n * n;
  const size_t pos_upper         = old.pos_upper;
  const size_t pos_lower         = old.pos_lower;

  // Unpack matrix elements involved in calculating delta S,
  // where t = H_new - H_old (at upper triangular position)
  const double t_Re = creal( Matrices[ offset + pos_upper ] - old.matrix_element );
  const double t_Im = cimag( Matrices[ offset + pos_upper ] - old.matrix_element );
  const double H_new_ij_Re = creal( Matrices[ offset + pos_upper ] );
  const double H_new_ij_Im = cimag( Matrices[ offset + pos_upper ] );

  // Off-diagonal case:
  // ------------------
  if( pos_upper != pos_lower )
  {
    return 4 * k * n * (
        2 * ( t_Re * H_new_ij_Re + t_Im * H_new_ij_Im ) - ( t_Re * t_Re + t_Im * t_Im )
        );
  }
  // Diagonal case:
  // --------------
  else
  {
    // H_TYPE:
    if( type == H_TYPE )
    {
      const double tr_H_new = tr_H( &Matrices[offset], n );
      return 2 * k * (
          n * ( 2 * t_Re * H_new_ij_Re - t_Re * t_Re ) + ( 2 * t_Re * tr_H_new - t_Re * t_Re )
          );
    }
    // L_TYPE:
    else
    {
      return 2 * k * n * (
          2 * ( t_Re * H_new_ij_Re + t_Im * H_new_ij_Im ) - ( t_Re * t_Re + t_Im * t_Im )
          );
    }
  }
}


/*********************************************************************************/
/********                                                                *********/
/********                      ACTION S = Tr D^4                         *********/
/********                                                                *********/
/*********************************************************************************/

double traceD4(
    REAL complex const *Matrices,        // Array of all matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    )
{
  //-----------------------------------------------------------------------------//
  // Calculating the action:                                                     //
  //                                                                             //
  // S = Tr D^4                                                                  //
  //   = 2 k \sum_A [ N Tr(A^4) + 4 Tr(A) Tr(A^3) + 3 (Tr(A^2))^2 ]  ---> part I //
  //                                                                             //
  //-----------------------------------------------------------------------------//

  // Unpack property parameters
  size_t const num_h = prop.num_h;
  size_t const num_l = prop.num_l;
  size_t const num_m = num_h + num_l;
  size_t const n     = prop.n;
  size_t const k     = prop.k;

  double partI = 0.0;
  for( size_t posA = 0; posA < num_m; ++posA )
  {
    double const TrA1 = tr1( Matrices, num_h, num_l, n, posA );
    double const TrA2 = tr2( Matrices, num_h, num_l, n, posA, posA );
    double const TrA3 = tr3( Matrices, num_h, num_l, n, posA, posA, posA );
    double const TrA4 = tr4( Matrices, num_h, num_l, n, posA, posA, posA, posA );

    partI += n * TrA4 + 4.0 * TrA1 * TrA3 + 3.0 * TrA2 * TrA2;
  }

  double partII = 0.0;
  for( size_t posB = 0; posB < num_m; ++posB )
  {
    for( size_t posA = 0; posA < posB; ++posA )
    {
      double const TrA    = tr1( Matrices, num_h, num_l, n, posA );
      double const TrB    = tr1( Matrices, num_h, num_l, n, posB );
      double const TrAA   = tr2( Matrices, num_h, num_l, n, posA, posA );
      double const TrAB   = tr2( Matrices, num_h, num_l, n, posA, posB );
      double const TrBB   = tr2( Matrices, num_h, num_l, n, posB, posB );
      double const TrAAB  = tr3( Matrices, num_h, num_l, n, posA, posA, posB );
      double const TrABB  = tr3( Matrices, num_h, num_l, n, posA, posB, posB );
      double const TrAABB = tr4( Matrices, num_h, num_l, n, posA, posA, posB, posB );
      double const TrABAB = tr4( Matrices, num_h, num_l, n, posA, posB, posA, posB );

      int const sgnA = posA < num_h ? 1 : -1;
      int const sgnB = posB < num_h ? 1 : -1;
      int const sAB = sigmaAB[ posA * num_m + posB ];

      partII += n * ( 2.0 * TrAABB + sAB * TrABAB )
              + ( 2.0 + sAB )
              * ( 2.0 * TrA * TrABB + 2.0 * TrB * TrAAB + 2.0 * sgnA * sgnB * TrAB * TrAB + TrAA * TrBB );
    }
  }

  double const action = 2.0 * k * ( partI + 2.0 * partII );
  return action;
}

// Calculate the change of the action S = Tr D^2,
// if in one matrix one element is changed ( t = a + I * b ),
// where the matrix is in its already changed state
double delta_traceD4(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    const struct Matrix_State old_state, // struct containing information about old state
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    )
{
  /*******************************************************
   *                                                     *
   *                         INIT                        *
   *                                                     *
   *******************************************************/

  // Unpack property parameters
  const size_t num_h = prop.num_h;
  const size_t num_l = prop.num_l;
  const size_t num_m = num_h + num_l;
  const size_t n     = prop.n;
  const size_t k     = prop.k;
  const size_t positionA = old_state.matrix; // number of the matrix that has changed
  const size_t pos_upper = old_state.pos_upper; // index in uppper triangle where element was changed
  // FIXME: Did this come in a specific order????
  const size_t pos_x = old_state.pos_x;
  const size_t pos_y = old_state.pos_y;
  const REAL complex old = old_state.matrix_element; // the old matrix element at that position

  const size_t off = positionA * n * n;
  const int sgnA = positionA < num_h ? 1 : -1;
  const size_t num_trABCD = sigmaABCD[positionA][0];

  // FIXME:
  // HERE WAS DEFINED
  // delta A = A_new - A_old == a + b*i
  //       A = A_old
  //       A_ij = A_old_ij == c + d*i
  // TODO:
  // Add the relevant part of the documentation as comment

  // Precalculate some variables
  const double a = creal( Matrices[ off + pos_upper ] - old );
  const double b = cimag( Matrices[ off + pos_upper ] - old );
  const double c = creal( Matrices[ off + pos_upper ] );
  const double d = cimag( Matrices[ off + pos_upper ] );

  const double trA  = tr1(Matrices, num_h, num_l, n, positionA);
  const double trA2 = tr2(Matrices, num_h, num_l, n, positionA, positionA);
  const double A_ii = creal( Matrices[ off + pos_x * n + pos_x ] );
  const double A_jj = creal( Matrices[ off + pos_y * n + pos_y ] );

  const double abs2    = a * a + b * b;
  const double mixsum  = a * c + b * d;

  // Initialise the accumulation variables
  double delta1 = 0.0;
  double delta2 = 0.0;
  double delta3 = 0.0;
  double sum[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  /*******************************************************
   *                                                     *
   *                    OFF-DIAGONAL                     *
   *                                                     *
   *******************************************************/

  if( pos_x != pos_y )
  {

    /*******************************/
    /**********  PART I  ***********/
    /*******************************/

    sum[0] = tr3_real_ij(Matrices, positionA, positionA, positionA, pos_x, pos_y);
    sum[1] = tr3_imag_ij(Matrices, positionA, positionA, positionA, pos_x, pos_y);
    sum[2] = row_norm_squared(Matrices, positionA, pos_x);
    sum[3] = row_norm_squared(Matrices, positionA, pos_y);
    if( sgnA == 1 )
    {
      sum[4] = tr2_real_ij(Matrices, positionA, positionA, pos_x, pos_y);
      sum[5] = tr2_imag_ij(Matrices, positionA, positionA, pos_x, pos_y);
    }
    else
    {
      sum[4] = 0.0;
      sum[5] = 0.0;
    }

    const double F = 4.0 * mixsum - 2.0 * abs2;

    const double part_Ia = 8.0 * ( a * sum[0] + b * sum[1] )
                         - 4.0 * (
                                   ( a * a - b * b ) * ( c * c - d * d )
                                   + 4.0 * a * b * c * d
                                   + abs2 * A_ii * A_jj
                                 )
                         - 4.0 * abs2 * ( sum[2] + sum[3] )
                         + 8.0 * abs2 * mixsum
                         - 2.0 * abs2 * abs2;
    const double part_Ib = 12.0 * trA * ( 2.0 * a * sum[4] + 2.0 * b * sum[5] - abs2 * ( A_ii + A_jj ) );
    const double part_Ic = 3.0 * F * ( 2.0 * trA2 - F );

    delta1 = n * part_Ia + part_Ib + part_Ic;

    /*******************************/
    /**********  PART II  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES B != A */
    for( size_t positionB = 0; positionB < num_m; positionB++ )
    {

      if( positionB == positionA ) continue;

      const int sAB = sigmaAB[ positionA * num_m + positionB ];
      const int sgnB = positionB < num_h ? 1 : -1;
      const size_t offB = positionB * n * n;

      const double trB  = tr1(Matrices, num_h, num_l, n, positionB);
      const double trB2 = tr2(Matrices, num_h, num_l, n, positionB, positionB);
      const double trAB = tr2(Matrices, num_h, num_l, n, positionA, positionB);

      const double c_B  = creal( Matrices[ offB + pos_upper ] );
      const double d_B  = cimag( Matrices[ offB + pos_upper ] );
      const double B_ii = creal( Matrices[ offB + pos_x * n + pos_x ] );
      const double B_jj = creal( Matrices[ offB + pos_y * n + pos_y ] );

      sum[0] = tr3_real_ij( Matrices, positionA, positionB, positionB, pos_x, pos_y );
      sum[0]+= tr3_real_ij( Matrices, positionB, positionB, positionA, pos_x, pos_y );
      sum[1] = tr3_imag_ij( Matrices, positionA, positionB, positionB, pos_x, pos_y );
      sum[1]+= tr3_imag_ij( Matrices, positionB, positionB, positionA, pos_x, pos_y );
      sum[2] = row_norm_squared( Matrices, positionB, pos_x );
      sum[3] = row_norm_squared( Matrices, positionB, pos_y );
      sum[4] = tr3_real_ij( Matrices, positionB, positionA, positionB, pos_x, pos_y );
      sum[5] = tr3_imag_ij( Matrices, positionB, positionA, positionB, pos_x, pos_y );
      if(sgnA==1)
      {
        sum[6] = tr2_real_ij(Matrices, positionB, positionB, pos_x, pos_y);
        sum[7] = tr2_imag_ij(Matrices, positionB, positionB, pos_x, pos_y);
      }
      else
      {
        sum[6] = 0.0;
        sum[7] = 0.0;
      }
      if(sgnB==1)
      {
        sum[8] = tr2_real_ij(Matrices, positionA, positionB, pos_x, pos_y);
        sum[8]+= tr2_real_ij(Matrices, positionB, positionA, pos_x, pos_y);
        sum[9] = tr2_imag_ij(Matrices, positionA, positionB, pos_x, pos_y);
        sum[9]+= tr2_imag_ij(Matrices, positionB, positionA, pos_x, pos_y);
      }
      else
      {
        sum[8] = 0.0;
        sum[9] = 0.0;
      }

      const double part_IIa = 2 * a * sum[0] + 2 * b * sum[1] + abs2 * ( sum[2] + sum[3] );
      const double part_IIb = 4 * a * sum[4] + 4 * b * sum[5] +
                 2 * ( (a * a - b * b) * (c_B * c_B - d_B * d_B) + 4 * a * b * c_B * d_B + abs2 * B_ii * B_jj );
      const double part_IIc = 2 * trA * ( a * sum[6] + b * sum[7] ) +
                 trB * ( 2 * a * sum[8] + 2 * b * sum[9] + abs2 * ( B_ii + B_jj ) ) +
                 4 * sgnA * sgnB * ( a * c_B + b * d_B ) * ( trAB + a * c_B + b * d_B ) +
                 trB2 * ( 2 * mixsum + abs2 );

      delta2 += 4 * n * part_IIa + 2 * sAB * n * part_IIb + ( 8 + 4 * sAB ) * part_IIc;

    } /* END LOOP OVER MATRICES B */

    /*******************************/
    /*********  PART III  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES D != C != B != A pairwise, with   *
     * non-vanishing (pre-calculated) traces Tr(A*B*C*D) != 0   */
    for(int element=0; element<num_trABCD; element++) {

      const size_t positionB = sigmaABCD[positionA][9*element+1];
      const size_t positionC = sigmaABCD[positionA][9*element+2];
      const size_t positionD = sigmaABCD[positionA][9*element+3];
      const int b1           = sigmaABCD[positionA][9*element+4];
      const int b2           = sigmaABCD[positionA][9*element+5];
      const int c1           = sigmaABCD[positionA][9*element+6];
      const int c2           = sigmaABCD[positionA][9*element+7];
      const int d1           = sigmaABCD[positionA][9*element+8];
      const int d2           = sigmaABCD[positionA][9*element+9];
      const int s = b1 + b2 + c1 + c2 + d1 + d2;

      const size_t offB = positionB * n * n;
      const size_t offC = positionC * n * n;
      const size_t offD = positionD * n * n;

      const int sgnB = positionB < num_h ? 1 : -1;
      const int sgnC = positionC < num_h ? 1 : -1;
      const int sgnD = positionD < num_h ? 1 : -1;

      const double c_B = creal( Matrices[pos_upper+offB] );
      const double d_B = cimag( Matrices[pos_upper+offB] );
      const double c_C = creal( Matrices[pos_upper+offC] );
      const double d_C = cimag( Matrices[pos_upper+offC] );
      const double c_D = creal( Matrices[pos_upper+offD] );
      const double d_D = cimag( Matrices[pos_upper+offD] );

      /* Traces over one matrix: */
      const double trB = tr1(Matrices, num_h, num_l, n, positionB);
      const double trC = tr1(Matrices, num_h, num_l, n, positionC);
      const double trD = tr1(Matrices, num_h, num_l, n, positionD);

      /* Traces over two matrices: */
      const double trABtrCD = (a * c_B + b * d_B) * tr2(Matrices, num_h, num_l, n, positionC, positionD);
      const double trACtrBD = (a * c_C + b * d_C) * tr2(Matrices, num_h, num_l, n, positionB, positionD);
      const double trADtrBC = (a * c_D + b * d_D) * tr2(Matrices, num_h, num_l, n, positionB, positionC);

      /* Traces over one and three matrices: */
      double trF2 = 0.0;
      double trF3 = 0.0;
      double trF4 = 0.0;
      if( sgnB == 1 )
      {
        trF2 = a*tr2_real_ij(Matrices, positionC, positionD, pos_x, pos_y);
        trF2+= a*tr2_real_ij(Matrices, positionD, positionC, pos_x, pos_y);
        trF2+= b*tr2_imag_ij(Matrices, positionC, positionD, pos_x, pos_y);
        trF2+= b*tr2_imag_ij(Matrices, positionD, positionC, pos_x, pos_y);
      }
      if( sgnC == 1 )
      {
        trF3 = a*tr2_real_ij(Matrices, positionB, positionD, pos_x, pos_y);
        trF3+= a*tr2_real_ij(Matrices, positionD, positionB, pos_x, pos_y);
        trF3+= b*tr2_imag_ij(Matrices, positionB, positionD, pos_x, pos_y);
        trF3+= b*tr2_imag_ij(Matrices, positionD, positionB, pos_x, pos_y);
      }
      if( sgnD == 1 )
      {
        trF4 = a*tr2_real_ij(Matrices, positionB, positionC, pos_x, pos_y);
        trF4+= a*tr2_real_ij(Matrices, positionC, positionB, pos_x, pos_y);
        trF4+= b*tr2_imag_ij(Matrices, positionB, positionC, pos_x, pos_y);
        trF4+= b*tr2_imag_ij(Matrices, positionC, positionB, pos_x, pos_y);
      }

      /* Traces over four matrices: */
      double trAE = (b1+d2)*a*tr3_real_ij(Matrices, positionB, positionC, positionD, pos_x, pos_y);
      trAE+= (b1+d2)*a*tr3_real_ij(Matrices, positionD, positionC, positionB, pos_x, pos_y);
      trAE+= (b2+c2)*a*tr3_real_ij(Matrices, positionB, positionD, positionC, pos_x, pos_y);
      trAE+= (b2+c2)*a*tr3_real_ij(Matrices, positionC, positionD, positionB, pos_x, pos_y);
      trAE+= (c1+d1)*a*tr3_real_ij(Matrices, positionC, positionB, positionD, pos_x, pos_y);
      trAE+= (c1+d1)*a*tr3_real_ij(Matrices, positionD, positionB, positionC, pos_x, pos_y);

      trAE+= (b1+d2)*b*tr3_imag_ij(Matrices, positionB, positionC, positionD, pos_x, pos_y);
      trAE+= (b1+d2)*b*tr3_imag_ij(Matrices, positionD, positionC, positionB, pos_x, pos_y);
      trAE+= (b2+c2)*b*tr3_imag_ij(Matrices, positionB, positionD, positionC, pos_x, pos_y);
      trAE+= (b2+c2)*b*tr3_imag_ij(Matrices, positionC, positionD, positionB, pos_x, pos_y);
      trAE+= (c1+d1)*b*tr3_imag_ij(Matrices, positionC, positionB, positionD, pos_x, pos_y);
      trAE+= (c1+d1)*b*tr3_imag_ij(Matrices, positionD, positionB, positionC, pos_x, pos_y);

#ifdef LARGE_N
      delta3 +=   n*trAE;
#else
      delta3 += ( n*trAE + s*(trB*trF2 + trC*trF3 + trD*trF4) +
                  2*s*(sgnA*sgnB*trABtrCD + sgnA*sgnC*trACtrBD + sgnA*sgnD*trADtrBC) );
#endif
    } /* END LOOP OVER MATRICES */
  } /* END OFF-DIAGONAL CASE */

  /*******************************************************
   *                                                     *
   *                       DIAGONAL                      *
   *                                                     *
   *******************************************************/

  else
  {

    /*******************************/
    /**********  PART I  ***********/
    /*******************************/

    sum[0] = row_norm_squared(Matrices, positionA, pos_x);
    sum[1] = tr3_real_ij(Matrices, positionA, positionA, positionA, pos_x, pos_x);
    if( sgnA == 1 )
    {
      sum[2] = tr3(Matrices, num_h, num_l, n, positionA, positionA, positionA);
    }
    else
    {
      sum[2] = 0.0;
    }

    const double F = 2.0 * a * c - a * a;

    const double part_Ia = 4.0 * a * sum[1]
                         - 2.0 * a * a * c * c
                         - 4.0 * a * a * sum[0]
                         + 4.0 * a * a * a * c
                         - a * a * a * a;
    double part_Ib = 0.0;
    if( sgnA == 1 )
    {
      part_Ib = 4.0 * ( ( trA - a ) * ( 3.0 * a * sum[0] - 3.0 * a * a * c + a * a * a ) + a * sum[2] );
    }
    const double part_Ic = 3.0 * F * ( 2.0 * trA2 - F );

    delta1 = n * part_Ia + part_Ib + part_Ic;

    /*******************************/
    /**********  PART II  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES B != A */
    for(size_t positionB = 0; positionB < num_m; ++positionB )
    {

      if( positionB == positionA ) continue;

      const int sAB = sigmaAB[ positionA * num_m + positionB ];
      const int sgnB = positionB < num_h ? 1 : -1;
      const size_t offB = positionB * n * n;

      const double trB  = tr1( Matrices, num_h, num_l, n, positionB );
      const double trB2 = tr2( Matrices, num_h, num_l, n, positionB, positionB );
      const double trAB = tr2( Matrices, num_h, num_l, n, positionA, positionB );
      double trAB2 = 0.0;
      if( sgnA == 1 )
      {
        trAB2 = tr3( Matrices, num_h, num_l, n, positionA, positionB, positionB );
      }

      const double c_B = creal( Matrices[pos_upper+offB] );

      sum[0] = tr3_real_ij(Matrices, positionA, positionB, positionB, pos_x, pos_x);
      sum[1] = row_norm_squared(Matrices, positionB, pos_x);
      sum[2] = tr3_real_ij(Matrices, positionB, positionA, positionB, pos_x, pos_x);
      sum[3] = tr2_real_ij(Matrices, positionA, positionB, pos_x, pos_x);

      const double part_IIa = 2*a*sum[0] + a*a*sum[1];
      const double part_IIb = 2*a*sum[2] + a*a*c_B*c_B;
      double part_IIc;
      if( sgnA == 1 )
      {
        part_IIc = 2 * ( a*trAB2 + (a*trA+a*a)*sum[1] ) +
                   2 * ( a*a*c_B + 2*a*sum[3] ) * trB +
                   2 * sgnA * sgnB * ( 2*a*c_B*trAB + a*a*c_B*c_B ) +
                   trB2 * ( 2*a*c + a*a );
      }
      else
      {
        part_IIc = 2 * ( a*a*c_B + 2*a*sum[3] ) * trB +
		               2 * sgnA * sgnB * ( 2*a*c_B*trAB + a*a*c_B*c_B ) +
                   trB2 * ( 2*a*c + a*a );
      }

      delta2 += 4*n*part_IIa + 2*sAB*n*part_IIb + (4+2*sAB)*part_IIc;

    } /* END LOOP OVER MATRICES B */

    /*******************************/
    /*********  PART III  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES D != C != B != A pairwise, with   *
     * non-vanishing (pre-calculated) traces Tr(A*B*C*D) != 0   */
    for(int element=0; element<num_trABCD; element++) {

      const size_t positionB = sigmaABCD[positionA][9*element+1];
      const size_t positionC = sigmaABCD[positionA][9*element+2];
      const size_t positionD = sigmaABCD[positionA][9*element+3];
      const int b1           = sigmaABCD[positionA][9*element+4];
      const int b2           = sigmaABCD[positionA][9*element+5];
      const int c1           = sigmaABCD[positionA][9*element+6];
      const int c2           = sigmaABCD[positionA][9*element+7];
      const int d1           = sigmaABCD[positionA][9*element+8];
      const int d2           = sigmaABCD[positionA][9*element+9];
      const int s = b1 + b2 + c1 + c2 + d1 + d2;

      const size_t offB = positionB * n * n;
      const size_t offC = positionC * n * n;
      const size_t offD = positionD * n * n;

      const int sgnB = positionB < num_h ? 1 : -1;
      const int sgnC = positionC < num_h ? 1 : -1;
      const int sgnD = positionD < num_h ? 1 : -1;

      const double c_B = creal( Matrices[pos_upper+offB] );
      const double c_C = creal( Matrices[pos_upper+offC] );
      const double c_D = creal( Matrices[pos_upper+offD] );

      /* Traces over one matrix: */
      const double trB = tr1(Matrices, num_h, num_l, n, positionB);
      const double trC = tr1(Matrices, num_h, num_l, n, positionC);
      const double trD = tr1(Matrices, num_h, num_l, n, positionD);

      /* Traces over two matrices: */
      const double trABtrCD = c_B * tr2(Matrices, num_h, num_l, n, positionC, positionD);
      const double trACtrBD = c_C * tr2(Matrices, num_h, num_l, n, positionB, positionD);
      const double trADtrBC = c_D * tr2(Matrices, num_h, num_l, n, positionB, positionC);

      /* Traces over one and three matrices: */
      double trF1 = 0.0;
      double trF2 = 0.0;
      double trF3 = 0.0;
      double trF4 = 0.0;
      if( sgnA == 1 )
      {
        trF1 = (b1+c2+d1)*tr3(Matrices, num_h, num_l, n, positionB, positionC, positionD);
        trF1+= (b2+c1+d2)*tr3(Matrices, num_h, num_l, n, positionC, positionB, positionD);
      }
      if( sgnB == 1 )
      {
        trF2 = tr2_real_ij(Matrices, positionC, positionD, pos_x, pos_x);
      }
      if( sgnC == 1 )
      {
        trF3 = tr2_real_ij(Matrices, positionB, positionD, pos_x, pos_x);
      }
      if( sgnD == 1 )
      {
        trF4 = tr2_real_ij(Matrices, positionB, positionC, pos_x, pos_x);
      }

      /* Traces over four matrices: */
      double trAE = (b1+d2)*tr3_real_ij(Matrices, positionB, positionC, positionD, pos_x, pos_x);
      trAE+= (b2+c2)*tr3_real_ij(Matrices, positionB, positionD, positionC, pos_x, pos_x);
      trAE+= (c1+d1)*tr3_real_ij(Matrices, positionC, positionB, positionD, pos_x, pos_x);

#ifdef LARGE_N
      delta3 += a * n * trAE;
#else
      delta3 += a * ( n*trAE + trF1 + s*(trB*trF2 + trC*trF3 + trD*trF4) +
                      s*(sgnA*sgnB*trABtrCD + sgnA*sgnC*trACtrBD + sgnA*sgnD*trADtrBC) );
#endif
    } /* END LOOP OVER MATRICES */
  } /* END DIAGONAL CASE */

  return 2 * k * ( delta1 + delta2 + 4 * delta3 );
}

// Wrapper function for the full action:
// calculates S = g2 * Tr(D^2) + g4 * Tr(D^4)
double Calculate_Action(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD                      // pre-calculated Clifford products of 4 Gamma matrices
    )
{
  return prop.g2 * traceD2( Matrices, prop )
       + prop.g4 * traceD4( Matrices, prop, sigmaAB, sigmaABCD );
}

// Wrapper function for the change in the action:
// calculates delta_S = g2 * Tr(delta_D^2) + g4 * Tr(delta_D^4)
double Calculate_Delta_Action(
    REAL complex const *Matrices,        // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    int *sigmaAB,                        // pre-calculated Clifford products of 2 Gamma matrices
    int **sigmaABCD,                     // pre-calculated Clifford products of 4 Gamma matrices
    const struct Matrix_State old        // struct containing information about old state
    )
{
  return prop.g2 * delta_traceD2( Matrices, prop, old )
       + prop.g4 * delta_traceD4( Matrices, prop, old, sigmaAB, sigmaABCD );
}
