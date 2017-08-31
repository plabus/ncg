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

double tr2_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
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


double tr2_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
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


double tr3_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
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


double tr3_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
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


double matrix_norm_squared(REAL complex *Matrices, int pos1) {
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


double row_norm_squared(REAL complex *Matrices, int pos1, int pos_row) {
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
double deltaS_traceD2(
    REAL complex *Matrices,              // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    const struct Matrix_State old        // old state
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
    REAL complex const *Matrices,       // Array of all matrices
    struct Matrix_Properties const prop // includes num_h, num_l, n and k
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
  size_t const n     = prop.n;
  size_t const k     = prop.k;

  double partI = 0.0;

  for( size_t posA = 0; posA < num_h + num_l; ++posA )
  {
    double const TrA1 = tr1( Matrices, num_h, num_l, n, posA );
    double const TrA2 = tr2( Matrices, num_h, num_l, n, posA, posA );
    double const TrA3 = tr3( Matrices, num_h, num_l, n, posA, posA, posA );
    double const TrA4 = tr4( Matrices, num_h, num_l, n, posA, posA, posA, posA );

    partI += n * TrA4 + 4 * TrA1 * TrA3 + 3 * TrA2 * TrA2;
  }

  partI *= 2 * k;

  double const action = partI;
  return action;
}

double delta_action_traceD4(
    REAL complex *Matrices,
    int positionA,
    REAL complex temp,
    int pos_x,
    int pos_y,
    int *sigmaAB,
    int **sigmaABCD,
    int NUM_H,
    int NUM_L
    )
{
  int off = positionA * SWEEP;
  int offB, offC, offD;
  int pos_upper = pos_x<=pos_y ? pos_x*N+pos_y : pos_y*N+pos_x;
  int sgnA = positionA < NUM_H ? 1 : -1;
  int sgnB, sgnC, sgnD, sAB;
  int b1, b2, c1, c2, d1, d2, s;
  int num_trABCD = sigmaABCD[positionA][0];

  /* Save the old value of the matrix where a new element was generated. */
  REAL complex old = Matrices[pos_upper+off];

  double delta1;
  double delta2;
  double delta3;

  double trA;
  double trB;
  double trC;
  double trD;
  double trA2;
  double trB2;
  double trAB;
  double trAB2;

  double trABtrCD;
  double trACtrBD;
  double trADtrBC;
  double trF1;
  double trF2;
  double trF3;
  double trF4;
  double trAE;

  double a, b, c, d;
  double abs2, mixsum;
  double A_ii, A_jj;
  double c_B, d_B;
  double c_C, d_C;
  double c_D, d_D;
  double B_ii, B_jj;
  double F;
  double sum[10] = {0,0,0,0,0,0,0,0,0,0};
  double part_Ia, part_Ib, part_Ic;
  double part_IIa, part_IIb, part_IIc;


  /*******************************************************
   *                                                     *
   *                         INIT                        *
   *                                                     *
   *******************************************************/

  a = creal(temp) - creal(old);
  b = cimag(temp) - cimag(old);
  c = creal(old);
  d = cimag(old);

  trA  = tr1(Matrices, NUM_H, NUM_L, N, positionA);
  trA2 = tr2(Matrices, NUM_H, NUM_L, N, positionA, positionA);
  A_ii = creal( Matrices[pos_x*N+pos_x+off] );
  A_jj = creal( Matrices[pos_y*N+pos_y+off] );

  abs2    = a * a + b * b;
  mixsum  = a * c + b * d;

  delta1 = 0.f;
  delta2 = 0.f;
  delta3 = 0.f;

  /*******************************************************
   *                                                     *
   *                    OFF-DIAGONAL                     *
   *                                                     *
   *******************************************************/

  if( pos_x != pos_y ) {

    /*******************************/
    /**********  PART I  ***********/
    /*******************************/

    sum[0] = tr3_real_ij(Matrices, positionA, positionA, positionA, pos_x, pos_y);
    sum[1] = tr3_imag_ij(Matrices, positionA, positionA, positionA, pos_x, pos_y);
    sum[2] = row_norm_squared(Matrices, positionA, pos_x);
    sum[3] = row_norm_squared(Matrices, positionA, pos_y);
    if(sgnA==1) {
      sum[4] = tr2_real_ij(Matrices, positionA, positionA, pos_x, pos_y);
      sum[5] = tr2_imag_ij(Matrices, positionA, positionA, pos_x, pos_y);
    } else {
      sum[4] = 0;
      sum[5] = 0;
    }

    F = 4*mixsum + 2*abs2;

    part_Ia = 8 * ( a*sum[0] + b*sum[1] ) +
              4 * abs2 * ( sum[2] + sum[3] ) +
              4 * ( (a*a-b*b)*(c*c-d*d) + 4*a*b*c*d + abs2*A_ii*A_jj ) +
              8 * abs2 * mixsum +
              2 * abs2 * abs2;
    part_Ib = 3 * trA * ( 2*a*sum[4] + 2*b*sum[5] + abs2*(A_ii+A_jj) );
    part_Ic = F * ( F + 2*trA2 );

#ifdef LARGE_N
    delta1 = N*part_Ia;
#else
    delta1 = N*part_Ia + 4*part_Ib + 3*part_Ic;
#endif

    /*******************************/
    /**********  PART II  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES B != A */
    for(int positionB=0; positionB<NUM_M; positionB++) {
      if(positionB==positionA) continue;

      sAB = sigmaAB[positionA*NUM_M+positionB];
      sgnB = positionB < NUM_H ? 1 : -1;
      offB = positionB * SWEEP;

      trB  = tr1(Matrices, NUM_H, NUM_L, N, positionB);
      trB2 = tr2(Matrices, NUM_H, NUM_L, N, positionB, positionB);
      trAB = tr2(Matrices, NUM_H, NUM_L, N, positionA, positionB);

      c_B  = creal( Matrices[pos_upper+offB] );
      d_B  = cimag( Matrices[pos_upper+offB] );
      B_ii = creal( Matrices[pos_x*N+pos_x+offB] );
      B_jj = creal( Matrices[pos_y*N+pos_y+offB] );

      sum[0] = tr3_real_ij(Matrices, positionA, positionB, positionB, pos_x, pos_y);
      sum[0]+= tr3_real_ij(Matrices, positionB, positionB, positionA, pos_x, pos_y);
      sum[1] = tr3_imag_ij(Matrices, positionA, positionB, positionB, pos_x, pos_y);
      sum[1]+= tr3_imag_ij(Matrices, positionB, positionB, positionA, pos_x, pos_y);
      sum[2] = row_norm_squared(Matrices, positionB, pos_x);
      sum[3] = row_norm_squared(Matrices, positionB, pos_y);
      sum[4] = tr3_real_ij(Matrices, positionB, positionA, positionB, pos_x, pos_y);
      sum[5] = tr3_imag_ij(Matrices, positionB, positionA, positionB, pos_x, pos_y);
      if(sgnA==1) {
        sum[6] = tr2_real_ij(Matrices, positionB, positionB, pos_x, pos_y);
        sum[7] = tr2_imag_ij(Matrices, positionB, positionB, pos_x, pos_y);
      } else {
        sum[6] = 0;
        sum[7] = 0;
      }
      if(sgnB==1) {
        sum[8] = tr2_real_ij(Matrices, positionA, positionB, pos_x, pos_y);
        sum[8]+= tr2_real_ij(Matrices, positionB, positionA, pos_x, pos_y);
        sum[9] = tr2_imag_ij(Matrices, positionA, positionB, pos_x, pos_y);
        sum[9]+= tr2_imag_ij(Matrices, positionB, positionA, pos_x, pos_y);
      } else {
        sum[8] = 0;
        sum[9] = 0;
      }

      part_IIa = 2*a*sum[0] + 2*b*sum[1] + abs2*(sum[2]+sum[3]);
      part_IIb = 4*a*sum[4] + 4*b*sum[5] +
                 2 * ( (a*a-b*b)*(c_B*c_B-d_B*d_B) +
	         4*a*b*c_B*d_B + abs2*B_ii*B_jj );
      part_IIc = 2 * trA * ( a*sum[6] + b*sum[7] ) +
                 trB * ( 2*a*sum[8] + 2*b*sum[9] + abs2*(B_ii+B_jj) ) +
                 4 * sgnA * sgnB * (a*c_B + b*d_B) * (trAB + a*c_B + b*d_B) +
                 trB2 * ( 2 * mixsum + abs2 );

#ifdef LARGE_N
      delta2 += 4*N*part_IIa + 2*sAB*N*part_IIb;
#else
      delta2 += 4*N*part_IIa + 2*sAB*N*part_IIb + (8+4*sAB)*part_IIc;
#endif
    } /* END LOOP OVER MATRICES B */

    /*******************************/
    /*********  PART III  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES D != C != B != A pairwise, with   *
     * non-vanishing (pre-calculated) traces Tr(A*B*C*D) != 0   */
    for(int element=0; element<num_trABCD; element++) {

      int positionB = sigmaABCD[positionA][9*element+1];
      int positionC = sigmaABCD[positionA][9*element+2];
      int positionD = sigmaABCD[positionA][9*element+3];
      b1            = sigmaABCD[positionA][9*element+4];
      b2            = sigmaABCD[positionA][9*element+5];
      c1            = sigmaABCD[positionA][9*element+6];
      c2            = sigmaABCD[positionA][9*element+7];
      d1            = sigmaABCD[positionA][9*element+8];
      d2            = sigmaABCD[positionA][9*element+9];
      s = b1 + b2 + c1 + c2 + d1 + d2;

      offB = positionB * SWEEP;
      offC = positionC * SWEEP;
      offD = positionD * SWEEP;

      sgnB = positionB < NUM_H ? 1 : -1;
      sgnC = positionC < NUM_H ? 1 : -1;
      sgnD = positionD < NUM_H ? 1 : -1;

      c_B = creal( Matrices[pos_upper+offB] );
      d_B = cimag( Matrices[pos_upper+offB] );
      c_C = creal( Matrices[pos_upper+offC] );
      d_C = cimag( Matrices[pos_upper+offC] );
      c_D = creal( Matrices[pos_upper+offD] );
      d_D = cimag( Matrices[pos_upper+offD] );

      /* Traces over one matrix: */
      trB = tr1(Matrices, NUM_H, NUM_L, N, positionB);
      trC = tr1(Matrices, NUM_H, NUM_L, N, positionC);
      trD = tr1(Matrices, NUM_H, NUM_L, N, positionD);

      /* Traces over two matrices: */
      trABtrCD = (a * c_B + b * d_B) * tr2(Matrices, NUM_H, NUM_L, N, positionC, positionD);
      trACtrBD = (a * c_C + b * d_C) * tr2(Matrices, NUM_H, NUM_L, N, positionB, positionD);
      trADtrBC = (a * c_D + b * d_D) * tr2(Matrices, NUM_H, NUM_L, N, positionB, positionC);

      /* Traces over one and three matrices: */
      if(sgnB==1) {
        trF2 = a*tr2_real_ij(Matrices, positionC, positionD, pos_x, pos_y);
        trF2+= a*tr2_real_ij(Matrices, positionD, positionC, pos_x, pos_y);
        trF2+= b*tr2_imag_ij(Matrices, positionC, positionD, pos_x, pos_y);
        trF2+= b*tr2_imag_ij(Matrices, positionD, positionC, pos_x, pos_y);
      } else {
        trF2 = 0;
      }

      if(sgnC==1) {
        trF3 = a*tr2_real_ij(Matrices, positionB, positionD, pos_x, pos_y);
        trF3+= a*tr2_real_ij(Matrices, positionD, positionB, pos_x, pos_y);
        trF3+= b*tr2_imag_ij(Matrices, positionB, positionD, pos_x, pos_y);
        trF3+= b*tr2_imag_ij(Matrices, positionD, positionB, pos_x, pos_y);
      } else {
        trF3 = 0;
      }

      if(sgnD==1) {
        trF4 = a*tr2_real_ij(Matrices, positionB, positionC, pos_x, pos_y);
        trF4+= a*tr2_real_ij(Matrices, positionC, positionB, pos_x, pos_y);
        trF4+= b*tr2_imag_ij(Matrices, positionB, positionC, pos_x, pos_y);
        trF4+= b*tr2_imag_ij(Matrices, positionC, positionB, pos_x, pos_y);
      } else {
        trF4 = 0;
      }

      /* Traces over four matrices: */
      trAE = (b1+d2)*a*tr3_real_ij(Matrices, positionB, positionC, positionD, pos_x, pos_y);
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
      delta3 +=   N*trAE;
#else
      delta3 += ( N*trAE + s*(trB*trF2 + trC*trF3 + trD*trF4) +
                  2*s*(sgnA*sgnB*trABtrCD + sgnA*sgnC*trACtrBD + sgnA*sgnD*trADtrBC) );
#endif
    } /* END LOOP OVER MATRICES */
  } /* END OFF-DIAGONAL CASE */

  /*******************************************************
   *                                                     *
   *                       DIAGONAL                      *
   *                                                     *
   *******************************************************/

  else {

    /*******************************/
    /**********  PART I  ***********/
    /*******************************/

    sum[0] = row_norm_squared(Matrices, positionA, pos_x);
    sum[1] = tr3_real_ij(Matrices, positionA, positionA, positionA, pos_x, pos_x);
    if(sgnA==1) sum[2] = tr3(Matrices, NUM_H, NUM_L, N, positionA, positionA, positionA);
    else        sum[2] = 0;

    F = a*a + 2*a*c;

    part_Ia = a*a*a*a + 4*a*a*a*c + 2*a*a*c*c + 4*a*a*sum[0] + 4*a*sum[1];
    if(sgnA==1) part_Ib = a*sum[2] + ( a  + trA ) * ( a*a*a + 3*a*a*c + 3*a*sum[0] );
    else        part_Ib = 0;
    part_Ic = F * ( F + 2*trA2 );

#ifdef LARGE_N
    delta1 = N*part_Ia;
#else
    delta1 = N*part_Ia + 4*part_Ib + 3*part_Ic;
#endif

    /*******************************/
    /**********  PART II  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES B != A */
    for(int positionB=0; positionB<NUM_M; positionB++) {

      if(positionB==positionA) continue;

      sAB = sigmaAB[positionA*NUM_M+positionB];
      sgnB = positionB < NUM_H ? 1 : -1;
      offB = positionB * SWEEP;

      trB  = tr1(Matrices, NUM_H, NUM_L, N, positionB);
      trB2 = tr2(Matrices, NUM_H, NUM_L, N, positionB, positionB);
      trAB = tr2(Matrices, NUM_H, NUM_L, N, positionA, positionB);
      if(sgnA==1) trAB2 = tr3(Matrices, NUM_H, NUM_L, N, positionA, positionB, positionB);
      else        trAB2 = 0;

      c_B = creal( Matrices[pos_upper+offB] );

      sum[0] = tr3_real_ij(Matrices, positionA, positionB, positionB, pos_x, pos_x);
      sum[1] = row_norm_squared(Matrices, positionB, pos_x);
      sum[2] = tr3_real_ij(Matrices, positionB, positionA, positionB, pos_x, pos_x);
      sum[3] = tr2_real_ij(Matrices, positionA, positionB, pos_x, pos_x);

      part_IIa = 2*a*sum[0] + a*a*sum[1];
      part_IIb = 2*a*sum[2] + a*a*c_B*c_B;
      if(sgnA==1) {
        part_IIc = 2 * ( a*trAB2 + (a*trA+a*a)*sum[1] ) +
                   2 * ( a*a*c_B + 2*a*sum[3] ) * trB +
                   2 * sgnA * sgnB * ( 2*a*c_B*trAB + a*a*c_B*c_B ) +
                   trB2 * ( 2*a*c + a*a );
      } else {
        part_IIc = 2 * ( a*a*c_B + 2*a*sum[3] ) * trB +
		   2 * sgnA * sgnB * ( 2*a*c_B*trAB + a*a*c_B*c_B ) +
                   trB2 * ( 2*a*c + a*a );
      }

#ifdef LARGE_N
      delta2 += 4*N*part_IIa + 2*sAB*N*part_IIb;
#else
      delta2 += 4*N*part_IIa + 2*sAB*N*part_IIb + (4+2*sAB)*part_IIc;
#endif
    } /* END LOOP OVER MATRICES B */

    /*******************************/
    /*********  PART III  **********/
    /*******************************/

    /* LOOP OVER ALL MATRICES D != C != B != A pairwise, with   *
     * non-vanishing (pre-calculated) traces Tr(A*B*C*D) != 0   */
    for(int element=0; element<num_trABCD; element++) {

      int positionB = sigmaABCD[positionA][9*element+1];
      int positionC = sigmaABCD[positionA][9*element+2];
      int positionD = sigmaABCD[positionA][9*element+3];
      b1            = sigmaABCD[positionA][9*element+4];
      b2            = sigmaABCD[positionA][9*element+5];
      c1            = sigmaABCD[positionA][9*element+6];
      c2            = sigmaABCD[positionA][9*element+7];
      d1            = sigmaABCD[positionA][9*element+8];
      d2            = sigmaABCD[positionA][9*element+9];
      s = b1 + b2 + c1 + c2 + d1 + d2;

      offB = positionB * SWEEP;
      offC = positionC * SWEEP;
      offD = positionD * SWEEP;

      sgnB = positionB < NUM_H ? 1 : -1;
      sgnC = positionC < NUM_H ? 1 : -1;
      sgnD = positionD < NUM_H ? 1 : -1;

      c_B = creal( Matrices[pos_upper+offB] );
      c_C = creal( Matrices[pos_upper+offC] );
      c_D = creal( Matrices[pos_upper+offD] );

      /* Traces over one matrix: */
      trB = tr1(Matrices, NUM_H, NUM_L, N, positionB);
      trC = tr1(Matrices, NUM_H, NUM_L, N, positionC);
      trD = tr1(Matrices, NUM_H, NUM_L, N, positionD);

      /* Traces over two matrices: */
      trABtrCD = c_B * tr2(Matrices, NUM_H, NUM_L, N, positionC, positionD);
      trACtrBD = c_C * tr2(Matrices, NUM_H, NUM_L, N, positionB, positionD);
      trADtrBC = c_D * tr2(Matrices, NUM_H, NUM_L, N, positionB, positionC);

      /* Traces over one and three matrices: */
      if(sgnA==1) {
        trF1 = (b1+c2+d1)*tr3(Matrices, NUM_H, NUM_L, N, positionB, positionC, positionD);
        trF1+= (b2+c1+d2)*tr3(Matrices, NUM_H, NUM_L, N, positionC, positionB, positionD);
      } else {
        trF1 = 0;
      }

      if(sgnB==1) {
        trF2 = tr2_real_ij(Matrices, positionC, positionD, pos_x, pos_x);
      } else {
        trF2 = 0;
      }

      if(sgnC==1) {
        trF3 = tr2_real_ij(Matrices, positionB, positionD, pos_x, pos_x);
      } else {
        trF3 = 0;
      }

      if(sgnD==1) {
        trF4 = tr2_real_ij(Matrices, positionB, positionC, pos_x, pos_x);
      } else {
        trF4 = 0;
      }

      /* Traces over four matrices: */
      trAE = (b1+d2)*tr3_real_ij(Matrices, positionB, positionC, positionD, pos_x, pos_x);
      trAE+= (b2+c2)*tr3_real_ij(Matrices, positionB, positionD, positionC, pos_x, pos_x);
      trAE+= (c1+d1)*tr3_real_ij(Matrices, positionC, positionB, positionD, pos_x, pos_x);

#ifdef LARGE_N
      delta3 += a * N * trAE;
#else
      delta3 += a * ( N*trAE + trF1 + s*(trB*trF2 + trC*trF3 + trD*trF4) +
                      s*(sgnA*sgnB*trABtrCD + sgnA*sgnC*trACtrBD + sgnA*sgnD*trADtrBC) );
#endif
    } /* END LOOP OVER MATRICES */
  } /* END DIAGONAL CASE */

  return 2*K*(delta1+delta2+4*delta3);
}

// Wrapper function for the full action:
// calculates S = g2 * Tr(D^2) + g4 * Tr(D^4)
double Calculate_Action(
    REAL complex const *Matrices,       // array of matrices
    struct Matrix_Properties const prop // includes num_h, num_l, n and k
    )
{
  return prop.g2 * traceD2( Matrices, prop ) + prop.g4 * traceD4( Matrices, prop );
}
