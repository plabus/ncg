#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_lapacke.h"
#include "Actions.h"
#include "Constants.h"

/**********************************************************************************************/

float tr1(MKL_Complex8 *Matrices, int pos1, int NUM_H) {
  float trace = 0.f;
  int off1 = pos1*SWEEP;

  if(pos1<NUM_H) {
    for(int a=0;a<N;++a) {
      trace += Matrices[a*N+a+off1].real;
    }
  }

  return trace;
}


float tr2(MKL_Complex8 *Matrices, int pos1, int pos2) {
  float trace = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;

  for(int a=0;a<N;a++) {
    for(int b=0;b<N;b++) {
      trace += Matrices[a*N+b+off1].real * Matrices[b*N+a+off2].real;
      trace -= Matrices[a*N+b+off1].imag * Matrices[b*N+a+off2].imag;
    }  
  }

  return trace;
}


float tr3(MKL_Complex8 *Matrices, int pos1, int pos2, int pos3) {
  float trace = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;
  int off3 = pos3*SWEEP;

  for(int a=0;a<N;a++) {
    for(int b=0;b<N;b++) {
      for(int c=0;c<N;c++) {
        trace += Matrices[a*N+b+off1].real * Matrices[b*N+c+off2].real * Matrices[c*N+a+off3].real;
        trace -= Matrices[a*N+b+off1].imag * Matrices[b*N+c+off2].imag * Matrices[c*N+a+off3].real;
        trace -= Matrices[a*N+b+off1].real * Matrices[b*N+c+off2].imag * Matrices[c*N+a+off3].imag;
        trace -= Matrices[a*N+b+off1].imag * Matrices[b*N+c+off2].real * Matrices[c*N+a+off3].imag;
      }
    }  
  }

  return trace;
}

/**********************************************************************************************/

float tr2_real_ij(MKL_Complex8 *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
  float sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    sum += Matrices[i*N+b+off1].real * Matrices[b*N+j+off2].real;
    sum -= Matrices[i*N+b+off1].imag * Matrices[b*N+j+off2].imag;
  }

  return sum;
}


float tr2_imag_ij(MKL_Complex8 *Matrices, int pos1, int pos2, int pos_x, int pos_y) {
  float sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    sum += Matrices[i*N+b+off1].real * Matrices[b*N+j+off2].imag;
    sum += Matrices[i*N+b+off1].imag * Matrices[b*N+j+off2].real;
  }

  return sum;
}


float tr3_real_ij(MKL_Complex8 *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
  float sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;
  int off3 = pos3*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    for(int c=0;c<N;c++) {
      sum += Matrices[i*N+b+off1].real * Matrices[b*N+c+off2].real * Matrices[c*N+j+off3].real;
      sum -= Matrices[i*N+b+off1].imag * Matrices[b*N+c+off2].imag * Matrices[c*N+j+off3].real;
      sum -= Matrices[i*N+b+off1].real * Matrices[b*N+c+off2].imag * Matrices[c*N+j+off3].imag;
      sum -= Matrices[i*N+b+off1].imag * Matrices[b*N+c+off2].real * Matrices[c*N+j+off3].imag;
    }
  }  

  return sum;
}


float tr3_imag_ij(MKL_Complex8 *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y) {
  float sum = 0.f;
  int off1 = pos1*SWEEP;
  int off2 = pos2*SWEEP;
  int off3 = pos3*SWEEP;

  /* We always assume i < j */
  int i = pos_x<=pos_y ? pos_x : pos_y;
  int j = pos_x<=pos_y ? pos_y : pos_x;

  for(int b=0;b<N;b++) {
    for(int c=0;c<N;c++) {
      sum -= Matrices[i*N+b+off1].imag * Matrices[b*N+c+off2].imag * Matrices[c*N+j+off3].imag;
      sum += Matrices[i*N+b+off1].real * Matrices[b*N+c+off2].real * Matrices[c*N+j+off3].imag;
      sum += Matrices[i*N+b+off1].real * Matrices[b*N+c+off2].imag * Matrices[c*N+j+off3].real;
      sum += Matrices[i*N+b+off1].imag * Matrices[b*N+c+off2].real * Matrices[c*N+j+off3].real;
    }
  }  

  return sum;
}


float matrix_norm_squared(MKL_Complex8 *Matrices, int pos1) {
  float sum_off_diag = 0.f;
  float trace_squared = 0.f;
  int off1 = pos1*SWEEP;

  for(int a=0;a<N;a++) {
    for(int b=a+1;b<N;b++) {
      sum_off_diag += Matrices[a*N+b+off1].real * Matrices[a*N+b+off1].real;
      sum_off_diag += Matrices[a*N+b+off1].imag * Matrices[a*N+b+off1].imag;
    }
  }  

  for(int a=0;a<N;a++) trace_squared += Matrices[a*N+a+off1].real * Matrices[a*N+a+off1].real;

  return 2*sum_off_diag+trace_squared;
}


float row_norm_squared(MKL_Complex8 *Matrices, int pos1, int pos_row) {
  float sum = 0.f;
  int off = pos1*SWEEP + pos_row*N;

  for(int a=0;a<N;a++) {
    sum += Matrices[off+a].real * Matrices[off+a].real;
    sum += Matrices[off+a].imag * Matrices[off+a].imag;
  }

  return sum;
}

/**********************************************************************************************/

float traceD2(MKL_Complex8 *Matrices, int NUM_H, int NUM_L)
{
  /* Calculating the action S = Tr D^2 */
  int itimesN;
  int Nplus1 = N+1;
  int offset;
  float trace_H2 = 0.f;
  float trace_H = 0.f; 

  /* For the Tr H^2 part */
  for(int n=0;n<NUM_M;++n) { /* This is for all matrices H and L */
    offset = n*SWEEP;
    for (int i=0;i<N;++i) {
      itimesN = i*N;
      for (int j=0;j<N;++j) {
        trace_H2 += Matrices[itimesN+j+offset].real*Matrices[itimesN+j+offset].real;
        trace_H2 += Matrices[itimesN+j+offset].imag*Matrices[itimesN+j+offset].imag;
      }
    }
  }

  /* For the (Tr H)^2 part */
  for(int n=0;n<NUM_H;++n) { /* This is for matrices H only */
    offset = n*SWEEP;
    for (int i=0;i<N;++i) {
      trace_H += Matrices[i*Nplus1+offset].real;
    }
  }

  return 2*K*(N*trace_H2+trace_H*trace_H);
}


float traceD4(MKL_Complex8 *Matrices, int NUM_H, int NUM_L)
{
  return 0;
}

/**********************************************************************************************/

/*********************************************************************************/
/********                                                                *********/
/********                      ACTION S = Tr D^2                         *********/
/********                                                                *********/
/*********************************************************************************/

float delta_action_traceD2(MKL_Complex8 *Matrices, int position, MKL_Complex8 temp, int pos_x, int pos_y, int NUM_H, int NUM_L)
{

  float delta;
  float trace_new, trace_old;
  int pos_upper = pos_x<=pos_y ? pos_x*N+pos_y : pos_y*N+pos_x;
  int offset = position*SWEEP;

  /* Generate a real number for diagonal element or else a complex number in upper half */
  if (pos_x == pos_y) { /* DIAGONAL CASE */

    if(position<NUM_H) { /* Case of the matrices H */
      trace_old = 0.f;
      for(int l=0;l<N;++l) trace_old += Matrices[l*N+l+offset].real;
      trace_new = trace_old - Matrices[pos_upper+offset].real + temp.real;

#ifdef LARGE_N
      delta = 2*K*N*(temp.real*temp.real - Matrices[pos_upper+offset].real*Matrices[pos_upper+offset].real);
      //printf("Yeah, yeah, I'm large N!\n");
#else
      delta = 2*K*N*(temp.real*temp.real - Matrices[pos_upper+offset].real*Matrices[pos_upper+offset].real) +
              2*K*(trace_new*trace_new - trace_old*trace_old);
#endif
    } else { /* Case of the traceless matrices L */
      delta = 2*K*N*(temp.real*temp.real - Matrices[pos_upper+offset].real*Matrices[pos_upper+offset].real);
    }

  } else { /* OFF-DIAGONAL CASE */
    /* If an off-diagonal element is changed the trace remains the same, *
     * which implies that the two cases for H and L are the same.        */
    delta = 4*K*N*(temp.real*temp.real + temp.imag*temp.imag -  
        Matrices[pos_upper+offset].real*Matrices[pos_upper+offset].real - Matrices[pos_upper+offset].imag*Matrices[pos_upper+offset].imag);
  }

  return delta;
}


/*********************************************************************************/
/********                                                                *********/
/********                      ACTION S = Tr D^4                         *********/
/********                                                                *********/
/*********************************************************************************/

float delta_action_traceD4(MKL_Complex8 *Matrices, int positionA, MKL_Complex8 temp, int pos_x, int pos_y, int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L)
{
  int off = positionA * SWEEP;
  int offB, offC, offD;
  int pos_upper = pos_x<=pos_y ? pos_x*N+pos_y : pos_y*N+pos_x;
  int pos_lower = pos_x<=pos_y ? pos_y*N+pos_x : pos_x*N+pos_y;
  int sgnA = positionA < NUM_H ? 1 : -1;
  int sgnB, sgnC, sgnD, sAB;
  int b1, b2, c1, c2, d1, d2, s;
  int num_trABCD = sigmaABCD[positionA][0];

  /* Save the old value of the matrix where a new element was generated. */
  MKL_Complex8 old = Matrices[pos_upper+off];

  float delta1;
  float delta2;
  float delta3;

  float trA;
  float trB;
  float trC;
  float trD;
  float trA2;
  float trB2;
  float trAB;
  float trAB2;
  
  float trABtrCD;
  float trACtrBD;
  float trADtrBC;
  float trF1;
  float trF2;
  float trF3;
  float trF4;
  float trAE;

  float a, b, c, d;
  float abs2, mixsum, mixdiff;
  float A_ii, A_jj;
  float c_B, d_B;
  float c_C, d_C;
  float c_D, d_D;
  float B_ii, B_jj;
  float F;
  float sum[10] = {0,0,0,0,0,0,0,0,0,0};
  float part_Ia, part_Ib, part_Ic;
  float part_IIa, part_IIb, part_IIc;


  /*******************************************************
   *                                                     *
   *                         INIT                        *
   *                                                     *
   *******************************************************/

  a = temp.real-old.real;
  b = temp.imag-old.imag;
  c = old.real;
  d = old.imag;

  trA = tr1(Matrices, positionA, NUM_H);
  trA2 = tr2(Matrices, positionA, positionA);
  A_ii = Matrices[pos_x*N+pos_x+off].real;
  A_jj = Matrices[pos_y*N+pos_y+off].real;

  abs2 = a*a + b*b;
  mixsum = a*c + b*d;
  mixdiff = a*d - b*c;

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

      trB  = tr1(Matrices, positionB, NUM_H);
      trB2 = tr2(Matrices, positionB, positionB);
      trAB = tr2(Matrices, positionA, positionB);

      c_B  = Matrices[pos_upper+offB].real;
      d_B  = Matrices[pos_upper+offB].imag;
      B_ii = Matrices[pos_x*N+pos_x+offB].real;
      B_jj = Matrices[pos_y*N+pos_y+offB].real;

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

      c_B = Matrices[pos_upper+offB].real;
      d_B = Matrices[pos_upper+offB].imag;
      c_C = Matrices[pos_upper+offC].real;
      d_C = Matrices[pos_upper+offC].imag;
      c_D = Matrices[pos_upper+offD].real;
      d_D = Matrices[pos_upper+offD].imag;

      /* Traces over one matrix: */
      trB = tr1(Matrices, positionB, NUM_H);
      trC = tr1(Matrices, positionC, NUM_H);
      trD = tr1(Matrices, positionD, NUM_H);

      /* Traces over two matrices: */
      trABtrCD = (a*c_B+b*d_B)*tr2(Matrices, positionC, positionD);
      trACtrBD = (a*c_C+b*d_C)*tr2(Matrices, positionB, positionD);
      trADtrBC = (a*c_D+b*d_D)*tr2(Matrices, positionB, positionC);

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
    if(sgnA==1) sum[2] = tr3(Matrices, positionA, positionA, positionA);
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

      trB   = tr1(Matrices, positionB, NUM_H);
      trB2  = tr2(Matrices, positionB, positionB);
      trAB  = tr2(Matrices, positionA, positionB);
      if(sgnA==1) trAB2 = tr3(Matrices, positionA, positionB, positionB);
      else        trAB2 = 0;

      c_B = Matrices[pos_upper+offB].real;

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

      c_B = Matrices[pos_upper+offB].real;
      c_C = Matrices[pos_upper+offC].real;
      c_D = Matrices[pos_upper+offD].real;

      /* Traces over one matrix: */
      trB = tr1(Matrices, positionB, NUM_H);
      trC = tr1(Matrices, positionC, NUM_H);
      trD = tr1(Matrices, positionD, NUM_H);

      /* Traces over two matrices: */
      trABtrCD = c_B * tr2(Matrices, positionC, positionD);
      trACtrBD = c_C * tr2(Matrices, positionB, positionD);
      trADtrBC = c_D * tr2(Matrices, positionB, positionC);

      /* Traces over one and three matrices: */
      if(sgnA==1) {
        trF1 = (b1+c2+d1)*tr3(Matrices, positionB, positionC, positionD);
        trF1+= (b2+c1+d2)*tr3(Matrices, positionC, positionB, positionD);
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
