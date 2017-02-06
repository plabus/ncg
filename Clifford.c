#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "Clifford.h"
#include "Constants.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define sgn(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

void printMint(int *matrix, int n) {
  for(int i=0;i<n;++i) {  
    for(int j=0;j<n;++j) {  
      if(matrix[i*n+j]==1) printf("  1");
      else if(matrix[i*n+j]==-1) printf(" -1");
      else if(matrix[i*n+j]==2) printf("  i");
      else if(matrix[i*n+j]==-2) printf(" -i");
      else printf("  0");
    }
    printf("\n");
  }
}

void printMcomplex(float complex *matrix, int n) {
  for(int i=0;i<n;++i) {  
    for(int j=0;j<n;++j) {  
      if( abs(creal(matrix[i*n+j]))<1e-3 && abs(cimag(matrix[i*n+j]))<1e-3 ) printf("  0");
      else if( abs(creal(matrix[i*n+j])-1.)<1e-3 && abs(cimag(matrix[i*n+j]))<1e-3 ) printf("  1");
      else if( abs(creal(matrix[i*n+j])+1.)<1e-3 && abs(cimag(matrix[i*n+j]))<1e-3 ) printf(" -1");
      else if( abs(creal(matrix[i*n+j]))<1e-3 && abs(cimag(matrix[i*n+j])-1.)<1e-3 ) printf("  i");
      else if( abs(creal(matrix[i*n+j]))<1e-3 && abs(cimag(matrix[i*n+j])+1.)<1e-3 ) printf(" -i");
      else printf("  %.2f%+.2fi", creal(matrix[i*n+j]), cimag(matrix[i*n+j]));
    }
    printf("\n");
  }
}

int factorial(int n) {
  if(n==0) return 1;
  else if(n==1) return 1;
  else return (n*factorial(n-1));
}

int binomial(int a, int b) {
  return factorial(a)/(factorial(b)*factorial(a-b));
}

void commutator(float complex *A, float complex *B, float complex *C, int k) {
  /* Calculate C = A*B-B*A */

  /* Initialise the matrix C to zero */
  for(int i=0;i<k*k;++i) C[i] = 0. + 0.*I;

  float complex *buffer = (float complex*) calloc(k*k,sizeof(float complex));

  /* A*B */
  for(int i=0;i<k;++i) {
    for(int j=0;j<k;++j) {
      for(int l=0;l<k;++l) {
        buffer[i*k+j] += A[i*k+l]*B[l*k+j];
      }
    }
  }

  /* B*A */
  for(int i=0;i<k;++i) {
    for(int j=0;j<k;++j) {
      for(int l=0;l<k;++l) {
        C[i*k+j] += B[i*k+l]*A[l*k+j];
      }
    }
  }

  /* Adding the contributions up *
   * C = buffer - C              */
  for(int i=0;i<k*k;++i) {
    C[i] = (buffer[i]-C[i])/2;
  }

  free(buffer);
}

/** [combination c n p x]
 *  * get the [x]th lexicographically ordered set of [p] elements in [n]
 *   * output is in [c], and should be sizeof(int)*[p] */
void combination(int* c,int n,int p, int x){
  int i,r,k = 0;
  for(i=0;i<p-1;i++){
    c[i] = (i != 0) ? c[i-1] : 0;
    do {
      c[i]++;
      r = binomial(n-c[i],p-(i+1));
      k = k + r;
    } while(k < x);
    k = k - r;
  }
  c[p-1] = c[p-2] + x - k;
}

void Count_Hs_and_Ls(int p, int q, int *seq, int n, int *num_h, int *num_l) {
  /* Assume n>=1 */
  int exp;

  if(n==1) {
    if(seq[0]<=p) *num_h += 1;
    else          *num_l += 1;
  }
  else {
    exp = (n-1)*n/2;
    for(int i=0;i<n;++i) if(seq[i]>p) exp++;

    if(exp%2==0) *num_h += 1;
    else         *num_l += 1;
  }

}

void antisymmetrise(float complex *gammas, int dim, int k, int num_indices, int *sequence, float complex *matrix) {

  int n1, n2;
  int index;

  /* Initialise the matrix to zero */
  for(int i=0;i<k*k;++i) matrix[i] = 0. + 0.*I;

  /* Check that there aren't more indices than dimension */
  if(num_indices>dim) {
    printf("  ERROR: Number of indices bigger than dimension!\n");
    return;
  }
  /* Check that there aren't any elements in sequence bigger than dimension */
  for(int i=0;i<num_indices;++i) {
    if(sequence[i]>dim) {
      printf("  ERROR: Sequence out of range!\n");
      return;
    }
  }

  /* If there aren't any indices set the matrix to unity */
  if(num_indices==0) {
    for(int i=0;i<k;++i) {
      matrix[i*k+i] = 1;
    }
  }
  /* If there is one index copy relevant gamma to matrix */
  else if(num_indices==1) {
    /* Note n1 in [1,2,...,dim] */
    n1 = sequence[0]-1;
    for(int i=0;i<k*k;++i) {
      matrix[i] = gammas[i+k*k*n1];
    }
  }
  /* If there are two indeces calculate commutator and save in matrix */
  else if(num_indices==2) {
    n1 = sequence[0]-1;
    n2 = sequence[1]-1;
    commutator(&gammas[n1*k*k], &gammas[n2*k*k], matrix, k);
  }
  /* If there are more than two indices anti-symmetrise recursively */
  else if(num_indices>2) {
    /* Iterate over all elements in sequence:               *
     * 0. choose element n                                  *
     * 1. seq_new = [a1, a2, ..., (an), ..., a_num_indices] *
     * 2. call antisymmetrise recursively and save in buff1 *
     * 3. take the product between buff1 and the remaining  *
     *    gamma matrix and save in buff2                    *
     * 4. matrix += (-1)^(pos-1) * buff2                    *
     * 5. after iteration matrix = matrix / num_indices     */
     int new_seq[num_indices-1];
     float complex *buff1 = (float complex*) malloc(k*k*sizeof(float complex));
     float complex *buff2 = (float complex*) malloc(k*k*sizeof(float complex));

     for(int i=0;i<num_indices;++i) { /* Iterator over sequence */
       /* 0. Initialise buff2 to zero */
       for(int ii=0;ii<k*k;++ii) buff2[ii] = 0. + 0.*I;
       /* 1. Copy new sequence w/o ith element */
       index = 0;
       for(int j=0;j<num_indices;++j) {
         if(i!=j) {
           new_seq[index] = sequence[j];
           index++;
         }
       }
       /* 2. Recursive step */
       antisymmetrise(gammas, dim, k, num_indices-1, new_seq, buff1);
       /* 3. Take the product */
       for(int ii=0;ii<k;++ii) {
         for(int jj=0;jj<k;++jj) {
           for(int ll=0;ll<k;++ll) {
             buff2[ii*k+jj] += buff1[ii*k+ll]*gammas[ll*k+jj + (sequence[i]-1)*k*k];
           }
         }
       }
       /* 4. Add to the result matrix */
       for(int ii=0;ii<k*k;++ii)
         if(i%2==0) matrix[ii] += buff2[ii];
         else       matrix[ii] -= buff2[ii];
     }
     /* 5. Divide by #indices */
     for(int i=0;i<k*k;++i) matrix[i] /= num_indices;

     free(buff1);
     free(buff2);
  }

}

void Generate_Gammas(int p, int q, float complex *gammas) {
  /* We generate gamma matrices for the case (p+q,0) first *
   * and in the end multiply the last q matrices with i.   */

  /* Initial Matrices:                             *
   * (1,0): gamma_1 = (1)                          *
   * (2,0): gamma_1 = sigma_1,                     *
   *        gamma_2 = sigma_2.                     *
   *                                               *
   * GAMMA MATRICES IN D-DIM (D even):             *
   * (gamma_mu x sigma1, 1 x sigma2, 1 x sigma3 )  *
   * where gamma_mu are the matrices for d = dim-2 * 
   *                                               *
   * GAMMA MATRICES IN D-DIM (D odd):              *
   * (gamma_mu, gamma_5)                           *
   * where gamma_mu are the matrices in d = dim-1  *
   *                                               */

  int dim = p+q;
  int s = (8*8*8-dim)%8;
  int k = (dim%2==0) ? (int)pow(2,dim/2) : (int)pow(2,(dim-1)/2);
  int size = dim*k*k; /* space for all gamma matrices */
  int small_k = k/2; /* for dimensional recursion */
  int small_size = (dim-2)*small_k*small_k;

  int offset, offset2;

  /* PAULI MATRICES */
  float complex sigma1[4] = {0,1,1,0};
  float complex sigma2[4] = {0,-I,I,0};
  float complex sigma3[4] = {1,0,0,-1};

  /* Initialising the gamma matrices recursively for type (d,0) */
  if(dim==1) {
    gammas[0] = 1;
  }
  else if(dim==2) {
    for(int i=0;i<k*k;++i) {
      gammas[i] = sigma1[i];
      gammas[i+k*k] = sigma2[i];
    }
  }

  /* The EVEN DIMENSIONAL case */
  else if(dim>2 && dim%2==0) { 
    /* Allocate memory for the (d-2) dimensional matrices */
    float complex *small_gammas = (float complex*) calloc(small_size,sizeof(float complex));

    /* First generate the gamma matrices for dim-2 */
    Generate_Gammas(dim-2, 0, small_gammas);

    /* gamma_mu (x) sigma1 */
    for(int n=0;n<dim-2;++n) { /* Loop over matrices */
      offset = n * small_k * small_k; /* From one small matrix to the next */
      offset2 = n * k * k; /* Offset big matrices */
      for(int i=0;i<small_k;++i) { /* Indices gamma matrices */
        for(int j=0;j<small_k;++j) {
          for(int ii=0;ii<2;++ii) { /* Indices sigma1 */
            for(int jj=0;jj<2;++jj) {
              gammas[i*4*small_k+2*j+(ii*2*small_k+jj)+offset2] = small_gammas[i*small_k+j+offset] * sigma1[ii*2+jj];
            }
          }
        }
      }
    }
    /* 1I_(k-2) (x) sigma2 */
    for(int i=0;i<small_k;++i) { /* Diagonal unity matrix */
      for(int ii=0;ii<2;++ii) { /* Indices sigma1 */
        for(int jj=0;jj<2;++jj) {
          gammas[i*4*small_k+2*i+(ii*2*small_k+jj)+((dim-2)*k*k)] = sigma2[ii*2+jj];
        }
      }
    }
    /* 1I_(k-2) (x) sigma3 */
    for(int i=0;i<small_k;++i) { /* Diagonal unity matrix */
      for(int ii=0;ii<2;++ii) { /* Indices sigma1 */
        for(int jj=0;jj<2;++jj) {
          gammas[i*4*small_k+2*i+(ii*2*small_k+jj)+((dim-1)*k*k)] = sigma3[ii*2+jj];
        }
      }
    }
    free(small_gammas);
  }

  /* The ODD DIMENSONAL case */
  else {

    /* First generate the gamma matrices for dim-1 *
     * which is even dimensional                   */
    Generate_Gammas(dim-1, 0, gammas);

    /* Buffer for matrix multiplication */
    float complex *buffer = (float complex*) calloc(k*k,sizeof(float complex));

    /* Add g_5 to the end of the array */
    offset = (dim-1)*k*k;

    /* Multiply gamma1 with gamma2 and save in buffer */
    for(int a=0;a<k;++a) {
      for(int b=0;b<k;++b) {
        for(int c=0;c<k;++c) {
          buffer[a*k+b] += gammas[a*k+c] * gammas[c*k+b+k*k];
        }
      }
    }

    /* Multiply the product of first (n-1) gamma matrices with the nth and save in gamma_5 */
    for(int n=2;n<dim-1;++n) {
      for(int a=0;a<k;++a) {
        for(int b=0;b<k;++b) {
          for(int c=0;c<k;++c) {
            gammas[a*k+b+offset] += buffer[a*k+c] * gammas[c*k+b+n*k*k];
          }
        }
      }
      /* Copy new gamma_5 back to buffer and set gamma_5 to zero */
      for(int i=0;i<k*k;++i) {
        buffer[i] = gammas[i+offset];
        gammas[i+offset] = 0; 
      }
    }

    /* Calculate prefactor i^s(s+1)/2 (this is still for type (d,0) */
    int exponent = (s*(s+1)/2)%4;
    float complex factor;
    if(exponent==0) factor= 1; /* i^4n */
    if(exponent==1) factor= I; /* i^4n+1 */
    if(exponent==2) factor=-1; /* i^4n+2 */
    if(exponent==3) factor=-I; /* i^4n+3 */

    /* Multipy with prefactor */
    for(int i=0;i<k*k;++i) {
      gammas[i+offset] = factor * buffer[i];
    }
    
    free(buffer);
  }

  /* Multiply the last q matrices with i */
  for(int i=p*k*k;i<size;++i)
    gammas[i] = I * gammas[i];
}

void Reshuffle_Clifford_Group(int p, int q, float complex *big_gammas, int num_h, int num_l, int ODD) {
  /* Order the big gammas in such a way, that the hermitian matrices come first */

  int dim = p+q;
  int k = (dim%2==0) ? (int)pow(2,dim/2) : (int)pow(2,(dim-1)/2);
  int num = (int)pow(2,dim);
  if(ODD==1) num = (int)pow(2,dim-1);
  float complex first;

  float complex *buffer = (float complex*) malloc(num*k*k*sizeof(float complex));

  int off1 = 0, off2 = num_h*k*k;
  
  for(int n=0;n<num;++n) {
    /* Claculate the first element of the gamma matrix squared */
    first = 0;
    for(int i=0;i<k;++i) {
      first += big_gammas[i+n*k*k] * big_gammas[i*k+n*k*k];
    }
    /* If it squares to 1, copy it to the beginning of the buffer */
    if(abs(creal(first)-1)<1e-5) {
      for(int i=0;i<k*k;++i) buffer[i+off1] = big_gammas[i+n*k*k];
      off1 += k*k;
    }
    /* If it squares to -1 it's anti-hermitian */
    else if(abs(creal(first)+1)<1e-5) {
      for(int i=0;i<k*k;++i) buffer[i+off2] = big_gammas[i+n*k*k];
      off2 += k*k;
    }
  }

  /* Copy the whole buffer back to big_gammas */
      for(int i=0;i<num*k*k;++i) big_gammas[i] = buffer[i];

  free(buffer);

}

void Generate_Clifford_Group(int p, int q, float complex *big_gammas, int *num_h, int *num_l) {

  int dim = p+q;
  int k = (dim%2==0) ? (int)pow(2,dim/2) : (int)pow(2,(dim-1)/2);
  int num_mat; /* Number of gammas with fixed number of indices */
  int offset = k*k; /* Offset to iterate over big gammas */

  /* Initialise number of (anit-)hermitian matrices to zero */
  *num_h = 1; /* The unity matrix 1I is hermitian */
  *num_l = 0;

  /* Allocate space for a k*k buffer matrix */
  float complex *matrix = (float complex*) calloc(k*k,sizeof(float complex));
  /* Allocate space for index sequence */
  int *sequence = (int*) malloc(dim*sizeof(int));
  /* Allocate space for gamma matrices and generate them */
  int size_gammas = dim*k*k; /* space for all gamma matrices */
  float complex *gammas = (float complex*) calloc(size_gammas,sizeof(float complex));
  Generate_Gammas(p,q,gammas);

  /* Now generate (anti-symmetric) products of 0 up to d
   * gamma matrices:
   * Gamma_M = {1I, gamma_mu, gamma_mu,nu, ..., gamma_5} */

  /* Write the unit matrix 1I first */
  for(int i=0;i<k;++i) big_gammas[i*k+i] = 1;

  for(int n=1;n<=dim;++n) { /* Number of indices in the gammas */
    num_mat = binomial(dim, n); /* Calculate number of matrices with fixed number of indices */
    for(int m=0;m<num_mat;++m) { /* Iterations over gammas with fixed number of indices */
      combination(sequence, dim, n, m+1);
      Count_Hs_and_Ls(p, q, sequence, n, num_h, num_l);
      antisymmetrise(gammas, dim, k, n, sequence, matrix);
      for(int i=0;i<k*k;++i) {
        big_gammas[i+offset] = matrix[i];
      }
      offset += k*k;
    }
  }

  Reshuffle_Clifford_Group(p, q, big_gammas, *num_h, *num_l, 0);

  free(matrix);
  free(sequence);
  free(gammas);

}

void Generate_Clifford_Odd_Group(int p, int q, float complex *big_gammas, int *num_h, int *num_l) {

  int dim = p+q;
  int k = (dim%2==0) ? (int)pow(2,dim/2) : (int)pow(2,(dim-1)/2);
  int num_mat; /* Number of gammas with fixed number of indices */
  int offset = 0; /* Offset to iterate over big gammas */

  /* Initialise number of (anit-)hermitian matrices to zero */
  *num_h = 0;
  *num_l = 0;

  /* (1) Allocate space for a k*k buffer matrix */
  /* (2) Allocate space for index sequence */
  /* (3) Allocate space for gamma matrices and generate them */
  float complex *matrix = (float complex*) calloc(k*k,sizeof(float complex));
  int *sequence = (int*) malloc(dim*sizeof(int));
  int size_gammas = dim*k*k; /* space for all gamma matrices */
  float complex *gammas = (float complex*) calloc(size_gammas,sizeof(float complex));
  Generate_Gammas(p,q,gammas);

  /* Now generate odd products gamma matrices:
   * Gamma_Odd__M = {gamma_mu, gamma_mu_nu_rho, ...} */

  for(int n=1;n<=dim;n+=2) { /* Number of indices in the gammas (odd) */
    num_mat = binomial(dim, n); /* Calculate number of matrices with fixed number of indices */
    for(int m=0;m<num_mat;++m) { /* Iterations over gammas with fixed number of indices */
      combination(sequence, dim, n, m+1);
      Count_Hs_and_Ls(p, q, sequence, n, num_h, num_l);
      antisymmetrise(gammas, dim, k, n, sequence, matrix);
      for(int i=0;i<k*k;++i) {
        big_gammas[i+offset] = matrix[i];
      }
      offset += k*k;
    }
  }

  /* Reshuffle the indices to have all hermitian matrices first */
  Reshuffle_Clifford_Group(p, q, big_gammas, *num_h, *num_l, 1);

  free(matrix);
  free(sequence);
  free(gammas);

}

void Calculate_Trace_Gamma_ABAB(float complex *Gamma_Matrices, int *sigmaAB, int NUM_H) {
  int offi, offj;
  int num_m = (int)pow(2,D-1);
  float complex trace;

  for(int i=0;i<num_m;++i) {
    for(int j=i+1;j<num_m;++j) {

      trace = 0.f;
      offi = i*K*K;
      offj = j*K*K;

      for(int ii=0;ii<K;++ii)
        for(int jj=0;jj<K;++jj)
          for(int ll=0;ll<K;++ll)
            for(int mm=0;mm<K;++mm)
              trace += Gamma_Matrices[ii*K+jj+offi] * Gamma_Matrices[jj*K+ll+offj] *
                Gamma_Matrices[ll*K+mm+offi] * Gamma_Matrices[mm*K+ii+offj];

      if( abs(creal(trace))>0 ) {
        trace *= (i<NUM_H?1:-1)*(j<NUM_H?1:-1);
        sigmaAB[i*num_m+j] = sgn(creal(trace));
        sigmaAB[j*num_m+i] = sgn(creal(trace));
      }

    }
  }

}

int calculate_sigmaABCD(float complex *Gamma_Matrices, int a, int b, int c, int d, int NUM_H) {
  float trace = 0.f;
  int offa=a*K*K;
  int offb=b*K*K;
  int offc=c*K*K;
  int offd=d*K*K;
  
  /* TRACE OVER FOUR MATRICES */
  for(int ii=0;ii<K;++ii)
    for(int jj=0;jj<K;++jj)
      for(int ll=0;ll<K;++ll)
        for(int mm=0;mm<K;++mm)
          trace += Gamma_Matrices[ii*K+jj+offa] * Gamma_Matrices[jj*K+ll+offb] *
                   Gamma_Matrices[ll*K+mm+offc] * Gamma_Matrices[mm*K+ii+offd];
  /* Include i's coming from L's and cast trace to int. *
   * Note that sgn = +1 implies that the result is real */
  trace *= (a<NUM_H?1:I)*(b<NUM_H?1:I)*(c<NUM_H?1:I)*(d<NUM_H?1:I);
  return sgn(creal(trace));
}

void Calculate_Trace_Gamma_ABCD(float complex *Gamma_Matrices, int **sigmaABCD, int NUM_H) {

  int counter;
  int num_m = pow(2,D-1);
  int total_sign;
  int s1, s2, c1, c2, c3;
  int comb1;
  int comb2;
  int comb3;
  int comb4;
  int comb5;
  int comb6;

  /* LOOPING OVER FOUR PAIRWISE DIFFERENT MATRICES */
  for(int a=0;a<num_m;++a) {
    counter = 0;

    for(int b=0;b<num_m;++b) {
      if(b==a) continue;

      for(int c=0;c<num_m;++c) {
        if(c<=b || c==a) continue;

        for(int d=0;d<num_m;++d) {
          if(d<=b || d<=c || d==a) continue;

          /* (Anti-) commutators of Hermitian matrices vanish if there aren't *
           * pairs of (tracefree) non-tracefree matrices <=> sgn = -1.        */ 
	  total_sign = (a<NUM_H?1:-1)*(b<NUM_H?1:-1)*(c<NUM_H?1:-1)*(d<NUM_H?1:-1);
          if(total_sign==-1) continue;
          comb1 = calculate_sigmaABCD(Gamma_Matrices, a, b, c, d, NUM_H);

          /* SET VALUES IF THERE ARE NON-ZERO */
          if(comb1!=0) {
            comb2 = calculate_sigmaABCD(Gamma_Matrices, a, b, d, c, NUM_H);
            comb3 = calculate_sigmaABCD(Gamma_Matrices, a, c, b, d, NUM_H);
            comb4 = calculate_sigmaABCD(Gamma_Matrices, a, c, d, b, NUM_H);
            comb5 = calculate_sigmaABCD(Gamma_Matrices, a, d, b, c, NUM_H);
            comb6 = calculate_sigmaABCD(Gamma_Matrices, a, d, c, b, NUM_H);

            /* write to sigma s.t. sigmaABCD[A][] = {#, B, C, D, b1, b2, c1, c2, d1, d2, ... }  */
            sigmaABCD[a][9*counter+1] = b;
            sigmaABCD[a][9*counter+2] = c;
            sigmaABCD[a][9*counter+3] = d;
            sigmaABCD[a][9*counter+4] = comb1;
            sigmaABCD[a][9*counter+5] = comb2;
            sigmaABCD[a][9*counter+6] = comb3;
            sigmaABCD[a][9*counter+7] = comb4;
            sigmaABCD[a][9*counter+8] = comb5;
            sigmaABCD[a][9*counter+9] = comb6;
            counter++;
          }

        }
      }
    }

    /* Set the number of quadruples (a,b,c,d) with  *
     * fixed a to the first element of sigmaABCD[a] */
    sigmaABCD[a][0] = counter;

  } /* END LOOP over a */

}
