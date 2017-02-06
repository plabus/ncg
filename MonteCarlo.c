#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/time.h>
#include "mkl_lapacke.h"
#include "MonteCarlo.h"
#include "Random.h"
#include "Actions.h"
#include "Constants.h"



// Returns time in seconds (double)
double cclock() {
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

// Print the Matrix (for debugging)
void printM(const int L, MKL_Complex8 * M) {
  for(int i = 0; i < L; i++) {
    for(int j = 0; j < L; j++)
      if(M[i*L+j].imag>=0) printf( " %.2f + %.2fi\t ", M[i*L+j].real, M[i*L+j].imag );
      else                 printf( " %.2f - %.2fi\t ", M[i*L+j].real, -M[i*L+j].imag );
    printf("\n");
  }
}

// Two compare two floats for sorting
int compare_floats(const void *a, const void *b) {
  const float *da = (const float *) a;
  const float *db = (const float *) b;
  return (*da > *db) - (*da < *db);
}

/* Delta approximation for EV distribtion */
float delta_approx(float x) {
  float abs = x > 0 ? 1.f-x/EPSILON : 1.f+x/EPSILON;
  return ( abs > 0 ? abs/EPSILON : 0 );
}

// Erase the content of a complex square matrix with side length l
void nullify(const int l, MKL_Complex8 * m) {
  for(int i=0;i<l*l;++i) {
    m[i].real = 0.;
    m[i].imag = 0.;
  }
}

// Given a (NxN) matrix m compute {m,.}
void Arrange_Anticommutator(MKL_Complex8 * m, MKL_Complex8 * acomm)
{
  nullify(SWEEP, acomm);

  /* Copy Hx1I */
  for(int i=0;i<N;++i) {
    for(int j=0;j<N;++j) {
      for(int k=0;k<N;++k) {
        acomm[(i*SWEEP+j)*N+k*(SWEEP+1)].real += m[i*N+j].real;
        acomm[(i*SWEEP+j)*N+k*(SWEEP+1)].imag += m[i*N+j].imag;
      }
    }
  }
  /* Copy +1IxH^T */
  for(int k=0;k<N;++k) {
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        acomm[(j*SWEEP+i)+k*N*(SWEEP+1)].real += m[i*N+j].real;
        acomm[(j*SWEEP+i)+k*N*(SWEEP+1)].imag += m[i*N+j].imag;
      }
    }
  }
}


// Given a (NxN) matrix m compute [m,.]
void Arrange_Commutator(MKL_Complex8 * m, MKL_Complex8 * comm)
{
  nullify(SWEEP, comm);

  /* Copy Lx1I */
  for(int i=0;i<N;++i) {
    for(int j=0;j<N;++j) {
      for(int k=0;k<N;++k) {
        comm[(i*SWEEP+j)*N+k*(SWEEP+1)].real += m[i*N+j].real;
        comm[(i*SWEEP+j)*N+k*(SWEEP+1)].imag += m[i*N+j].imag;
      }
    }
  }
  /* Copy -1IxL^T */
  for(int k=0;k<N;++k) {
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        comm[(j*SWEEP+i)+k*N*(SWEEP+1)].real -= m[i*N+j].real;
        comm[(j*SWEEP+i)+k*N*(SWEEP+1)].imag -= m[i*N+j].imag;
      }
    }
  }
}


// Arrange Matrix D (kn^2 x kn^2) from NUM_H Matrices H and NUM_L Matrices L of type (P,Q)
void Arrange_Dirac_Matrix(float complex *gamma_passed, MKL_Complex8 *Matrices, MKL_Complex8 *Matrix_Operators, MKL_Complex8 *Dirac, int NUM_H, int NUM_L)
{

  /* Cast gammas to use old code */
  /* HAS TO BE CHANGED! */
  int size_gamma = pow(2,D-1)*K*K;
  int *gamma = (int*) calloc(size_gamma,sizeof(int));
  for(int i=0;i<size_gamma;++i) {
    if(abs(creal(gamma_passed[i])-1)<TOLERANCE) gamma[i]=1;
    else if(abs(creal(gamma_passed[i])+1)<TOLERANCE) gamma[i]=-1;
    else if(abs(cimag(gamma_passed[i])-1)<TOLERANCE) gamma[i]=2;
    else if(abs(cimag(gamma_passed[i])+1)<TOLERANCE) gamma[i]=-2;
  }


  int offset; /* Offset to scan through array of Matrix Operators */
  int off_gamma; /* Offset to start with the right Gamma Matrix in case of Matrices L */

  /* Set the Dirac Matrix to zero and calculate (Anti-)Commutators from the Matrices (L) M */
  nullify(K*SWEEP,Dirac);
  for(int i=0;i<NUM_H;++i) {
    Arrange_Anticommutator(&Matrices[i*SWEEP],&Matrix_Operators[i*SWEEP*SWEEP]);
  }
  for(int i=0;i<NUM_L;++i) {
    Arrange_Commutator(&Matrices[i*SWEEP+NUM_H*SWEEP],&Matrix_Operators[i*SWEEP*SWEEP+NUM_H*SWEEP*SWEEP]);
  }

  /* COPY {H,.}'s */
  for(int n=0;n<NUM_H;++n) { /* Matrices H */
    offset = n*SWEEP*SWEEP;

    for(int II=0;II<K;++II) { /* Elements of Gamma Matrix, Blocks of Dirac Matrix -- Rows */
      for(int JJ=0;JJ<K;++JJ) { /* Elements of Gamma Matrix, Blocks of Dirac Matrix -- Columns */
        if(gamma[n*K*K+II*K+JJ]==0) {
            continue;
        } else if (gamma[n*K*K+II*K+JJ]==1) {
            for(int i=0;i<SWEEP;++i) { /* Rows */
              for(int j=0;j<SWEEP;++j) { /* Columns */
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real += Matrix_Operators[i*SWEEP+j+offset].real;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag += Matrix_Operators[i*SWEEP+j+offset].imag;
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==-1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real -= Matrix_Operators[i*SWEEP+j+offset].real;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag -= Matrix_Operators[i*SWEEP+j+offset].imag;
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real -= Matrix_Operators[i*SWEEP+j+offset].imag;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag += Matrix_Operators[i*SWEEP+j+offset].real;
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==-2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real += Matrix_Operators[i*SWEEP+j+offset].imag;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag -= Matrix_Operators[i*SWEEP+j+offset].real;
              }
            }
        } /* End II loop */
      } /* End JJ loop */
    } /* End loop over Gamma Matrix Elements */
  } /* End loop over Matrices */

  
  /* COPY [L,.]'s */
  for(int n=0;n<NUM_L;++n) { /* Matrices L */
    offset = n*SWEEP*SWEEP+NUM_H*SWEEP*SWEEP;
    off_gamma = NUM_H*K*K;

    for(int II=0;II<K;++II) { /* Elements of Gamma Matrix, Blocks of Dirac Matrix -- Rows */
      for(int JJ=0;JJ<K;++JJ) { /* Elements of Gamma Matrix, Blocks of Dirac Matrix -- Columns */
        if(gamma[off_gamma+n*K*K+II*K+JJ]==0) {
            continue;
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==2) {
            for(int i=0;i<SWEEP;++i) { /* Rows */
              for(int j=0;j<SWEEP;++j) { /* Columns */
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real += Matrix_Operators[i*SWEEP+j+offset].real;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag += Matrix_Operators[i*SWEEP+j+offset].imag;
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==-2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real -= Matrix_Operators[i*SWEEP+j+offset].real;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag -= Matrix_Operators[i*SWEEP+j+offset].imag;
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==-1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real -= Matrix_Operators[i*SWEEP+j+offset].imag;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag += Matrix_Operators[i*SWEEP+j+offset].real;
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].real += Matrix_Operators[i*SWEEP+j+offset].imag;
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)].imag -= Matrix_Operators[i*SWEEP+j+offset].real;
              }
            }
        } /* End II loop */
      } /* End JJ loop */
    } /* End loop over Gamma Matrix Elements */
  } /* End loop over Matrices */

  free(gamma);

}


// Initialise all Matrices, their Eigenvalues and the action
void Matrices_Initialisation(struct pcg32_random_t *rng, MKL_Complex8 *Matrices, float *action, int NUM_H, int NUM_L)
{
  /* Set zero-temperature initial state for all *
   * matrices and caluclate initial action      */
   
  int itimesN;
  int Nplus1 = N+1;
  int offset;

  for (int n=0;n<NUM_M;++n) {
    offset = n*SWEEP;
    for (int i=0;i<N;++i) {
      Matrices[i*Nplus1+offset].real = MAX_ELEMENT * (pcg32_boundedrand_r(&rng[n],2)?-1:1) * ldexp(pcg32_random_r(&rng[n]),-32);
      Matrices[i*Nplus1+offset].imag = 0.0f;
      itimesN = i*N;
      for (int j=i+1;j<N;j++) {
        Matrices[itimesN+j+offset].real = MAX_ELEMENT * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
        Matrices[itimesN+j+offset].imag = MAX_ELEMENT * (pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
        Matrices[j*N+i+offset].real =  Matrices[itimesN+j+offset].real;
        Matrices[j*N+i+offset].imag = -Matrices[itimesN+j+offset].imag; /* This is hermitian! */
      }
    }
  }

  //*action = traceD2(Matrices, NUM_H, NUM_L);
  //*action = G2 * traceD2(Matrices, NUM_H, NUM_L) + G4 * traceD4(Matrices, NUM_H, NUM_L);
}


// Creates a new Markov chain element
void Get_Next_MCMC_Element(struct pcg32_random_t *rng, MKL_Complex8 *Matrices, float *action,
                           int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L, int *acc) 
{
  int pos_x, pos_y;
  int pos_upper, pos_lower;
  float p;                    // Random float for accepting MC elemement
  MKL_Complex8 temp;          // To save random value that is changed
  float delta_action;

  int itimesN;
  int Nplus1 = N+1;
  int offset;

  /* For each Matrix change a value in the upper triangle randomly *
   * calculate the the change of the action and finally decide if  *
   * the new matrix should be accpted.                             */
  for(int n=0;n<NUM_M;++n)
  {
    /* Set the offset to write to the right matrix */
    offset = n*SWEEP;

    /* Calculate random float in [0,1) for Monte Carlo Move Decision */
    p = ldexp(pcg32_random_r(&rng[n]),-32);

    /* Calculate two random integers and generate position in upper and lower half */
    pos_x = pcg32_boundedrand_r(&rng[n],N);
    pos_y = pcg32_boundedrand_r(&rng[n],N);
    pos_upper = pos_x<=pos_y ? pos_x*N+pos_y : pos_y*N+pos_x;
    pos_lower = pos_x>pos_y ? pos_x*N+pos_y : pos_y*N+pos_x;

    temp.real = STEP_SIZE*(pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
    temp.imag = STEP_SIZE*(pcg32_boundedrand_r(&rng[n],2)?-1:1)*ldexp(pcg32_random_r(&rng[n]),-32);
    temp.real+= Matrices[pos_upper+offset].real;
    temp.imag+= Matrices[pos_upper+offset].imag;
    if(temp.real>MAX_ELEMENT || temp.real<-MAX_ELEMENT) continue;
    if(temp.imag>MAX_ELEMENT || temp.imag<-MAX_ELEMENT) continue;

    //delta_action = delta_action_traceD2(Matrices, n, temp, pos_x, pos_y, NUM_H, NUM_L);
    delta_action  = G2 * delta_action_traceD2(Matrices, n, temp, pos_x, pos_y, NUM_H, NUM_L);
    delta_action += G4 * delta_action_traceD4(Matrices, n, temp, pos_x, pos_y, sigmaAB, sigmaABCD, NUM_H, NUM_L);

    /* Finally test if new action is smaller or except randomly if exp supressed, *
     * if yes write new element in upper and lower half and copy new eigenvalues  */
    if(delta_action<=0 || expf(-delta_action)>p) {
      *acc += 1;
      *action += delta_action;
      if(pos_x != pos_y) {
        Matrices[pos_upper+offset].real = temp.real;
        Matrices[pos_upper+offset].imag = temp.imag;
        Matrices[pos_lower+offset].real = temp.real;
        Matrices[pos_lower+offset].imag = -temp.imag;
      } else {
        Matrices[pos_upper+offset].real = temp.real;
      }
    }
  }

}


void Measure_Eigenvalues_Dirac(float complex *Gamma_Matrices, MKL_Complex8 *Matrices,
		MKL_Complex8 *Matrix_Operators, MKL_Complex8 *Dirac,
		float *evs_D, float *evs_D_avrg, float *evs_D_avrg2, int NUM_H, int NUM_L)
{
  Arrange_Dirac_Matrix(Gamma_Matrices,Matrices,Matrix_Operators,Dirac,NUM_H,NUM_L);
  LAPACKE_cheev(LAPACK_ROW_MAJOR,'N','U',K*SWEEP,Dirac,K*SWEEP,evs_D);

  for (int i=0;i<K*SWEEP;i++) {
    evs_D_avrg[i]  += evs_D[i];
    evs_D_avrg2[i] += evs_D[i]*evs_D[i];
  }
}


void Measure_Eigenvaluedistribution_Dirac(float *support_points, float *evs_D, float *dist_evs_D_avrg)
{
  for (int i = 0; i < POINTS; i++)
    for (int j = 0; j < K*SWEEP; j++)
      dist_evs_D_avrg[i] += delta_approx(support_points[i]-evs_D[j]);
}

void Measure_Orderparameter_Frac(MKL_Complex8 *Matrices, double *frac,
                                 double *frac_squared, int NUM_H)
{
  double numerator = 0.0;
  double denominator = 0.0;

  double *traceH = (double*) calloc(NUM_H,sizeof(double));
  double *traceH2 = (double*) calloc(NUM_H,sizeof(double));

  /* Calculate Tr H_i for all i */
  for(int num=0;num<NUM_H;++num) {
    int offset = num*SWEEP;
    for(int i=0;i<N;++i) traceH[num] += (double) Matrices[i*(N+1)+offset].real;
  }


  /* Calculate Tr H^2_i for all i */
  for(int num=0;num<NUM_H;++num) {
    int offset = num*SWEEP;
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        traceH2[num] += (double) Matrices[i*N+j+offset].real*Matrices[j*N+i+offset].real;
        traceH2[num] -= (double) Matrices[i*N+j+offset].imag*Matrices[j*N+i+offset].imag;
      }
    }
  }

  /* Calculate the numerator and denominator */
  for(int num=0;num<NUM_H;++num) {
    numerator += traceH[num]*traceH[num];
    denominator += traceH2[num];
  }

  *frac         += numerator/denominator;
  *frac_squared += numerator*numerator/(denominator*denominator);

  free(traceH);
  free(traceH2);
}
