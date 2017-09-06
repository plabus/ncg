#include <complex.h>
#include "Measurements.h"

// Given a (NxN) matrix m compute {m,.}
void Arrange_Anticommutator(float complex * m, float complex * acomm)
{
  nullify(SWEEP, acomm);

  /* Copy Hx1I */
  for(int i=0;i<N;++i) {
    for(int j=0;j<N;++j) {
      for(int k=0;k<N;++k) {
        acomm[(i*SWEEP+j)*N+k*(SWEEP+1)] += m[i*N+j];
      }
    }
  }
  /* Copy +1IxH^T */
  for(int k=0;k<N;++k) {
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        acomm[(j*SWEEP+i)+k*N*(SWEEP+1)] += m[i*N+j];
      }
    }
  }
}


// Given a (NxN) matrix m compute [m,.]
void Arrange_Commutator(float complex * m, float complex * comm)
{
  nullify(SWEEP, comm);

  /* Copy Lx1I */
  for(int i=0;i<N;++i) {
    for(int j=0;j<N;++j) {
      for(int k=0;k<N;++k) {
        comm[(i*SWEEP+j)*N+k*(SWEEP+1)] += m[i*N+j];
      }
    }
  }
  /* Copy -1IxL^T */
  for(int k=0;k<N;++k) {
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        comm[(j*SWEEP+i)+k*N*(SWEEP+1)] -= m[i*N+j];
      }
    }
  }
}


// Arrange Matrix D (kn^2 x kn^2) from NUM_H Matrices H and NUM_L Matrices L of type (P,Q)
void Arrange_Dirac_Matrix(float complex *gamma_passed, float complex *Matrices, float complex *Matrix_Operators, float complex *Dirac, int NUM_H, int NUM_L)
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
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] += Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==-1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] -= Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] += I * Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[n*K*K+II*K+JJ]==-2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] -= I * Matrix_Operators[i*SWEEP+j+offset];
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
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] += Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==-2) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] -= Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==-1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] += I * Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } else if (gamma[off_gamma+n*K*K+II*K+JJ]==1) {
            for(int i=0;i<SWEEP;++i) {
              for(int j=0;j<SWEEP;++j) {
                Dirac[i*K*SWEEP+j+(II*K*SWEEP*SWEEP+JJ*SWEEP)] -= I * Matrix_Operators[i*SWEEP+j+offset];
              }
            }
        } /* End II loop */
      } /* End JJ loop */
    } /* End loop over Gamma Matrix Elements */
  } /* End loop over Matrices */

  free(gamma);

}




void Measure_Eigenvalues_Dirac(float complex *Gamma_Matrices, float complex *Matrices,
		float complex *Matrix_Operators, float complex *Dirac,
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

void Measure_Orderparameter_Frac(float complex *Matrices, double *frac,
                                 double *frac_squared, int NUM_H)
{
  double numerator   = 0.0;
  double denominator = 0.0;

  double *traceH  = (double*) calloc( NUM_H, sizeof(double) );
  double *traceH2 = (double*) calloc( NUM_H, sizeof(double) );

  /* Calculate Tr H_i for all i */
  for(int num=0;num<NUM_H;++num) {
    int offset = num*SWEEP;
    for(int i=0;i<N;++i) traceH[num] += (double) creal( Matrices[i*(N+1)+offset] );
  }


  /* Calculate Tr H^2_i for all i */
  for(int num=0;num<NUM_H;++num) {
    int offset = num*SWEEP;
    for(int i=0;i<N;++i) {
      for(int j=0;j<N;++j) {
        traceH2[num] += (double) creal( Matrices[i*N+j+offset] ) * creal( Matrices[j*N+i+offset] );
        traceH2[num] -= (double) cimag( Matrices[i*N+j+offset] ) * cimag( Matrices[j*N+i+offset] );
      }
    }
  }

  /* Calculate the numerator and denominator */
  for(int num=0;num<NUM_H;++num) {
    numerator   += traceH[num] * traceH[num];
    denominator += traceH2[num];
  }

  *frac         += numerator / denominator;
  *frac_squared += numerator * numerator / ( denominator * denominator );

  free(traceH);
  free(traceH2);
}
