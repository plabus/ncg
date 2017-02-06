#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <mpi.h>
#include <omp.h>
#include "mkl_lapacke.h"
#include "Constants.h"
#include "Clifford.h"
#include "Random.h"
#include "MonteCarlo.h"

int main() {

/////////////////////////////////////////////////////////////////////
//                                                                 //
//             DECLARATION / ALLOCATION / SEEDING RNG              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  char buff[100];
  time_t now = time(0);
  strftime(buff,100,"%Y-%m-%d %H:%M:%S",localtime(&now));

  float action;
  float incr;
  int accepted;
  double start;
  int rank, nproc, rest;
  int loc_msrments, loc_t_max;
  int stepsize, stepsize_burn;
  int print_steps = 10; /* Number of steps for percentage printout */
  int NUM_H, NUM_L;


  /* Generate the ODD Clifford Group and number of (anti-)hermitian matrices *
   * Hermitian Matrices are stored first, antihermitian matrices second      */
  int size_gammas = (int) pow(2, D-1) * K * K;
  float complex *Gamma_Matrices = (float complex*)
                                  malloc(size_gammas*sizeof(float complex));
  Generate_Clifford_Odd_Group(P, Q, Gamma_Matrices, &NUM_H, &NUM_L);

  /* Allocate matrix SigmaAB and calulate its values from the Gammas */
  int *SigmaAB = (int*) calloc(NUM_M*NUM_M,sizeof(int));
  Calculate_Trace_Gamma_ABAB(Gamma_Matrices, SigmaAB, NUM_H);

  /* Allocate matrix sigmaABCD and calulate its values from the Gammas */
  int **SigmaABCD = (int **)malloc(NUM_M*sizeof(int*));
  // Adjust that size for less memory use and higher D:
  for(int i=0; i<NUM_M; i++) SigmaABCD[i] = (int *)malloc(7*D*D*D*D * sizeof(int));
  Calculate_Trace_Gamma_ABCD(Gamma_Matrices, SigmaABCD, NUM_H);


  /* Melissa O'NEILL's RNG lib seeded with sys time and mem address *
   * http://www.pcg-random.org/                                     */
  struct pcg32_random_t rngs[NUM_M];
  for(int i=0;i<NUM_M;++i) {
    pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
  }


  /* Array allocations */
  MKL_Complex8 *Matrices;
  Matrices = (MKL_Complex8 *) calloc(NUM_M*SWEEP,sizeof(MKL_Complex8));

#if MEASURE_EVS
  MKL_Complex8 *Matrix_Operators;
  MKL_Complex8 *Dirac_Matrix;

  float *Eigenvalues_Dirac;
  float *Eigenvalues_Dirac_Average;
  float *Eigenvalues_Dirac_Average_Squared;

  Matrix_Operators = (MKL_Complex8 *) malloc(NUM_M*SWEEP*SWEEP*sizeof(MKL_Complex8));
  Dirac_Matrix  = (MKL_Complex8 *) malloc(K*K*SWEEP*SWEEP*sizeof(MKL_Complex8));

  Eigenvalues_Dirac = (float *) malloc(K*SWEEP*sizeof(float));
  Eigenvalues_Dirac_Average = (float *) calloc(K*SWEEP,sizeof(float));
  Eigenvalues_Dirac_Average_Squared = (float *) calloc(K*SWEEP,sizeof(float));

  FILE *fp_evs;
  fp_evs  = fopen("evs", "w+");
#endif

#if MEASURE_DIST
  float *Support_Points;
  float *Distribution_Eigenvalues_Average;

  Support_Points = (float *) malloc(POINTS*sizeof(float));
  Distribution_Eigenvalues_Average = (float *) calloc(POINTS,sizeof(float));

  FILE *fp_dist;
  fp_dist = fopen("dist", "w+");
#endif

#if MEASURE_FRAC
  double frac = 0.0;
  double frac_squared = 0.0;
  FILE *fp_frac;
  fp_frac = fopen("frac", "a+");
#endif

  /* Note that this is only the memory allocated here. That does not include 
   * memory allocate in any of the routines, in particular CHEEV */
  int memory = (NUM_M*SWEEP + NUM_M*SWEEP*SWEEP + K*K*SWEEP*SWEEP) *
                 sizeof(MKL_Complex8) +
               (3*K*SWEEP + 2*POINTS) * sizeof(float) +
	       size_gammas * sizeof(float complex);
  #if MEASURE_FRAC
  #if MEASURE_EVS==0
  #if MEASURE_DIST==0
  memory = NUM_M*SWEEP*sizeof(MKL_Complex8);
  #endif
  #endif
  #endif

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                  MPI INIT AND LOCAL VARIABLES                   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Support points for distribution of EVs */
  incr = 2.f * MAX_DIST / POINTS;
#if MEASURE_DIST
  for (int i = 0; i < POINTS; i++) {
    Support_Points[i] = - MAX_DIST + i*incr;
  }
#endif

  loc_msrments = MEASUREMENTS/nproc;
  rest = MEASUREMENTS%nproc;
  if(rank<rest) loc_msrments++;

  loc_t_max = (2*TAU*(loc_msrments-1)+1);
  stepsize_burn = BURN_IN/print_steps;
  stepsize = loc_t_max/print_steps;

/////////////////////////////////////////////////////////////////////
//                                                                 //
//     INIT OF RANDOM MATRICES AND EXP VALUES OF OBSERVABLES       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  start = cclock();
  Matrices_Initialisation(rngs, Matrices, &action, NUM_H, NUM_L);

  if(rank==0) {
    printf("===============================================\n\n");
    fprintf(stdout, "  Type: (%d,%d)\n", P, Q);
    fprintf(stdout, "  Dimension: %d, (K = %d)\n", D, K);
    fprintf(stdout, "  Signature: %d\n", S);
    fprintf(stdout, "  Action: S = Tr(%.1f*D^2 + %.1f*D^4)\n", G2, G4);
    fprintf(stdout, "\n");
    fprintf(stdout, "  Number Hermitian Matrices H: %d\n", NUM_H);
    fprintf(stdout, "  Number Anti-Hermitian Matrices L: %d\n", NUM_L);
    fprintf(stdout, "  Matrix Size N: %d\n", N);
    fprintf(stdout, "\n");
    fprintf(stdout, "  Measurements: %d\n", MEASUREMENTS);
    fprintf(stdout, "  Points for Distribution: %d\n", POINTS);
    fprintf(stdout, "  Max x-value for Distribution: %.1f\n", MAX_DIST);
    fprintf(stdout, "\n");
    fprintf(stdout, "  MPI Processes: %d\n", nproc);
    fprintf(stdout, "  OMP Threads: %d\n", omp_get_max_threads());
    fprintf(stdout, "  Memory used: %.3g MiB\n", memory/1048576.);
    fprintf(stdout, "\n");
    fprintf(stdout, "  Starting Calculation...\n");
    fprintf(stdout, "  Starting time: %s\n", buff);
    fprintf(stdout, "\n");
    printf("===============================================\n\n");
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                          BURN IN                                //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  accepted = 0;
  if( rank==0) fprintf(stdout, "  BURN IN\n");

  /* Create burn in chain */
  for (int t = 0; t < BURN_IN; t++) {
    Get_Next_MCMC_Element(rngs, Matrices, &action, SigmaAB, SigmaABCD,
                          NUM_H, NUM_L, &accepted);

    /* Print percentage of progress at root process */
    if( t%stepsize_burn==0 && rank==0 )
      fprintf(stdout, "  After %7.3f min: ~ %5.0f%% done\n",
              (cclock()-start)/60., (double)t/BURN_IN*100.);
  }

  if(rank==0) fprintf(stdout, "  After %7.3f min: ~ %5.0f%% done\n",
                      (cclock()-start)/60., 100.);

  /* Some user information */
  if (rank==0) {
    fprintf(stdout, "  Acceptance rate during burn in: %.3f\n",
                              (float)accepted/(BURN_IN*NUM_M));
    fprintf(stdout, "\n");
    fprintf(stdout, "===============================================\n\n");
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                        MEASUREMENTS                             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  /* Measure every 2*TAU time steps after burn in:        */
  /* Eigenvalues Dirac Operator 		          */
  /* Eigenvalue Distribution			          */
  /* Orderparameter F = Sum(Tr H_i)^2 / Tr(Sum(H_i^2) / N */

  accepted = 0;
  if( rank==0) fprintf(stdout, "  MEASUREMENTS\n");

  for(int t=0;t<loc_t_max;t++) {

    /* Create the actual Markov chain */
    Get_Next_MCMC_Element(rngs, Matrices, &action, SigmaAB, SigmaABCD,
                          NUM_H, NUM_L, &accepted);
			   
    /* Measurements */
    if( t%(2*TAU)==0 ) {

      #if MEASURE_EVS
      Measure_Eigenvalues_Dirac(Gamma_Matrices, Matrices, Matrix_Operators, Dirac_Matrix,
		                Eigenvalues_Dirac, Eigenvalues_Dirac_Average,
				Eigenvalues_Dirac_Average_Squared, NUM_H, NUM_L);
      #endif

      #if MEASURE_DIST
      Measure_Eigenvaluedistribution_Dirac(Support_Points, Eigenvalues_Dirac,
                                           Distribution_Eigenvalues_Average);
      #endif

      #if MEASURE_FRAC
      Measure_Orderparameter_Frac(Matrices, &frac, &frac_squared, NUM_H);
      #endif

    } // Measurement

    /* Print percentage of progress at root process */
    if( t%stepsize==0 && rank==0 )
      fprintf(stdout, "  After %7.3f min: ~ %5.0f%% done\n",
              (cclock()-start)/60., (double)t/loc_t_max*100.);

  }

  /* Some user information */
  if (rank==0) {
    fprintf(stdout, "  Acceptance rate during measurements: %.3f\n",
            (float)accepted/(loc_t_max*NUM_M));
    fprintf(stdout, "\n");
    fprintf(stdout, "===============================================\n");
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                        POSTPROCESSING                           //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#if MEASURE_EVS
  for (int i=0;i<K*SWEEP;i++) {
    Eigenvalues_Dirac_Average[i] /= (float) MEASUREMENTS;
    Eigenvalues_Dirac_Average_Squared[i] /= (float) MEASUREMENTS;
  }

  MPI_Allreduce(MPI_IN_PLACE,Eigenvalues_Dirac_Average,K*SWEEP,MPI_FLOAT,
                MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,Eigenvalues_Dirac_Average_Squared,K*SWEEP,
                MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
#endif

#if MEASURE_DIST
  for (int i=0;i<POINTS;++i) {
    Distribution_Eigenvalues_Average[i] /= (float) MEASUREMENTS*K*SWEEP;
  }

  MPI_Allreduce(MPI_IN_PLACE,Distribution_Eigenvalues_Average,POINTS,
                MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
#endif

#if MEASURE_FRAC
  frac = frac/(double)(MEASUREMENTS*N);
  frac_squared = frac_squared/(double)(MEASUREMENTS*N*N);
  MPI_Allreduce(MPI_IN_PLACE,&frac,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&frac_squared,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                   OUTPUT / WRITE TO FILE                        //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  if(rank==0)
  {

#if MEASURE_EVS
    /* Printing EVs of Dirac Operator with yerrors */
    fprintf(stdout, "\n  Average Eigenvalues of Dirac Operator:\n");
    for(int i=0;i<K*SWEEP;i++) {
      fprintf(fp_evs, "%d\t%.8f\t%.8f\n", i+1, Eigenvalues_Dirac_Average[i],
              sqrtf( (Eigenvalues_Dirac_Average_Squared[i] - 
              Eigenvalues_Dirac_Average[i]*Eigenvalues_Dirac_Average[i])/MEASUREMENTS) );
    }
    fprintf(stdout, "  Written to file...\n");
#endif

#if MEASURE_DIST
    /* Printing the distribution of EVs of Dirac Operator */
    fprintf(stdout, "\n  Eigenvalue Distribution of Dirac Operator:\n");
    for(int i=0;i<POINTS;i++) {
      fprintf(fp_dist, "%.5f\t%.9f\n", Support_Points[i], 
              Distribution_Eigenvalues_Average[i]);
    }
    fprintf(stdout, "  Written to file...\n");
#endif

#if MEASURE_FRAC
    /* Appending the orderparameter value F to file */
    fprintf(stdout, "\n  Order parameter F:\n");
    fprintf(fp_frac, "%.3f\t%.9f\t%.9f\n", G2, frac,
            sqrt((frac_squared - frac*frac)/MEASUREMENTS) );
    fprintf(stdout, "  Written to file...\n");
#endif

  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                       FINALISATION                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  if(rank==0) {
    fprintf(stdout, "\n");
    fprintf(stdout, "  Time to solution: %.3g min\n\n", ((cclock()-start)/60.));
    fprintf(stdout, "===============================================\n");
  }

  free(SigmaAB);
  for(int i=0; i<NUM_M; i++) free(SigmaABCD[i]);
  free(Gamma_Matrices);
  free(Matrices);

#if MEASURE_EVS
  fclose(fp_evs);
  free(Matrix_Operators);
  free(Dirac_Matrix);
  free(Eigenvalues_Dirac);
  free(Eigenvalues_Dirac_Average);
  free(Eigenvalues_Dirac_Average_Squared);
#endif

#if MEASURE_DIST
  fclose(fp_dist);
  free(Distribution_Eigenvalues_Average);
  free(Support_Points);
#endif

#if MEASURE_FRAC
  fclose(fp_frac);
#endif

  MPI_Finalize();
  return 0;
}
