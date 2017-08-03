#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <mpi.h>
#include <omp.h>
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
  float complex *Gamma_Matrices = (float complex*) malloc(size_gammas*sizeof(float complex));
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
  float complex *Matrices;
  Matrices = (float complex *) calloc(NUM_M*SWEEP,sizeof(float complex));

  /* Memory needed for H's and L's, and gamma matrices */
  int memory = NUM_M * SWEEP * sizeof(float complex) + size_gammas * sizeof(float complex);

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                  MPI INIT AND LOCAL VARIABLES                   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  if(rank==0)
  {
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
  if( rank==0 ) fprintf(stdout, "  BURN IN\n");

  /* Create burn in chain */
  for( int t = 0; t < BURN_IN; t++ ) {
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

  MPI_Finalize();
  return 0;
}
