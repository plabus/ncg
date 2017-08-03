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
  strftime( buff, 100, "%Y-%m-%d %H:%M:%S", localtime(&now) );

  float action;
  double start_time;
  int start_sweep = 0;
  int rank, nproc;
  int NUM_H, NUM_L;


  /* Generate the ODD Clifford Group and number of (anti-)hermitian matrices *
   * Hermitian Matrices are stored first, anti-hermitian matrices second      */
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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//     INIT OF RANDOM MATRICES AND EXP VALUES OF OBSERVABLES       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  start_time = cclock();
  Matrices_Initialisation(rngs, Matrices, &action, NUM_H, NUM_L);

  if(rank==0)
  {
    printf("===============================================\n\n");

    fprintf(stdout, "  Simulation details:\n");
    fprintf(stdout, "  -------------------\n\n");

    fprintf(stdout, "  MPI Processes                    : %d\n", nproc);
    fprintf(stdout, "  OMP Threads                      : %d\n", omp_get_max_threads());
    fprintf(stdout, "  Memory used                      : %.3g MiB\n", memory/1048576.);
    fprintf(stdout, "  Starting sweep                   : %d\n", start_sweep);
    fprintf(stdout, "  Chain elements to be produced    : %d sweep\n", CHAIN_LENGTH/SWEEP);
    fprintf(stdout, "  Starting time                    : %s\n\n", buff);

    fprintf(stdout, "  Geometry:\n");
    fprintf(stdout, "  ---------\n\n");

    fprintf(stdout, "  Type                             : (%d,%d)\n", P, Q);
    fprintf(stdout, "  Matrix Size N                    : %d\n", N);
    fprintf(stdout, "  Action                           : S = Tr( %.1f * D^2 + %.1f * D^4 )\n\n", G2, G4);

    fprintf(stdout, "  Dimension                        : %d, (K = %d)\n", D, K);
    fprintf(stdout, "  Signature                        : %d\n", S);
    fprintf(stdout, "  Number      Hermitian Matrices H : %d\n", NUM_H);
    fprintf(stdout, "  Number Anti-Hermitian Matrices L : %d\n\n", NUM_L);

    printf("===============================================\n\n");

    fprintf(stdout, "  GENERATING CHAIN:\n\n");
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                      CHAIN CREATION                             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  int accepted     = 0;   // counter for total accepted steps
  int accepted_old = 0;   // buffer to calculate accepted steps per sweep
  float step_size  = 0.5; // initial step length

  for( int t = 0; t < CHAIN_LENGTH; ++t )
  {
    // Generate new chain element
    Get_Next_MCMC_Element(rngs, Matrices, &action, SigmaAB, SigmaABCD, NUM_H, NUM_L, &accepted, step_size);

    // Print some diagnostics at each SWEEP
    if( t % SWEEP == 0 && t != 0 )
    {
      // Get the time of date as string
      time_t now = time(0);
      strftime( buff, 100, "%Y-%m-%d %H:%M:%S", localtime(&now) );

      // Calculate number of accepts in last sweep, as well as
      // the accumulated and recent (last sweep) acceptance rates,
      // finally tune the step size according to recent acceptance
      // (Remember, we are counting accepts for each matrix and each step t)
      int accepted_sweep    = accepted - accepted_old;
      double acc_rate_accum = (double) accepted / ( t * NUM_M );
      double acc_rate_sweep = (double) accepted_sweep / ( SWEEP * NUM_M );
      accepted_old = accepted;
      tune_step_size( acc_rate_sweep, &step_size );

      fprintf(
          stdout,
          "  %s, rank %3d, sweep %5d, S = %.3f, \
          acceptance = %4.2f%% (accumulated), acceptance = %4.2f%% (last sweep)\n",
          buff, rank, t/SWEEP, action, 100.0*acc_rate_accum, 100.0*acc_rate_sweep
          );
    }
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                       FINALISATION                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  if(rank==0) {
    fprintf(stdout, "\n");
    fprintf(stdout, "  Time to solution: %.3g min\n\n", ((cclock()-start_time)/60.));
    fprintf(stdout, "===============================================\n");
  }

  free(SigmaAB);
  for(int i=0; i<NUM_M; i++) free(SigmaABCD[i]);
  free(Gamma_Matrices);
  free(Matrices);

  MPI_Finalize();
  return 0;
}
