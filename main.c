#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <string.h>
// #include <mpi.h>
// #include <omp.h>
#include "Constants.h"
#include "Utilities.h"
#include "Clifford.h"
#include "Random.h"
#include "MonteCarlo.h"
#include "Actions.h"
#include "Matrix_Properties.h"

// Double escape in order to print
// precision REAL as string
#define xstr(a) str(a)
#define str(a) #a

int main()
{

/////////////////////////////////////////////////////////////////////
//                                                                 //
//             DECLARATION / ALLOCATION / SEEDING RNG              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  char buff[100];
  time_t now = time(0);
  strftime( buff, 100, "%Y-%m-%d %H:%M:%S", localtime(&now) );

  int start_sweep = 0;
  int rank = 0;
  // int nproc = 0;
  int NUM_H, NUM_L;


  // Generate the ODD Clifford Group and number of (anti-)hermitian matrices
  // Hermitian Matrices are stored first, anti-hermitian matrices second
  int size_gammas = (int) pow(2, D-1) * K * K;
  float complex *Gamma_Matrices = (float complex*) malloc(size_gammas*sizeof(float complex));
  Generate_Clifford_Odd_Group(P, Q, Gamma_Matrices, &NUM_H, &NUM_L);

  // Allocate matrix SigmaAB and calulate its values from the Gammas
  int *SigmaAB = (int*) calloc(NUM_M*NUM_M,sizeof(int));
  Calculate_Trace_Gamma_ABAB(Gamma_Matrices, SigmaAB, NUM_H);

  // Allocate matrix sigmaABCD and calulate its values from the Gammas
  int **SigmaABCD = (int **)malloc(NUM_M*sizeof(int*));
  // Adjust that size for less memory use and higher D:
  for(int i=0; i<NUM_M; i++) SigmaABCD[i] = (int *)malloc(7*D*D*D*D * sizeof(int));
  Calculate_Trace_Gamma_ABCD(Gamma_Matrices, SigmaABCD, NUM_H);


  // Melissa O'NEILL's RNG lib seeded with sys time and mem address
  // http://www.pcg-random.org/
  struct pcg32_random_t rngs[NUM_M];
  for(int i=0;i<NUM_M;++i) {
    pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
  }


  // Array allocations
  REAL complex *Matrices;
  Matrices = (REAL complex *) calloc(NUM_M*SWEEP,sizeof(REAL complex));

  // Memory needed for H's and L's, and gamma matrices
  const int memory = NUM_M * SWEEP * sizeof(REAL complex) + size_gammas * sizeof(float complex);

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                  MPI INIT AND LOCAL VARIABLES                   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  // MPI_Init(NULL, NULL);
  // MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

/////////////////////////////////////////////////////////////////////
//                                                                 //
//     INIT OF RANDOM MATRICES AND EXP VALUES OF OBSERVABLES       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  const double start_time = cclock();

  const struct Matrix_Properties parameters = {
    .num_h = NUM_H,
    .num_l = NUM_L,
    .n = N,
    .k = K
  };
  double action = Matrices_Initialisation( rngs, Matrices, parameters );

  if(rank==0)
  {
    printf("===============================================\n\n");

    fprintf(stdout, "  Simulation details:\n");
    fprintf(stdout, "  -------------------\n\n");

    fprintf(stdout, "  Compiled in precision            : %s\n", xstr(REAL));
    // fprintf(stdout, "  MPI Processes                    : %d\n", nproc);
    // fprintf(stdout, "  OMP Threads                      : %d\n", omp_get_max_threads());
    fprintf(stdout, "  Memory used                      : %.3g MiB\n", memory/1048576.);
    fprintf(stdout, "  Starting sweep                   : %d\n", start_sweep);
    fprintf(stdout, "  Chain elements to be produced    : %d sweep\n", CHAIN_LENGTH/SWEEP);
    fprintf(stdout, "  Starting time                    : %s\n\n", buff);

    fprintf(stdout, "  Geometry:\n");
    fprintf(stdout, "  ---------\n\n");

    fprintf(stdout, "  Type                             : (%d,%d)\n", P, Q);
    fprintf(stdout, "  Matrix Size N                    : %d\n", N);
    fprintf(stdout, "  Action                           : S = Tr( %.1f * D^2 + %.1f * D^4 )\n", G2, G4);
    fprintf(stdout, "  Initial value                    : S = %f\n\n", action);

    fprintf(stdout, "  Dimension                        : %d, (K = %d)\n", D, K);
    fprintf(stdout, "  Signature                        : %d\n", S);
    fprintf(stdout, "  Number      Hermitian Matrices H : %d\n", NUM_H);
    fprintf(stdout, "  Number Anti-Hermitian Matrices L : %d\n\n", NUM_L);

    printf("===============================================\n\n");

    fprintf(stdout, "  GENERATING CHAIN:\n\n");
    fprintf(stdout,
        "  time \t rank \t sweep \t action S \t acceptance (accumulated) \t acceptance (last %d sweep)\n",
        WRITEOUT_FREQ/SWEEP
        );
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                      CHAIN CREATION                             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  uint64_t accepted     = 0;   // counter for total accepted steps
  uint64_t accepted_old = 0;   // buffer to calculate accepted steps per sweep
  double step_size      = 0.1; // initial step length

  for( uint64_t t = 0; t < CHAIN_LENGTH; ++t )
  {
    // Generate new chain element
    Get_Next_MCMC_Element( rngs, Matrices, parameters, SigmaAB, SigmaABCD, &accepted, step_size );

    // Print some diagnostics each WRITEOUT_FREQ SWEEPs
    if( t % WRITEOUT_FREQ == 0 && t != 0 )
    {
      // Get the time of the day (hrs, min, sec) as string
      time_t now = time(0);
      strftime( buff, 100, "%H:%M:%S", localtime(&now) );

      // Calculate number of accepts in last sweep, as well as
      // the accumulated and recent (last WRITEOUT_FREQ sweep) acceptance rates,
      // finally tune the step size according to recent acceptance
      // (Remember, we are counting accepts for each matrix and each step t)
      uint64_t accepted_sweep = accepted - accepted_old;
      double acc_rate_accum = (double) accepted / ( t * NUM_M );
      double acc_rate_sweep = (double) accepted_sweep / ( WRITEOUT_FREQ * NUM_M );
      accepted_old = accepted;
      tune_step_size( acc_rate_sweep, &step_size );

      // time, rank, sweep, action S, acceptance (accumulated), acceptance (last sweeps)
      action = Calculate_Action( Matrices, parameters );
      fprintf(
          stdout,
          "  %s \t %3d \t %5lu \t %.6f \t %4.2f \t %4.2f\n",
          buff, rank, t/SWEEP, action, 100.0*acc_rate_accum, 100.0*acc_rate_sweep
          );
    }
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                       FINALISATION                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  if(rank==0)
  {
    fprintf(stdout, "\n");
    fprintf(stdout, "  Total time of simulation %.3g min\n\n", ((cclock()-start_time)/60.));
    fprintf(stdout, "===============================================\n");
  }

  free(SigmaAB);
  for( int i = 0; i < NUM_M; ++i )
  {
    free(SigmaABCD[i]);
  }
  free(Gamma_Matrices);
  free(Matrices);

  // MPI_Finalize();
  return 0;
}
