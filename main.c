#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <string.h>
// #include <mpi.h>
// #include <omp.h>
#include "Precision.h"
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

// Fundamental Parameters
#define N 5
#define P 3
#define Q 1
#define G2 (1.0)
#define G4 (1.0)

// User Parameters
#define CHAIN_LENGTH (1000*N*N)
#define WRITEOUT_FREQ (1*N*N) // in SWEEPs

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

  // Initialise parameters
  // TODO: This needs NUM_H and NUM_L to be set, so Generate_Clifford_Odd_Group
  // need to have been called by now (otherwise we loose const-ness).
  // Can this be improved?
  struct Matrix_Properties parameters = {
    .n = N,
    .p = P,
    .q = Q,
    .d = P + Q,
    .s = ( Q - P + 64 ) % 8,
    .k = (P+Q) % 2 ? (int) pow( 2, ( P + Q - 1 ) / 2 ) : (int) pow( 2, ( P + Q ) / 2 ),
    .num_h = 0,
    .num_l = 0,
    .g2 = G2,
    .g4 = G4
  };

  // Generate the ODD Clifford Group and number of (anti-)hermitian matrices
  // Hermitian Matrices are stored first, anti-hermitian matrices second
  size_t const size_gammas = (int) pow( 2, parameters.d - 1 ) * parameters.k * parameters.k;
  float complex *Gamma_Matrices = (float complex *) malloc( size_gammas * sizeof(float complex) );
  Generate_Clifford_Odd_Group(Gamma_Matrices, &parameters);
  size_t const num_m = parameters.num_h + parameters.num_l;

  // Allocate matrix SigmaAB and calulate its values from the Gammas
  int *SigmaAB = (int*) calloc( num_m * num_m, sizeof(int) );
  Calculate_Trace_Gamma_ABAB(Gamma_Matrices, parameters, SigmaAB);

  // Allocate matrix sigmaABCD and calulate its values from the Gammas
  int **SigmaABCD = (int **) malloc( num_m * sizeof(int *) );
  // Adjust that size for less memory use and higher dimension d:
  for( size_t i = 0; i < num_m; ++i )
  {
    SigmaABCD[i] = (int *) malloc(
        7 * parameters.d * parameters.d * parameters.d * parameters.d * sizeof(int)
        );
  }
  Calculate_Trace_Gamma_ABCD(Gamma_Matrices, parameters, SigmaABCD);

  // Melissa O'NEILL's RNG lib seeded with sys time and mem address
  // http://www.pcg-random.org/
  struct pcg32_random_t rngs[num_m];
  for(int i=0;i<num_m;++i) {
    pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
  }

  // Array allocations
  REAL complex *Matrices;
  Matrices = (REAL complex *) calloc( num_m * parameters.n * parameters.n, sizeof(REAL complex) );

  // Memory needed for H's and L's, and gamma matrices
  const int memory = num_m * parameters.n * parameters.n * sizeof(REAL complex)
                   + size_gammas * sizeof(float complex);

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

  double action = Matrices_Initialisation( rngs, Matrices, parameters, SigmaAB, SigmaABCD );

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
    fprintf(stdout, "  Chain elements to be produced    : %zu sweep\n",
        CHAIN_LENGTH / parameters.n / parameters.n);
    fprintf(stdout, "  Starting time                    : %s\n\n", buff);

    fprintf(stdout, "  Geometry:\n");
    fprintf(stdout, "  ---------\n\n");

    fprintf(stdout, "  Type                             : (%zu,%zu)\n", parameters.p, parameters.q);
    fprintf(stdout, "  Matrix Size N                    : %d\n", N);
    fprintf(stdout, "  Action                           : S = Tr( %.1f * D^2 + %.1f * D^4 )\n", G2, G4);
    fprintf(stdout, "  Initial value                    : S = %f\n\n", action);

    fprintf(stdout, "  Dimension                        : %zu, (K = %zu)\n", parameters.d, parameters.k);
    fprintf(stdout, "  Signature                        : %zu\n", parameters.s);
    fprintf(stdout, "  Number      Hermitian Matrices H : %zu\n", parameters.num_h);
    fprintf(stdout, "  Number Anti-Hermitian Matrices L : %zu\n\n", parameters.num_l);

    printf("===============================================\n\n");

    fprintf(stdout, "  GENERATING CHAIN:\n\n");
    fprintf(stdout,
        "  time \t rank \t sweep \t action S \t acceptance (accumulated) \t acceptance (last %zu sweep)\n",
        WRITEOUT_FREQ / parameters.n / parameters.n
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
      double acc_rate_accum = (double) accepted / ( t * num_m );
      double acc_rate_sweep = (double) accepted_sweep / ( WRITEOUT_FREQ * num_m );
      accepted_old = accepted;
      tune_step_size( acc_rate_sweep, &step_size );

      // time, rank, sweep, action S, acceptance (accumulated), acceptance (last sweeps)
      action = Calculate_Action( Matrices, parameters, SigmaAB, SigmaABCD );
      fprintf(
          stdout,
          "  %s \t %3d \t %5lu \t %.6f \t %4.2f \t %4.2f\n",
          buff, rank, t / parameters.n / parameters.n, action, 100.0 * acc_rate_accum, 100.0 * acc_rate_sweep
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
  for( int i = 0; i < num_m; ++i )
  {
    free(SigmaABCD[i]);
  }
  free(Gamma_Matrices);
  free(Matrices);

  // MPI_Finalize();
  return 0;
}
