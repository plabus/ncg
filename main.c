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
#include "Matrix_Properties.h"
#include "Utilities.h"
#include "Random.h"
#include "IO.h"
#include "Clifford.h"
#include "Actions.h"
#include "MonteCarlo.h"

// Double escape in order to print
// precision REAL as string
#define xstr(a) str(a)
#define str(a) #a


int main(int argc, char *argv[])
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

  struct arguments args = parse_command_line_args(argc, argv);

  // Initialise simulation parameters
  struct Matrix_Properties parameters = {
    .n = args.N,
    .p = args.P,
    .q = args.Q,
    .d = args.P + args.Q,
    .s = ( args.Q - args.P + 64 ) % 8,
    .k = ( args.P + args.Q ) % 2 ?
      (int) pow( 2, ( args.P + args.Q - 1 ) / 2 ) : //  odd case
      (int) pow( 2, ( args.P + args.Q ) / 2 ),      // even case
    .num_h = 0,
    .num_l = 0,
    .num_m = 0,
    .g2 = args.G2,
    .g4 = args.G4
  };

  // TODO
  // generate file name from input parameters
  // need append or create mode, depending whether we "restart" or begin a chain
  // shall we keep the file open during the whole simulation or open it only when writing?
  FILE *file_ptr;
  file_ptr = fopen("configs", "w+");
  fprintf(file_ptr, "%zu %zu %zu %5.3f %5.3f\n",
      parameters.n,
      parameters.p,
      parameters.q,
      parameters.g2,
      parameters.g4
      );

  // Generate the ODD Clifford Group and number of (anti-)hermitian matrices,
  // where hermitian matrices are stored first, anti-hermitian matrices second.
  // Update the struct parameters during the process (num_h, num_l and num_m).
  size_t const size_gammas = (int) pow( 2, parameters.d - 1 ) * parameters.k * parameters.k;
  float complex *Gamma_Matrices = (float complex *) malloc( size_gammas * sizeof(float complex) );
  Generate_Clifford_Odd_Group(Gamma_Matrices, &parameters);

  // Allocate matrix SigmaAB and calculate its values from the Gammas
  int *SigmaAB = (int*) calloc( parameters.num_m * parameters.num_m, sizeof(int) );
  Calculate_Trace_Gamma_ABAB(Gamma_Matrices, parameters, SigmaAB);

  // Allocate matrix sigmaABCD and calculate its values from the Gammas
  int **SigmaABCD = (int **) malloc( parameters.num_m * sizeof(int *) );
  // Adjust that size for less memory use and higher dimension d:
  for( size_t i = 0; i < parameters.num_m; ++i )
  {
    SigmaABCD[i] = (int *) malloc(
        7 * parameters.d * parameters.d * parameters.d * parameters.d * sizeof(int)
        );
  }
  Calculate_Trace_Gamma_ABCD(Gamma_Matrices, parameters, SigmaABCD);

  // Melissa O'NEILL's RNG lib seeded with sys time and mem address
  // http://www.pcg-random.org/
  struct pcg32_random_t rngs[parameters.num_m];
  for( size_t i = 0; i < parameters.num_m; ++i )
  {
    pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
  }

  // Array allocations
  REAL complex *Matrices;
  Matrices = (REAL complex *) calloc( parameters.num_m * parameters.n * parameters.n, sizeof(REAL complex) );

  // Memory needed for H's and L's, and gamma matrices
  const int memory = parameters.num_m * parameters.n * parameters.n * sizeof(REAL complex)
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
        args.chain_length / parameters.n / parameters.n);
    fprintf(stdout, "  Starting time                    : %s\n\n", buff);

    fprintf(stdout, "  Geometry:\n");
    fprintf(stdout, "  ---------\n\n");

    fprintf(stdout, "  Type                             : (%zu,%zu)\n", parameters.p, parameters.q);
    fprintf(stdout, "  Matrix Size N                    : %zu\n", parameters.n);
    fprintf(stdout, "  Action                           : S = Tr( %.1f * D^2 + %.1f * D^4 )\n", parameters.g2, parameters.g4);
    fprintf(stdout, "  Initial value                    : S = %f\n\n", action);

    fprintf(stdout, "  Dimension                        : %zu, (K = %zu)\n", parameters.d, parameters.k);
    fprintf(stdout, "  Signature                        : %zu\n", parameters.s);
    fprintf(stdout, "  Number      Hermitian Matrices H : %zu\n", parameters.num_h);
    fprintf(stdout, "  Number Anti-Hermitian Matrices L : %zu\n\n", parameters.num_l);

    printf("===============================================\n\n");

    fprintf(stdout, "  GENERATING CHAIN:\n\n");
    fprintf(stdout,
        "  time       chain   sweep   action S \t acceptance (accumulated) \t acceptance (last %zu sweep)\n",
        args.writeout_freq / parameters.n / parameters.n
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

  for( uint64_t t = 0; t < args.chain_length; ++t )
  {
    // Generate new chain element
    Get_Next_MCMC_Element( rngs, Matrices, parameters, SigmaAB, SigmaABCD, &accepted, step_size );

    // Print some diagnostics each writeout_freq SWEEPs
    if( t % args.writeout_freq == 0 && t != 0 )
    {
      // Get the time of the day (hrs, min, sec) as string
      time_t now = time(0);
      strftime( buff, 100, "%H:%M:%S", localtime(&now) );

      // Calculate number of accepts in last sweep, as well as
      // the accumulated and recent (last writeout_freq sweep) acceptance rates,
      // finally tune the step size according to recent acceptance
      // (Remember, we are counting accepts for each matrix and each step t)
      uint64_t accepted_sweep = accepted - accepted_old;
      double acc_rate_accum = (double) accepted / ( t * parameters.num_m );
      double acc_rate_sweep = (double) accepted_sweep / ( args.writeout_freq * parameters.num_m );
      accepted_old = accepted;
      tune_step_size( acc_rate_sweep, &step_size );

      action = Calculate_Action( Matrices, parameters, SigmaAB, SigmaABCD );

      // Screen diagnostics from every process
      // time, rank, sweep, action S, acceptance (accumulated), acceptance (last sweeps)
      fprintf(
          stdout,
          "  %s %3d \t %5lu \t     %.6f \t %4.2f \t %4.2f\n",
          buff, rank+1, t / parameters.n / parameters.n, action, 100.0 * acc_rate_accum, 100.0 * acc_rate_sweep
          );

      // Write sweep and action for your chain to a file
      fprintf( file_ptr, "%zu %.12f\n",
          t / parameters.n / parameters.n,
          action
          );

      // Write matrices for your chain to a file
      write_matrices_to_file( Matrices, parameters, file_ptr );
    }
  }

/////////////////////////////////////////////////////////////////////
//                                                                 //
//                       FINALISATION                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  fclose(file_ptr);

  if(rank==0)
  {
    fprintf(stdout, "\n");
    fprintf(stdout, "  Total time of simulation %.3g min\n\n", ((cclock()-start_time)/60.));
    fprintf(stdout, "===============================================\n");
  }

  free(SigmaAB);
  for( size_t i = 0; i < parameters.num_m; ++i )
  {
    free(SigmaABCD[i]);
  }
  free(Gamma_Matrices);
  free(Matrices);

  // MPI_Finalize();
  return 0;
}
