#include <stdlib.h>
#include "IO.h"

const char *argp_program_version = "NCG 0.9";
const char *argp_program_bug_address = "<plabus@sissa.it>";
static char args_doc[] = "";
static char doc[] =
"\nMarkov Chain Monte Carlo configuration generator for non-commutative geometries.\n\n\
All arguments are optional and default values are written in [=brackets].";

static struct argp_option options[] =
{
  { "N", 1, "10", OPTION_ARG_OPTIONAL, "Side length of the matrices." },
  { "P", 2, "1", OPTION_ARG_OPTIONAL, "First value in Lorentz signiture (p,q)." },
  { "Q", 3, "0", OPTION_ARG_OPTIONAL, "Second value in Lorentz signiture (p,q)." },
  { "G2", 4, "1.0", OPTION_ARG_OPTIONAL, "Action parameter: S = g2 Tr D^2 + g4 Tr D^4." },
  { "G4", 5, "1.0", OPTION_ARG_OPTIONAL, "Action parameter: S = g2 Tr D^2 + g4 Tr D^4." },
  { "chain-length", 6, "100", OPTION_ARG_OPTIONAL, "Length (in sweeps = N * N) of the one single chain to be generated." },
  { "number-chains", 7, "1", OPTION_ARG_OPTIONAL,
    "Number of MCMC chains to be generate (equals number of MPI processes)." },
  { "write-out-freq", 8, "1", OPTION_ARG_OPTIONAL,
    "Number of sweep (= N * N) between every configuartion write out." },
  { 0 }
};


static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;
  switch (key)
  {
    case 1: arguments->N = atoi(arg); break;
    case 2: arguments->P = atoi(arg); break;
    case 3: arguments->Q = atoi(arg); break;
    case 4: arguments->G2 = atof(arg); break;
    case 5: arguments->G4 = atof(arg); break;
    case 6: arguments->chain_length = atoi(arg); break;
    case 7: arguments->number_chains = atoi(arg); break;
    case 8: arguments->writeout_freq = atoi(arg); break;
    case ARGP_KEY_ARG: return 0;
    default: return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

struct arguments parse_command_line_args(
    int argc,
    char *argv[]
    )
{
  struct arguments arguments;

  // Set defaut values for command line options
  arguments.N = 10;
  arguments.P = 1;
  arguments.Q = 0;
  arguments.G2 = 1.0;
  arguments.G4 = 1.0;
  arguments.number_chains = 1;
  arguments.chain_length = 100;
  arguments.writeout_freq = 1;

  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  arguments.chain_length  *= arguments.N * arguments.N;
  arguments.writeout_freq *= arguments.N * arguments.N;

  return arguments;
}

void write_matrices_to_file(
    REAL complex const *M,               // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    FILE* file_ptr                       // pointer to file for write to disk
    )
{
  for( size_t index = 0; index < prop.n * prop.n; ++index )
  {
    fprintf( file_ptr, "%.16f %.16f ",  creal(M[index]), cimag(M[index]) );
  }
  fprintf( file_ptr, "\n" );
}
