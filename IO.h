#ifndef _MY_IO_H_
#define _MY_IO_H_

#include <complex.h>
#include <argp.h>
#include "Precision.h"
#include "Matrix_Properties.h"

struct arguments
{
  size_t N;
  size_t P;
  size_t Q;
  double G2;
  double G4;
  size_t chain_length;
  size_t number_chains;
  size_t writeout_freq;
};

struct arguments parse_command_line_args(
    int argc,
    char *argv[]
    );

void write_matrices_to_file(
    REAL complex const *M,               // array of matrices
    struct Matrix_Properties const prop, // includes num_h, num_l, n and k
    FILE* file_ptr                       // pointer to file for write to disk
    );

#endif
