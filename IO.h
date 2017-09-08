#ifndef _MY_IO_H_
#define _MY_IO_H_

#include <argp.h>

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

#endif
