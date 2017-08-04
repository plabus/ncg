#ifndef _CONSTANTS_
#define _CONSTANTS_

// Floating point precision use
// (either double or float)
#define REAL double

// Fundamental Parameters
#define N 25
#define P 1
#define Q 0
#define G2 (1.0)
#define G4 (0.0)

// Derived Quantities
#define D (P+Q)
#define S ((Q-P+8)%8)
#define K ( (D%2)?(int)pow(2,(D-1)/2):(int)pow(2,D/2) )
#define SWEEP (N*N)
#define NUM_M ((int)pow(2,(D-1)))

// User Parameters
#define CHAIN_LENGTH (20*SWEEP)
#define WRITEOUT_FREQ 2 // in SWEEPs

#endif
