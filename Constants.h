#ifndef _CONSTANTS_
#define _CONSTANTS_

// Floating point precision use
// (either double or float)
#define REAL double

/* Fundamental Parameters */
#define N 5
#define P 1
#define Q 3
#define G2 0.f
#define G4 1.f

/* Measured Quantities */
#define MEASURE_EVS 1
#define MEASURE_DIST 1 /* Use only if MEASURE_EVS */
#define MEASURE_FRAC 1

/* Derived Quantities */
#define D (P+Q)
#define S ((Q-P+8)%8)
#define K ( (D%2)?(int)pow(2,(D-1)/2):(int)pow(2,D/2) )
#define SWEEP (N*N)
#define CHAIN_LENGTH (20*SWEEP)
#define TAU (2*SWEEP)
#define NUM_M (NUM_H+NUM_L)
#define TOLERANCE (1E-2)

/* Measurement Parameters */
#define MEASUREMENTS 2000
// This defines the maximum value for the
// bucketing used in the distribution measurement.
// This heuristic formula may be changed.
#define MAX_DIST (3.f+(D>4?D-4:0)*1.f)
#define POINTS (2000*2*(int)MAX_DIST)
#define EPSILON 1E-2
// The following determines if "Large N"
// is applied, i.e. if terms of O(N^0)
// are ignored when calculating \Delta S
// or not
// Large N   : NOT commented out (erase the
// leading "//"
// All terms : commented out
// #define LARGE_N

#endif
