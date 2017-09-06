#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "Precision.h"
#define EPSILON 1E-2

double cclock();
void printM(const int L, REAL complex * M);
int compare_floats(const void * a, const void * b);
float delta_approx(float);
void nullify(const int l, float complex * m);

#endif
