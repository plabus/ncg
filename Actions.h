#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include "Constants.h"

double tr1(REAL complex *Matrices, int pos1, int NUM_H);
double tr2(REAL complex *Matrices, int pos1, int pos2);
double tr3(REAL complex *Matrices, int pos1, int pos2, int pos3);

double tr2_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
double tr2_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
double tr3_real_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
double tr3_imag_ij(REAL complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
double row_norm_squared(REAL complex *Matrices, int pos1, int pos_row);

double traceD2(REAL complex *Matrices, int NUM_H, int NUM_L);
double traceD4(REAL complex *Matrices, int NUM_H, int NUM_L);

double delta_action_traceD2(REAL complex *Matrices, int position, REAL complex temp, int pos_x, int pos_y, int NUM_H, int NUM_L);
double delta_action_traceD4(REAL complex *Matrices, int positionA, REAL complex temp, int pos_x, int pos_y, int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L);
#endif
