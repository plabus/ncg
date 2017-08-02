#ifndef _ACTIONS_H_
#define _ACTIONS_H_

float tr1(float complex *Matrices, int pos1, int NUM_H);
float tr2(float complex *Matrices, int pos1, int pos2);
float tr3(float complex *Matrices, int pos1, int pos2, int pos3);

float tr2_real_ij(float complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
float tr2_imag_ij(float complex *Matrices, int pos1, int pos2, int pos_x, int pos_y);
float tr3_real_ij(float complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
float tr3_imag_ij(float complex *Matrices, int pos1, int pos2, int pos3, int pos_x, int pos_y);
float row_norm_squared(float complex *Matrices, int pos1, int pos_row);

float traceD2(float complex *Matrices, int NUM_H, int NUM_L);
float traceD4(float complex *Matrices, int NUM_H, int NUM_L);

float delta_action_traceD2(float complex *Matrices, int position, float complex temp, int pos_x, int pos_y, int NUM_H, int NUM_L);
float delta_action_traceD4(float complex *Matrices, int positionA, float complex temp, int pos_x, int pos_y, int *sigmaAB, int **sigmaABCD, int NUM_H, int NUM_L);
#endif
