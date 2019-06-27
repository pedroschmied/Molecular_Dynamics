#ifndef VISUALIZACION_H
#define VISUALIZACION_H

#include "math.h"

int save_lammpstrj(char *filename, double *x, double *v, int N, float L, int frame);
int load_frame(void *fp, double *x, double *v, int N, float *L);
int load_lammpstrj(char *filename, double *x, double *v, int N, float *L, int frame);

#endif
