#ifndef LAMBDA_H
#define LAMBDA_H
#include <math.h>
double lambda(double *x, double dl, int N);
double H_boltzmann (double *v, int N);
double distribucion_radial(double *x, double *g_r, double dr_g, double L, int N);
#endif
