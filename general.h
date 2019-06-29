#ifndef GENERAL_H
#define GENERAL_H
#include <math.h>
float aleatorio(float seed);
float gaussiana(float mu, float sigma, float seed);
double norma2(double *x);
int delta_x(double *x, int i, int j, double *delta_r);
#endif
