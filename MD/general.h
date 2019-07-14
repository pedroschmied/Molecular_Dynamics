#ifndef GENERAL_H
#define GENERAL_H
#include <math.h>
double aleatorio();
float gaussiana(double mu, double sigma);
double norma2(double *x);
int delta_x(double *x, int i, int j, double *delta_r);
double mean (double *v, int o, int n, int k);
double mean2 (double *v, int o, int n, int k);
double std2 (double *v, int o, int n, int k);
/*
long double mean (double *v, int o, int n, int k);
long double mean2 (double *v, int o, int n, int k);
long double std2 (double *v, int o, int n, int k);*/
#endif
