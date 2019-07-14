#ifndef AVANZAR_H
#define AVANZAR_H
#include <math.h>
int apply_PBC(double *x, int N, double L);
int position_verlet(double *x, double *v, int N, double h, double *F);
int velocity_verlet(double *v, int N, double h, double *F, double *F2);
int fuerza_PCB(double *delta_r, double L);
double fuerzas(double *tabla_F, double *tabla_V, double *F, double *F2, double *x, double rc2, double dr2, int N, double L);
double step_verlet(double *x, double *v, double *F, double *F2, double *tabla_F, double *tabla_V, double rc2, double dr2, double h, double L, int N);
double presion(double *tabla_F, double *tabla_V, double *x, double rc2, double dr2, int N, double L, double T, double rho);
#endif
