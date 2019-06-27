#ifndef AVANZAR_H
#define AVANZAR_H
#include <math.h>
int apply_PBC(double *x, int N, float L);
int position_verlet(double *x, double *v, int N, double h, double *F);
int velocity_verlet(double *v, int N, double h, double *F, double *F2);
int fuerza_PCB(double *delta_r, float L);
double fuerzas(double *tabla_F, double *tabla_V, double *F, double *F2, double *x, double rc2, double dr2, int N, float L);
#endif
