#ifndef INTERACCION_H
#define INTERACCION_H
#include <math.h>
double tablas (double *V,double *Fuerza, double rc2, int largo_tabla);
double pair_force (double *tabla_F, double *tabla_V, double rij2, double dr2, double *F_mod);
#endif
