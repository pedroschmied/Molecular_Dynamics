#include "interaccion.h"
double tablas (double *tabla_V,double *tabla_F, double rc2, int largo_tabla)
{
	double rf2 = 2.5 * 2.5; //sigma = 1
	double dr2 = rf2 / ((double)largo_tabla + 1.0), r02 = 0.0 + dr2; //así tengo 500 000 ptos y el r02 = dr2 para simplificar cuentas.
	int l, lf = (int) ((rf2 - r02) / dr2) + 1;
	double rij2 = 0.0, r6, rc6 = (1.0 / rc2) * (1.0 / rc2) * (1.0 / rc2);
	double Vrc = 4.0 * (rc6 * rc6 - rc6);

	for (l = 0; l < lf; l++)
	{
		rij2 += dr2;
		r6 = (1.0 / rij2) * (1.0 / rij2) * (1.0 / rij2);
		*(tabla_V + l) = 4.0 * (r6 * r6 - r6) - Vrc;
		*(tabla_F + l) = (24.0 / rij2) * (2.0 * r6 * r6 - r6);
	}
	return dr2;
}

double pair_force (double *tabla_F, double *tabla_V, double rij2, double dr2, double *F_mod)
{
	int k = (int)((rij2 - dr2) / dr2);
	if(k < 0)
	{
		k = 0;
	}
	*F_mod = (*(tabla_F + k + 1) - *(tabla_F + k)) * (rij2 - k * dr2) / dr2 + *(tabla_F + k);
	double Vij = (*(tabla_V + k + 1) - *(tabla_V + k)) * (rij2 - k * dr2) / dr2 + *(tabla_V + k);
	return Vij;
}
