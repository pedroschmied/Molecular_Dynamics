#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"
#include "termalizacion.h"
#include "magnitudes.h"

#define PI  3.14159

int main()
{
	int N = 512;
	double rho = 0.8442, L = cbrt(N / rho);

	double *x, *v;
	x = (double*) malloc(3 * N * sizeof(double));
	v = (double*) malloc(3 * N * sizeof(double));

//Tablas de Fuerza y Potencial
	int largo_tabla = 500000;
	double  *tabla_V, *tabla_F;
	tabla_V = (double*) malloc(largo_tabla * sizeof(double));
	tabla_F = (double*) malloc(largo_tabla * sizeof(double));
	double rc2 = 2.5 * 2.5 ;
	double dr2 = tablas (tabla_V, tabla_F, rc2, largo_tabla);
//----------------------
	double  *F;
	F = (double*) malloc(3 * N * sizeof(double));
	double  *F2;
	F2 = (double*) malloc(3 * N * sizeof(double));
//______set inicial______
	float T_gauss =  2.0;
	set_pos(x, N, L);
	double cinetica0 = set_vel(v, N, T_gauss);
	double temp0;
	fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
	temp0 = cinetica0 * 2.0 / (3.0 * (double)N);
//______datos variables__
	int pasos = 100;
	int loops_T = 20000;
	int termalizacion_inicial = 2000;
	int termalizacion = 100;
	int correlacion = 250;
	int T_medidas = 200;
	double Tf = 0.1;
	double h = 0.001;
//_______________________
	double T0 = temp0, T = T0, dT, pot, factor_T, P;
	float va;
	int i, t, c, n, l = 0;

	double  *potencial, *cinetica, *energia, *cinetica2;
	potencial = (double*) malloc((pasos)* sizeof(double));
	cinetica = (double*) malloc((pasos) * sizeof(double));
	cinetica2 = (double*) malloc((termalizacion / 2) * sizeof(double));
	energia = (double*) malloc((pasos)* sizeof(double));

	double *mean_V, *mean_K, *mean_E, *std2_V, *std2_K, *std2_E, *Cv, *Presion;
	mean_V = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	mean_K = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	mean_E = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	std2_V = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	std2_K = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	std2_E = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	Cv = (double*) malloc((loops_T / T_medidas)* sizeof(double));
	Presion = (double*) malloc((loops_T / T_medidas)* sizeof(double));

	dT = (Tf - T0) / (double)(loops_T);
//__termalización inicial
	for (n = 0; n < termalizacion_inicial; n++)
	{
		step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
	}
	printf("\nTermalizado (inicial)\n");
//_______________________
	for (t = 0; t < loops_T; t++)
	{
		va = (float) t * 100.0 / (float)loops_T;
		printf("Progreso %.2f", va);
		printf("%%\n");
//__termalizo___
		T0 = 0.0;
		for(n = 0; n < termalizacion / 2; n++)
		{
 			step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
		}
		for(n = 0; n < termalizacion / 2; n++)
		{
 			step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
			*(cinetica2 + n) = 0.0;
			for(i = 0; i < 3 * N; i++)
			{
				*(cinetica2 + n) += *(v + i) * *(v + i) / (2.0 * (double)N);
			}
		}
//______________
		if (t % T_medidas == 0)
		{
			P = 0.0;
			for (n = 0; n < pasos; n++)
			{
				*(cinetica + n) = 0.0;
	//__descorrelaciono_____
				for(c = 0; c < correlacion; c++)
				{
					pot = step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
				}
	//______________________
				for(i = 0; i < 3 * N; i++)
				{
					*(cinetica + n) += *(v + i) * *(v + i) / 2.0;
				}
				*(potencial + n) = pot / (double)N;
				*(cinetica + n) = *(cinetica + n) / (double)N;		
				*(energia + n) = *(potencial + n) + *(cinetica + n);
				P += presion(tabla_F, tabla_V, x, rc2, dr2, N, L, T, rho);
			}
			*(Presion + l) = P / pasos;
			*(mean_V + l) = mean(potencial, 0, pasos, 1);
			*(std2_V + l) = std2(potencial, 0, pasos, 1);
			*(mean_K + l) = mean(cinetica, 0, pasos, 1);
			*(std2_K + l) = std2(cinetica, 0, pasos, 1);
			*(mean_E + l) = mean(energia, 0, pasos, 1);
			*(std2_E + l) = std2(energia, 0, pasos, 1);
			*(Cv + l) = C_v(cinetica, pasos, N);

			T0 = *(mean_K + l) * 2.0 / 3.0;
			l += 1;

		}
		else
		{
			T0 = mean(cinetica2, 0, termalizacion / 2, 1);
			T0 = T0 * 2.0 / 3.0;
		}
		dT = (Tf - T0) / (double)(loops_T - t);
		T += dT;
		factor_T = (double)sqrt(T / T0);
		for (i = 0; i < 3 * N; i++)
		{
			*(v + i) = *(v + i) * factor_T;
		}
	}


///------------------------------------------

	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/new/prueba.txt");
	fp = fopen(filename, "w");
	l = 0;
	for (t = 0; t < loops_T; t++)
	{
		if (t % T_medidas == 0)
		{
	//		fprintf(fp, "%.15lf\t", *(mean_K + t) * 2.0 / 3.0);
			fprintf(fp, "%d\t", t);
			fprintf(fp, "%.15lf\t", *(mean_V + l));
			fprintf(fp, "%.15lf\t", *(mean_K + l));
			fprintf(fp, "%.15lf\t", *(mean_E + l));
			fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_V + l)));
			fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_K + l)));
			fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_E + l)));
			fprintf(fp, "%.15lf\t", *(Cv + l));
			fprintf(fp, "%.15lf\n", *(Presion + l));
			l += 1;
		}
	}
	fclose(fp);
	free(Presion);
	free(mean_V);
	free(mean_K);
	free(mean_E);
	free(std2_V);
	free(std2_K);
	free(std2_E);
	free(Cv);
	free(cinetica2);
	free(cinetica);
	free(potencial);
	free(energia);
	free(F);
	free(F2);
	free(tabla_F);
	free(tabla_V);	
	free(x);
	free(v);
	return 0;
}
#include "general.c"
#include "inicializar.c"
#include "visualizacion.c"
#include "interaccion.c"
#include "avanzar.c"
#include "termalizacion.c"
#include "magnitudes.c"
