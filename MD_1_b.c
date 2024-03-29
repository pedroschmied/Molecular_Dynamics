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
	float rho = 0.8442, L = cbrt(N / rho);

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

	int i, t, c, n;
	int pasos = 100, termalizacion = 2000, correlacion = 500;
	double h = 0.001, pot;
	double  *potencial, *cinetica, *energia;
	potencial = (double*) malloc((pasos)* sizeof(double));
	cinetica = (double*) malloc((pasos) * sizeof(double));
	energia = (double*) malloc((pasos)* sizeof(double));
	double *mean_V, *mean_K, *mean_E, *std2_V, *std2_K, *std2_E, *Cv;
	float va;
	double T0 = 1.0, dT = 0.05, Tf = 2.0;
	float T_gauss =  T0;
	int temperaturas = (int)((Tf - T0) / dT) + 1;

	mean_V = (double*) malloc((temperaturas)* sizeof(double));
	mean_K = (double*) malloc((temperaturas)* sizeof(double));
	mean_E = (double*) malloc((temperaturas)* sizeof(double));
	std2_V = (double*) malloc((temperaturas)* sizeof(double));
	std2_K = (double*) malloc((temperaturas)* sizeof(double));
	std2_E = (double*) malloc((temperaturas)* sizeof(double));
	Cv = (double*) malloc((temperaturas)* sizeof(double));

	for (t = 0; t < temperaturas; t++)
	{
	//-------set up inicial
		set_pos(x, N, L);
		set_vel(v, N, T_gauss);
		fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
	//---------------------
		for (n = 0; n < termalizacion; n++)
		{
			va = (float)(n) * 100.0 / (float)(pasos * correlacion + termalizacion);
			printf("Progreso %.2f", va);
			printf("%%\r");
			step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
		}
		for (n = 0; n < pasos; n++)
		{
			*(cinetica + n) = 0.0;
			*(potencial + n) = 0.0;
			for(c = 0; c < correlacion; c++)
			{
				va = (float)(n * correlacion + c + termalizacion) * 100.0 / (float)(pasos * correlacion + termalizacion);
				printf("Progreso %.2f", va);
				printf("%%\r");
				pot = step_verlet(x, v, F, F2, tabla_F, tabla_V, rc2, dr2, h, L, N);
			}
			*(potencial + n) = pot / (double)N;
			for(i = 0; i < 3 * N; i++)
			{
				*(cinetica + n) += *(v + i) * *(v + i) / 2.0;
			}
			*(cinetica + n) = *(cinetica + n) / (double)N;
			*(energia + n) = *(cinetica + n) + *(potencial + n);
		}
		*(mean_V + t) = mean(potencial, 0, pasos, 1);
		*(std2_V + t) = std2(potencial, 0, pasos, 1);
		*(mean_K + t) = mean(cinetica, 0, pasos, 1);
		*(std2_K + t) = std2(cinetica, 0, pasos, 1);
		*(mean_E + t) = mean(energia, 0, pasos, 1);
		*(std2_E + t) = std2(energia, 0, pasos, 1);
		*(Cv + t) = C_v(cinetica, pasos, N);
		printf("\n%d/%d\n", t, temperaturas - 1);
//		printf("T_gauss = %f\t", T_gauss);
//		printf("T real = %lf\n", *(mean_K + t) * 2.0 / 3.0);
		T_gauss += dT;
	}
///------------------------------------------

	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/MD_datos_1_b.txt");
	fp = fopen(filename, "w");
	for (t = T0; t < temperaturas; t++)
	{
		fprintf(fp, "%.15lf\t", *(mean_K + t) * 2.0 / 3.0);
		fprintf(fp, "%.15lf\t", *(mean_V + t));
		fprintf(fp, "%.15lf\t", *(mean_K + t));
		fprintf(fp, "%.15lf\t", *(mean_E + t));
		fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_V + t)));
		fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_K + t)));
		fprintf(fp, "%.15lf\t", (double)sqrt(*(std2_E + t)));
		fprintf(fp, "%.15lf\n", *(Cv + t));
	}
	fclose(fp);
	free(mean_V);
	free(mean_K);
	free(mean_E);
	free(std2_V);
	free(std2_K);
	free(std2_E);
	free(Cv);
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
