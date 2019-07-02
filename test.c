#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"

int main()
{
	int N = 512;
	float rho = 0.8442, L = cbrt(N/rho);
	float T_gauss =  2.0;
	double *x, *v;
	x = (double*) malloc(3 * N * sizeof(double));
	v = (double*) malloc(3 * N * sizeof(double));

//Tablas de Fuerza y Potencial
	int largo_tabla = 500000;
	double  *tabla_V, *tabla_F;
	tabla_V = (double*) malloc(largo_tabla * sizeof(double));
	tabla_F = (double*) malloc(largo_tabla * sizeof(double));
	double rc2 = 2.5 * 2.5 ;//L / 4.0;//?????????????????????;
	double dr2 = tablas (tabla_V, tabla_F, rc2, largo_tabla);
//----------------------
	double  *F;
	F = (double*) malloc(3 * N * sizeof(double));
	double  *F2;
	F2 = (double*) malloc(3 * N * sizeof(double));


	char filename2[255];
	sprintf(filename2, "/home/pedro/Documents/VMD_files/test.lammpstrj");
	save_lammpstrj(filename2, x, v, N, L, 0);
	double dl = set_pos(x, N, L);
	double cinetica0 = set_vel(v, N, T_gauss);
	double potencial0, temp0;
	potencial0 = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);

	temp0 = cinetica0 * 2.0 / (3.0 * (double)N);


	int i, k, t, , pasos = 100000;
	int p, loops_T = 1000;
	double T0 = temp0, T = T0, dT, Tf = 0.728;
	int tfinal = 100, pasos_T = loops_T * tfinal;
	double  *potencial;
	potencial = (double*) malloc((pasos_T + pasos)* sizeof(double));
	double  *cinetica;
	cinetica = (double*) malloc((pasos_T + pasos) * sizeof(double));
	double h = 0.001;
	float va;
	dT = (Tf - T0) / (double)(loops_T - 1);
	for (p = 0; p < loops_T; p++)
	{
		for (t = 0; t < tfinal; t++)
		{
			va = (float)(t + 1 + tfinal * p) * 100.0 / (float)(pasos_T + pasos);
			printf("Progreso %f", va);
			printf("%%\r");
			*(cinetica + p * tfinal + t) = 0.0;
			position_verlet(x, v, N, h, F);
			apply_PBC(x, N, L);
			*(potencial + p * tfinal + t) = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
			velocity_verlet(v, N, h, F, F2);
			for(i = 0; i < N; i++)
			{
				for(k = 0; k < 3; k++)
				{
					*(cinetica + p * tfinal + t) += *(v + 3 * i + k) * *(v + 3 * i + k) / 2.0;
				}
			}
			*(potencial + p * tfinal + t) = *(potencial + p * tfinal + t) / (double)N;
			*(cinetica + p * tfinal + t) = *(cinetica + p * tfinal + t) / (double)N;
			if (t % 10 == 0)
			{
				save_lammpstrj(filename2, x, v, N, L, p * tfinal + t + 1);
			}
		}
		T += dT;
		T0 = *(cinetica + p * tfinal + tfinal - 1) * 2.0 / 3.0;
		temp_change(v, N, T0, T);
	}

	for (t = 0; t < pasos + 1; t++)
	{
		va = (float)(pasos_T - 1 + t) * 100.0 / (float)(pasos_T + pasos);
		printf("Progreso %f", va);
		printf("%%\r");
		*(cinetica + pasos_T - 1 + t) = 0.0;
		position_verlet(x, v, N, h, F);
		apply_PBC(x, N, L);
		*(potencial + pasos_T - 1 + t) = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
		velocity_verlet(v, N, h, F, F2);
		for(i = 0; i < N; i++)
		{
			for(k = 0; k < 3; k++)
			{
				*(cinetica + pasos_T - 1 + t) += *(v + 3 * i + k) * *(v + 3 * i + k) / 2.0;
			}
		}
		*(potencial + pasos_T - 1 + t) = *(potencial + pasos_T - 1 + t) / (double)N;
		*(cinetica + pasos_T - 1 + t) = *(cinetica + pasos_T - 1 + t) / (double)N;
		if (t % 15 == 0)
		{
			save_lammpstrj(filename2, x, v, N, L, pasos_T - 1 + t);
		}
	}

	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/test_datos.txt");
	fp = fopen(filename, "w");
	fprintf(fp, "%d\t", 0);
	fprintf(fp, "%lf\t", potencial0 / (double)N);
	fprintf(fp, "%lf\t", cinetica0 / (double)N);
	fprintf(fp, "%lf\t", (cinetica0 + potencial0) / (double)N);
	fprintf(fp, "%lf\n", temp0);
	int n;
	for (n = 0; n < pasos_T + pasos ; n++)
	{
		fprintf(fp, "%d\t", n);
		fprintf(fp, "%lf\t", *(potencial + n));
		fprintf(fp, "%lf\t", *(cinetica + n));
		fprintf(fp, "%lf\t", *(cinetica + n) + *(potencial + n));
		fprintf(fp, "%lf\n", *(cinetica + n) * 2.0 / 3.0);
	}
	fclose(fp);
	free(cinetica);
	free(potencial);
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

