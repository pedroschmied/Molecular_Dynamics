#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"
#include "termalizacion.h"

#define PI  3.14159

int main()
{
	int N = 512;
	float rho = 0.8442, L = cbrt(N / rho);
	float T_gauss =  0.4;
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


	char filename2[255];
	sprintf(filename2, "/home/pedro/Documents/VMD_files/test.lammpstrj");
//condiciones iniciales
//	double dl = 
	set_pos(x, N, L);
//	double cinetica0 = 
	set_vel(v, N, T_gauss);
//	double potencial0, temp0;
//	potencial0 = 
	fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
//	temp0 = cinetica0 * 2.0 / (3.0 * (double)N);
//	double lambda0 = coeficiente_verlet(x, dl, N);
//	double H0 = H_boltzmann (v, N);

	save_lammpstrj(filename2, x, v, N, L, 0);

	int i, t, c;
	int pasos = 5000, termalizacion = 2000, correlacion = 500;
	double  *potencial;
	potencial = (double*) malloc((pasos)* sizeof(double));
	double  *cinetica;
	cinetica = (double*) malloc((pasos) * sizeof(double));
//	double *lambda, *H;
//	lambda = (double*) malloc((pasos) * sizeof(double));
//	H = (double*) malloc((pasos) * sizeof(double));

	double h = 0.001;
	float va;
	for (t = 0; t < termalizacion; t++)
	{
		va = (float)(t) * 100.0 / (float)(pasos * correlacion + termalizacion);
		printf("Progreso %.2f", va);
		printf("%%\r");
		position_verlet(x, v, N, h, F);
		apply_PBC(x, N, L);
		fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
		velocity_verlet(v, N, h, F, F2);
	}
	for (t = 0; t < pasos; t++)
	{
		*(cinetica + t) = 0.0;
		*(potencial + t) = 0.0;
		for(c = 0; c < correlacion; c++)
		{
			va = (float)(t * correlacion + c + termalizacion) * 100.0 / (float)(pasos * correlacion + termalizacion);
			printf("Progreso %.2f", va);
			printf("%%\r");

			position_verlet(x, v, N, h, F);
			apply_PBC(x, N, L);
			*(potencial + t) += fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
			velocity_verlet(v, N, h, F, F2);
			for(i = 0; i < 3 * N; i++)
			{
				*(cinetica + t) += *(v + i) * *(v + i) / 2.0;
			}
		}
		*(potencial + t) = *(potencial + t) / ((double)N * correlacion);
		*(cinetica + t) = *(cinetica + t) / ((double)N * correlacion);

		save_lammpstrj(filename2, x, v, N, L, t + 1);
//		*(lambda + t) = coeficiente_verlet(x, dl, N);
//		*(H + t) = H_boltzmann (v, N);
	}
	double temperatura = 0.0;
	for(t = 0; t < pasos; t++)
	{
		temperatura += *(cinetica + t) * 2.0 / 3.0;
	}
	temperatura = temperatura / pasos;
	printf("\n%lf\n", temperatura);

	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/MD_datos_1_a_T_%.2lf.txt", T_gauss);
	fp = fopen(filename, "w");
	int n;
	for (n = 0; n < pasos; n++)
	{
		fprintf(fp, "%d\t", n);
		fprintf(fp, "%lf\t", *(potencial + n));
		fprintf(fp, "%lf\t", *(cinetica + n));
		fprintf(fp, "%lf\t", *(cinetica + n) + *(potencial + n));
		fprintf(fp, "%lf\n", *(cinetica + n) * 2.0 / 3.0);
//		fprintf(fp, "%lf\t", *(lambda + n));
//		fprintf(fp, "%lf\n", *(H + n));
	}
	fclose(fp);
//	free(lambda);
//	free(H);
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
#include "termalizacion.c"
