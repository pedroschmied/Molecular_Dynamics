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
	int N = 216;
	float rho = 0.8442, L = cbrt(N/rho);
	float T =  2.0;
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

	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/MD_datos1.txt");
	fp = fopen(filename, "w");

	char filename2[255];
	sprintf(filename2, "/home/pedro/Documents/VMD_files/test.lammpstrj");

	set_pos(x, N, L);
	double cinetica = set_vel(v, N, T);
	double potencial, temp;
	potencial = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
	temp = cinetica * 2.0 / 3.0 * (double)N;
	fprintf(fp, "%d\t", 0);
	fprintf(fp, "%lf\t", potencial / (double)N);
	fprintf(fp, "%lf\t", cinetica / (double)N);
	fprintf(fp, "%lf\t", (cinetica + potencial) / (double)N);
	fprintf(fp, "%lf\n", temp);

	double h = 0.001;
	int i, k, t, tfinal = 100000;
	float va;
	for (t = 1; t < tfinal + 1; t++)
	{
		va = (float)t * 100.0 / ((float)tfinal);
		printf("Progreso\t");
		printf("%f%%\r", va);
//		printf("%%\r");
		cinetica = 0.0;
		position_verlet(x, v, N, h, F);
		apply_PBC(x, N, L);
		potencial = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
		velocity_verlet(v, N, h, F, F2);
		for(i = 0; i < N; i++)
		{
			for(k = 0; k < 3; k++)
			{
				cinetica += *(v + 3 * i + k) * *(v + 3 * i + k) / 2.0;
			}
		}
		potencial = potencial / (double)N;
		cinetica = cinetica / (double)N;
		temp = cinetica * 2.0 / 3.0;
		fprintf(fp, "%d\t", t);
		fprintf(fp, "%lf\t", potencial);
		fprintf(fp, "%lf\t", cinetica);
		fprintf(fp, "%lf\t", cinetica + potencial);
		fprintf(fp, "%lf\n", temp);
		save_lammpstrj(filename2, x, v, N, L, t-1);
	}
	fclose(fp);
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

