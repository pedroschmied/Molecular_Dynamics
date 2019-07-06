#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "general.h"
#include "inicializar.h"
#include "visualizacion.h"
#include "interaccion.h"
#include "avanzar.h"
#include "termalizacion.h"

int main()
{
	int N = 512;
	float rho = 0.8442, L = cbrt(N / rho);
	float T_gauss =  0.1;
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
	double dl = set_pos(x, N, L);
	double cinetica0 = set_vel(v, N, T_gauss);
	double potencial0, temp0;
	potencial0 = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
	temp0 = cinetica0 * 2.0 / (3.0 * (double)N);
	save_lammpstrj(filename2, x, v, N, L, 0);

	int i, t;
	int pasos = 3000;
	double  *potencial;
	potencial = (double*) malloc((pasos)* sizeof(double));
	double  *cinetica;
	cinetica = (double*) malloc((pasos) * sizeof(double));
	double h = 0.001;
	float va;
	double l, H;
	FILE * fp3;
	char filename3[500];
	sprintf (filename3,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/termalizacion%g.txt", T_gauss);
	fp3 = fopen(filename3, "w");

	for (t = 0; t < pasos; t++)
	{
		va = (float)(t) * 100.0 / (float)(pasos);
		printf("Progreso %f", va);
		printf("%%\r");
		*(cinetica + t) = 0.0;
		position_verlet(x, v, N, h, F);
		apply_PBC(x, N, L);
		*(potencial + t) = fuerzas(tabla_F, tabla_V, F, F2, x, rc2, dr2, N, L);
		velocity_verlet(v, N, h, F, F2);
		for(i = 0; i < 3 * N; i++)
		{
			*(cinetica + t) += *(v + i) * *(v + i) / 2.0;
		}
		*(potencial + t) = *(potencial + t) / (double)N;
		*(cinetica + t) = *(cinetica + t) / (double)N;
		if (t % 15 == 0)
		{
			save_lammpstrj(filename2, x, v, N, L, t + 1);
		}
		l = lambda(x, dl, N);
		H = H_boltzmann (v, N);
		fprintf(fp3, "%d\t", t);
		fprintf(fp3, "%lf\t", l);
		fprintf(fp3, "%lf\n", H);
	}
	fclose(fp3);
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/MD_datos.txt");
	fp = fopen(filename, "w");
	fprintf(fp, "%d\t", 0);
	fprintf(fp, "%lf\t", potencial0 / (double)N);
	fprintf(fp, "%lf\t", cinetica0 / (double)N);
	fprintf(fp, "%lf\t", (cinetica0 + potencial0) / (double)N);
	fprintf(fp, "%lf\n", temp0);
	int n;
	for (n = 0; n < pasos ; n++)
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
#include "termalizacion.c"
