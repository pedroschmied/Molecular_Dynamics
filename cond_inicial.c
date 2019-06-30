#include <stdlib.h>
#include <stdio.h>

#include "inicializar.h"
#include "visualizacion.h"
#include <unistd.h>
#include "interaccion.h"
#include "avanzar.h"


int main()
{
	int N = 512;
	float rho = 0.8442, L = cbrt(N/rho);
	float T =  2.0;
	int i;
	double *x, *v;
	x = (double*) malloc(3 * N * sizeof(double));
	v = (double*) malloc(3 * N * sizeof(double));
	set_pos(x, N, L);
	double E_cin = set_vel(v, N, T);
	printf("%lf", 2.0 * E_cin / (3.0 * (double)N));
	FILE * fp;
	char filename[500];
	sprintf (filename,"/home/pedro/Desktop/Universidad/Fisica_computacional/Datos_molecular_dynamics/MD/velocidad.txt");
	fp = fopen(filename, "w");

/*	for(i = 0; i < N; i++)
	{
		for(int k = 0; k < 3; k++)
		{
			E_cin = 
		}
	}
*/	for(i = 0; i < N; i++)
	{
		fprintf(fp, "%d\t", i);
		fprintf(fp, "%f\n", *(v + 3 * i + 0) * *(v + 3 * i + 0) + *(v + 3 * i + 1) * *(v + 3 * i + 1) + *(v + 3 * i + 2) * *(v + 3 * i + 2));
	}
	fclose(fp);

/*	char filename2[255];
	sprintf(filename2, "/home/pedro/Documents/VMD_files/test.lammpstrj");
	int N_frames = 100;
	float t = 0, dt = 0.1;

	for(int l = 0; l < N_frames; l++)
	{
		for(int i = 0; i < 3* N; i++)
		{
			*(x + i) += *(v + i)*t;
			*(x + i) -= ((int)(*(x + i)/L) % (int)L) * L;
		}
		t += dt;
		save_lammpstrj(filename2, x, v, N, L, l);
	}
*/
//	printf("Energía cinética = %f", E_cin);
	free(x);
	free(v);
	return 0;
}
#include "general.c"
#include "inicializar.c"
#include "visualizacion.c"
#include "interaccion.c"
#include "avanzar.c"

