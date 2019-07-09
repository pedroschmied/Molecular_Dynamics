#include "termalizacion.h"
#include "general.h"
#include "avanzar.h"

//#define PI  3.14159
double lambda(double *x, double dl, int N)
{
	double l = 0.0;
	int i;
	for (i = 0; i < 3 * N; i++)
	{
		l += cos((*(x + i) - dl / 2.0) * 2.0 * PI / dl);
	}
	l = l / (3.0 * (double)N); //promedio en las tres coordenadas lx, ly y lz
	return l;
}

double H_boltzmann (double *v, int N)
{
	int i, k;
	int columnas = 20;
	double *minimo, *maximo, *hist, *teo_H; // calculo los límites del histograma
	minimo = (double*) malloc(3 * sizeof(double));
	maximo = (double*) malloc(3 * sizeof(double));
	teo_H = (double*) malloc(3 * sizeof(double));
	hist = (double*) malloc(columnas * sizeof(double));
	for (k = 0; k < 3; k++)
	{
		*(minimo + k) = *(v + 5);
		*(maximo + k) = *(v + 5);
		for (i = 0; i < N; i++)
		{
			if(*(v + 3 * i + k) > *(maximo + k))
			{
				*(maximo + k) = *(v + 3 * i + k);
			}
			else if(*(v + 3 * i + k) < *(minimo + k))
			{
				*(minimo + k) = *(v + 3 * i + k);
			}
		}
	}
	double dh;
	int h;
	for(k = 0; k < 3; k++)
	{
		*(teo_H + k) = 0.0;
		dh = (*(maximo + k) - *(minimo + k)) / (double) (columnas - 1);
		for (i = 0; i < columnas; i++)
		{
			*(hist + i) = 0.0;
		}
		for (i = 0; i < N; i++)
		{
			h = (int)((*(v + 3 * i + k) - *(minimo + k)) / dh);

			*(hist + h) += 1.0 / (double)N;
		}
		for (i = 0; i < columnas; i++)
		{
			if (*(hist + i) != 0)
			{
				*(teo_H + k) += *(hist + i) * log(*(hist + i)) * dh;
			}
		}
	}
	double H = 0.0;
	for (k = 0; k < 3; k++)
	{
		H += *(teo_H + k) / 3.0;
	}
	free(hist);
	free(minimo);
	free(maximo);
	free(teo_H);
	return H;
}

double distribucion_radial(double *x, double *g_r, double dr_g, double L, int N)
{
	double rij, r;
	int i, j, h;
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));

	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			delta_x(x, i, j, delta_r);
			fuerza_PCB(delta_r, L);
			rij = norma2(delta_r);
			rij = sqrt(rij);
			if(rij <= L / 2.0)
			{
				h = (int)(rij / dr_g);
				r = (double) h * dr_g;	//como el int redondea para abajo rij > r siempre o igual
				if(rij - r > dr_g / 2.0)
				{				//estoy haciendo esto:
					h += 1;		// if(r-dr_g/2 < rij > r +dr_g/2)   {delta(r-rij) = 1}
				}
				*(g_r + h) += 1.0;
			}
		}
	}
	free(delta_r);
	return 0.0;
}
