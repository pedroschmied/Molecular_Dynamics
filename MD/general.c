#include "general.h"
double aleatorio()
{
	double x;
	x = (double) rand() / (double) RAND_MAX;
	return x;
}

float gaussiana(double mu, double sigma)
{
	int i;
	double z = 0.0, n = 10.0, x;
	for (i = 0; i < n; i++)
	{
		x = aleatorio();
		z+= x;
	}
	z = sqrt(12.0 * n) * (z / n - 0.5);
	double g = z * sigma + mu;
	return g;
}
// calculo la el r cuadrado de la partícula i
double norma2(double *x)
{
	int k = 0;
	double r2 = 0.0;
	for (k = 0; k < 3; k++)
	{
		r2 += *(x + k) * *(x + k);
	}
	return r2;
}
int delta_x(double *x, int i, int j, double *delta_r)
{
	int k;
	for (k = 0; k < 3; k++)
	{
		*(delta_r + k) = *(x + 3 * i + k) - *(x + 3 * j + k);
	}
	return 0;
}

double mean (double *v, int o, int n, int k)
{
	int i;
	double x = 0.0;
	for (i = o; i < (n + o);  i = i + k)
	{
		x = x + *(v + i);
	}
	x = x * (double)k / (double) n;
	return x;
}

double mean2 (double *v, int o, int n, int k)
{
	int i;
	double x = 0.0;
	for (i = o; i < (n + o);  i = i + k)
	{
		x = x + *(v + i) * *(v + i);
	}
	x = x * (double)k / (double) n;
	return x;
}

double std2 (double *v, int o, int n, int k)
{
	double x, prom, prom2;
	prom = mean(v, o, n, k);
	prom2 = mean2(v, o, n, k);
	x = prom2 - prom * prom;
	return x;
}
/*
long double mean (double *v, int o, int n, int k)
{
	int i;
	long double x = 0.0;
	for (i = o; i < (n + o);  i = i + k)
	{
		x = x + (long double)*(v + i);
	}
	x = x * (long double)k / (long double) n;
	return x;
}

long double mean2 (double *v, int o, int n, int k)
{
	int i;
	long double x = 0.0;
	for (i = o; i < (n + o);  i = i + k)
	{
		x = x + (long double)(*(v + i) * *(v + i));
	}
	x = x * (long double)k / (long double) n;
	return x;
}

long double std2 (double *v, int o, int n, int k)
{
	long double x, prom, prom2;
	prom = mean(v, o, n, k);
	prom2 = mean2(v, o, n, k);
	x = prom2 - prom * prom;
	return x;
}*/
