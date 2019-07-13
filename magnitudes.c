#include "magnitudes.h"

// f(x) q calcula y devuelve el Cv/(k*N) con k = 1
double C_v(double *cinetica, int pasos, int N)
{
	int t;
	double cv, T = 0.0, E_k = 0.0;
	for (t = 0; t < pasos; t++)
	{
		T += (double)*(cinetica + t) * 2.0 / 3.0;
	}
	T = T / (double)pasos;
	for (t = 0; t < pasos; t++)
	{
		E_k += 1.0 / (double) *(cinetica + t); //promedio de la inversa de la cinética
	}
	E_k = E_k / (double)pasos;

	cv = 1.0 / ((double)N - (double)N * T * ((double)N * 3.0 / 2.0 - 1.0) * E_k);
	return cv;
}
/*long double C_v(double *cinetica, int pasos, int N)
{
	int t;
	long double cv, T = 0.0, E_k = 0.0;
	for (t = 0; t < pasos; t++)
	{
		T += (long double)*(cinetica + t) * 2.0 / 3.0;
	}
	T = T / (long double)pasos;
	for (t = 0; t < pasos; t++)
	{
		E_k += 1.0 / (long double) *(cinetica + t); //promedio de la inversa de la cinética
	}
	E_k = E_k / (long double)pasos;

	cv = 1.0 / ((long double)N - (long double)N * T * ((long double)N * 3.0 / 2.0 - 1.0) * E_k);
	return cv;
}*/
