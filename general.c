float aleatorio()
{
	float x;
	x = (float) rand() / (float) RAND_MAX;
	return x;
}

float gaussiana(float mu, float sigma)
{
	int i;
	float z = 0.0, n = 10.0, x;
	for (i = 0; i < n; i++)
	{
		x = aleatorio();
		z+= x;
	}
	z = sqrt(12.0) * (z/n - 0.5);
	float g = z * sigma + mu;
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
