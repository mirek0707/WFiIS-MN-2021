#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double f (double x)
{
	return sin(x)-pow(x,2)/8;
}
double fpochodna (double x)
{
	return cos(x)-x/4;
}
void sieczne(double x0, double x1, int IT_MAX, FILE* file)
{
	double x2;
	for (int k=1; k<=IT_MAX; k++)
	{
		x2=x1-(f(x1)*(x1-x0))/(f(x1)-f(x0));
		fprintf(file, "%2d  %15f  %15f %15f \n", k, x2, f(x1), f(x0));
		printf("%2d  %15f  %15f %15f \n", k, x2, f(x1), f(x0));
		x0=x1;
		x1=x2;
	}
	fprintf(file, "\n");
	printf("\n");
}
void Newton(double x, int IT_MAX, FILE* file)
{
	for (int k=1; k<=IT_MAX; k++)
	{
		x=x-f(x)/fpochodna(x);
		fprintf(file, "%2d  %15f  %15f %15f \n", k, x, f(x), fpochodna(x));
		printf("%2d  %15f  %15f %15f \n", k, x, f(x), fpochodna(x));
	}
	fprintf(file, "\n");
	printf("\n");
}
int main()
{
	FILE *f1;
    f1 = fopen("Newton.txt", "w");
    FILE *f2;
    f2 = fopen("Sieczne.txt", "w");

    Newton(-8, 10, f1);
    Newton(8, 10, f1);

    sieczne(-8, -8.1, 15, f2);
    sieczne(8, 8.1, 15, f2);

    fclose(f1);
    fclose(f2);
	return 0;
}