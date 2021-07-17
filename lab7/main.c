#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double xa = -3.0;
const double xb = 8.0;

double f(double x)
{
	return x/(1+pow(x,2));
}

double Lagrange (double xn, int n, double *x, double *y)
{
	int j;
	double wielomian = 0.0, L, M;
	for (int i=0;i<=n;i++)
	{
		L=M=1.0;
		for (j=0;j<=n;j++)
		{
			
			if (i!=j)
			{
				L*=(xn-x[j]);
				M*=(x[i]-x[j]);
			}
			
		}
		wielomian+=y[i]*(L/M);
	}
	return wielomian;
}

double Czebyszew(int n, int i)
{
    return 0.5 * ((xb - xa) * cos(M_PI* (2.0 * (double)i + 1.0) / (2.0 * (double)n + 2.0)) + (xa + xb));
}

void calculate(int n ,FILE *p1)
{
    double x[n + 1];
    double y[n + 1];
    double cx[n + 1];
    double cy[n + 1];

    double deltax=(xb-xa)/(double)n;

    for (int i = 0; i <= n; i++)
    {
        x[i] = y[i] = cx[i] = cy[i] = 0.0;
        x[i] = xa + deltax * i;
        y[i] = f(x[i]);
        cx[i] = Czebyszew(n, i);
        cy[i] = f(cx[i]);
    }

    double k=((xb-xa)/(double)200);
    double wielomian = 0.0;
    double fx;

   
    for(int i=0; i<=200; i++)
    {
        double xn=xa+k*i;
        fx=f(xn);
        wielomian=Lagrange(xn,n,x,y);
        
        fprintf(p1, "%12.6g \t", xn);
        fprintf(p1, "%12.6g \t", fx);
        fprintf(p1, "%12.6g \t", wielomian);
        wielomian=Lagrange(xn,n,cx,cy);
		fprintf(p1, "%12.6g \n", wielomian);
    }
}

int main()
{

    int n_5 = 5;
    int n_10 = 10;
    int n_15 = 15;

    FILE *p1;
    p1 = fopen("p1.txt", "w");
    FILE *p2;
    p2 = fopen("p2.txt", "w");
    FILE *p3;
    p3 = fopen("p3.txt", "w");

    calculate(n_5, p1);

    calculate(n_10, p2);

    calculate(n_15, p3);

    fclose(p1);
    fclose(p2);
    fclose(p3);

    return 0;
}