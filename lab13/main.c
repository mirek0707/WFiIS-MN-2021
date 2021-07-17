#include"nrutil.c"
#include"nrutil.h"
#include"nr.h"
#include"gauher.c"
#include"gauleg.c"
#include"gaulag.c"
#include"gammln.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double fn1(double x)
{
    return log(x);
}

double fn2a(double x)
{
    return pow((x-10.0), 2.0) * sin(4.0 * x);
}

double fn2b(double x)
{
    return pow((x-10.0), 2.0) * sin(4.0 * x) * exp(-x);
}

double fn3a(double x)
{
    return pow(x, 7.0) * pow(2.0, (-pow(x, 2.0) +x+4.0));
}

double fn3b(double x)
{
    return pow(x, 7.0) * pow(2.0, (-pow(x, 2.0) +x+4.0)) * exp(-pow(x, 2.0));
}

int main()
{
	/////////1
	double a=10.0;
    double b=0.0;
    FILE *plik1 = fopen("plik1.txt", "w");
    double Cnum = a * log(a) - a;

    for(int n=5; n<=70; n++)
    {
        double Cdok=0.0;
        float *x = vector(1, n);
        float *w = vector(1, n);
        gauleg(b, a, x, w, n);
        for(int i=1; i<=n; i++)
        {
            Cdok+=w[i] * fn1(x[i]);
        }
        fprintf(plik1, " %d \t %f \t %e \n", n, Cdok, fabs(Cdok - Cnum));
    }

    /////////2a
    FILE *plik2a = fopen("plik2a.txt", "w");
    Cnum = 22.95461022;

    for(int n=5; n<=70; n++)
    {
        double Cdok = 0.0;
        float *x = vector(1, n);
        float *w = vector(1, n);
        gaulag(x, w, n, 0.0);
        for(int i=1; i<=n; i++)
        {
            Cdok+=w[i] * fn2a(x[i]);
        }
        fprintf(plik2a, " %d \t %f \t %e \n", n, Cdok, fabs(Cdok - Cnum));
    }

    /////////2b
    FILE *plik2b = fopen("plik2b.txt", "w");
    Cnum = 22.95461022;

    for(int n=5; n<=70; n++)
    {
        double Cdok = 0.0;
        float *x = vector(1,n);
        float *w = vector(1,n);
        gauleg(b, a, x, w, n);
        for(int i=1; i<=n; i++)
        {
            Cdok+=w[i] * fn2b(x[i]);
        }
        fprintf(plik2b, " %d \t %f \t %e \n", n, Cdok, fabs(Cdok - Cnum));
    }

    /////////3a
    FILE *plik3a = fopen("plik3a.txt", "w");
    Cnum = 14.83995751;

    for(int n=5; n<=70; n++)
    {
        double Cdok = 0.0;
        float *x = vector(1,n);
        float *w = vector(1,n);
        gauher(x, w, n);
        for(int i=1; i<=n; i++)
        {
            Cdok+=w[i] * fn3a(x[i]);
        }
        fprintf(plik3a, " %d \t %f \t %e \n", n, Cdok, fabs(Cdok - Cnum));
    }

    /////////3b
    FILE *plik3b = fopen("plik3b.txt", "w");
    Cnum = 14.83995751;
    double x1 = -10.0;
    double x2 = 15.0;

    for(int n=5; n<=70; n++)
    {
        double Cdok=0.0;
        float *x = vector(1,n);
        float *w = vector(1,n);
        gauleg(x1, x2, x, w, n);
        for(int i=1; i<=n; i++)
        {
            Cdok+=w[i] * fn3b(x[i]);
        }
        fprintf(plik3b, " %d \t %f \t %e \n", n,  Cdok, fabs(Cdok - Cnum));
    }
    ///////////////////
    FILE * other1, * other2, * other3;
    other1 = fopen("fun1.dat","w");
    other2 = fopen("fun2.dat","w");
    other3 = fopen("fun3.dat","w");

    for (double i = 0.1; i <=10; i+=0.1) 
        fprintf(other1,"%f\t%f\n",i, fn1(i));

    for (double i = 0; i <=10; i+=0.1) 
        fprintf(other2,"%f\t%f\n",i, fn2b(i));
        

    for (double i = -10; i <=15; i+=0.1) 
        fprintf(other3,"%f\t%f\n",i, fn3b(i));
    fclose(plik1);
    fclose(plik2a);
    fclose(plik2b);
    fclose(plik3a);
    fclose(plik3b);
    fclose(other1);
    fclose(other2);
    fclose(other3);
    return 0;
}
