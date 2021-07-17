#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include<math.h>
#include<stdio.h>

const double xmin=-5.0;
const double xmax=5.0;
const double alfa=0.0;
const double beta=0.0;

double f1(double x)
{
    return (1/(1+(pow(x,2))));
}

double f2(double x)
{
    return cos(2*x);
}

double deltax(int n)
{
    return (xmax-xmin)/(double)(n-1);
}

double di(double h, double h1, double y, double y1, double y2)
{
    return (6.0/(h+h1))*((y2-y1)/h1-(y1-y)/h);
}

void wyznacz_M(double *x, double *y, double *m, int n, double alfa, double beta)
{
	gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *d = gsl_vector_calloc(n);
    gsl_vector_set(d, 0, alfa);
    gsl_vector_set(d, n-1, beta);

    double lambda[n];
    double mi[n];

    for (int i=1; i<n-1; i++)
    {
    	gsl_vector_set(d, i, di((x[i]-x[i-1]), (x[i+1]-x[i]), y[i-1], y[i], y[i+1]));
    }

    for (int i=1; i<n-2; i++)
    {
    	lambda[i]=(x[i+1]-x[i])/(double)((x[i]-x[i-1])+(x[i+1]-x[i]));
    	mi[i]=1-lambda[i];
    }

    gsl_matrix_set(A, 0, 0, 1.0);
    gsl_matrix_set(A, n - 1, n - 1, 1.0);
    for (int i=1; i<n-1; i++)
    {
        gsl_matrix_set(A, i, i, 2.0);
    }
    for (int i=2; i<n; i++)
    {
        gsl_matrix_set(A, i - 1, i, lambda[i-1]);//////////////-2
    }

    for (int i=0; i<n-2; i++)
    {
        gsl_matrix_set(A, i + 1, i, mi[i+1]);//////////////0
    }

    gsl_linalg_HH_svx(A, d);

    for (int i=0; i<n; i++)
    {
        m[i] = gsl_vector_get(d, i);
    }

    gsl_vector_free(d);
    gsl_matrix_free(A);

}

double wyznacz_Sx(double *xw, double *yw, double *m, int n, double x)
{
	double Sx=0;
	int i;
    for (int j=1; j<n; j++)
    {
        if (x>=xw[j-1] && x<=xw[j])
        {
            i = j;
            break;
        }
    }
    double h = xw[i]-xw[i-1];
    double A = ((yw[i]-yw[i-1])/h)-(h/6)*(m[i]-m[i-1]);
    double B = yw[i-1]-(m[i-1]*pow(h,2))/6;
    Sx = m[i-1]*(pow((xw[i]-x), 3)/(6.0*h)) + m[i]*(pow((x-xw[i-1]),3)/(6.0*h)) + A*(x-xw[i-1])+B;

    return Sx;
}

void oblicz(int n, FILE *p1, FILE *p2)
{
    double x1[n], y1[n], m1[n], x2[n], y2[n], m2[n];

    for(int i=0; i<n; i++)
    {
    	x1[i] = x2[i]=y1[i]=y2[i]=0;
    	x1[i] = xmin+deltax(n)*i;
    	x2[i] = xmin+deltax(n)*i;
    	y1[i] = f1(x1[i]);
    	y2[i] = f2(x1[i]);
    }

    wyznacz_M(x1, y1, m1, n, alfa, beta);
    wyznacz_M(x2, y2, m2, n, alfa, beta);

    for (int i=0; i<=1000; i++)
    {
        double x = xmin+0.01*i;
        double Sx = wyznacz_Sx(x1, y1, m1, n, x);
        fprintf(p1, "%g %g\n", x, Sx);
    }
    for (int i=0; i<=1000; i++)
    {
        double x = xmin+0.01*i;
        double Sx = wyznacz_Sx(x2, y2, m2, n, x);
        fprintf(p2, "%g %g\n", x, Sx);
    }

    fprintf(p1, "\n\n");
    fprintf(p2, "\n\n");
}

int main()
{
	FILE *p1;
    p1 = fopen("f1.txt", "w");
    FILE *p2;
    p2 = fopen("f2.txt", "w");
    FILE *p3;
    p3 = fopen("pochodne.txt", "w");

    oblicz(5, p1, p2);
    oblicz(8, p1, p2);
    oblicz(21, p1, p2);

    int n = 10;
    double sigmax = 0.01;
    double x[n];
    double y[n];
    double m[n];
    double pochodna[n];

    for (int i=0; i<n; i++)
    {
        x[i]=y[i]=pochodna[i]=0.0;
        x[i] = xmin+deltax(n)*i;
        y[i] = f1(x[i]);
        pochodna[i] = (f1(x[i] - sigmax) - 2*f1(x[i]) + f1(x[i]+sigmax))/(pow(sigmax,2));
    }
    wyznacz_M(x, y, m, n, alfa, beta);

    for (int i=0;i<n;i++)
    {
        fprintf(p3, " %g   %g  %g \n", x[i], m[i], pochodna[i]);
    }

    fclose(p1);
    fclose(p2);
    fclose(p3);

    return 0;
}