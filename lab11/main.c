#include <stdio.h>
#include <math.h>
#include <time.h>
#include </usr/include/gsl/gsl_math.h>
#include </usr/include/gsl/gsl_linalg.h>
#include </usr/include/gsl/gsl_errno.h>
#include </usr/include/gsl/gsl_fft_complex.h>

double f(int N, int i)
{
    double omega = ((4*M_PI)/(N));
    return sin(omega*i) + sin(2*omega*i) + sin(3*omega*i);
}

double fszum(int N, int i)
{
	double szum = 2*(rand()/(RAND_MAX+1.0))-1;
	return (f(N,i)+szum);
}

void wypelnij(double *y, int n, FILE *p)
{
    for(int i=0; i<n; i++)
    {
        y[2*i] = fszum(n, i);
        y[2*i + 1] = 0.0;
        fprintf(p, "%d %e \n",i, y[2*i]);
    }
    fprintf(p, "\n\n");
}

void szukajMax(double *y, int n, double max, int l)
{
    double modul=0.0;
    for(int i=0; i<n; i++)
    {
        modul = sqrt((pow(y[2*i],2) +pow(y[2*i+1],2)));
        if(max < modul)
        {
            max = modul;
        }
    }
    for(int i=0; i<n; i++)
    {
        modul = sqrt((pow(y[2*i],2)+pow(y[2*i+1],2)));
        if(modul < (max/2))
        {
        	y[2*i] = 0.0;
        	y[2*i + 1] = 0.0;            
        }
    }
    printf("max dla k = %d: %g \n",l, max);
}

void norma(double *y, int n, FILE *p)
{
    for(int i = 0; i < n; i++)
    {
        y[2*i]/=n;
        y[2*i + 1]/=n;
        fprintf(p, "%d %e \n",i, y[2*i]);         
    }
}
int main()
{
	FILE *p1 = fopen("y8.txt", "w");
    FILE *p2 = fopen("y10.txt", "w");
    FILE *p3 = fopen("y12.txt", "w");
    FILE *p4 = fopen("fft8.txt", "w");

    const int k1 = 8;
    const int k2 = 10;
    const int k3 = 12;
    const int n1 = pow(2,k1);
    const int n2 = pow(2,k2);
    const int n3 = pow(2,k3);

    double y1[2*n1];
    double y2[2*n2];
    double y3[2*n3];

    //1
    wypelnij(y1, n1, p1);
    wypelnij(y2, n2, p2);
    wypelnij(y3, n3, p3);

    //2
    gsl_fft_complex_radix2_forward(y1, 1, n1);
    gsl_fft_complex_radix2_forward(y2, 1, n2);
    gsl_fft_complex_radix2_forward(y3, 1, n3);

    //3
    for(int i=0; i<n1; i++)
    {
        fprintf(p4, "%d %e %e\n",i, y1[2*i], y1[2*i+1]);
    }

    //4
    double max1=0.0;
    double max2=0.0;
    double max3=0.0;
    szukajMax(y1, n1, max1, k1);
    szukajMax(y2, n2, max2, k2);
    szukajMax(y3, n3, max3, k3);

    //5
    gsl_fft_complex_radix2_backward(y1, 1, n1);
    gsl_fft_complex_radix2_backward(y2, 1, n2);
    gsl_fft_complex_radix2_backward(y3, 1, n3);
    norma(y1, n1, p1);
    norma(y2, n2, p2);
    norma(y3, n3, p3);

    fclose(p1);
    fclose(p2);
    fclose(p3);
    fclose(p4);

    return 0;
}