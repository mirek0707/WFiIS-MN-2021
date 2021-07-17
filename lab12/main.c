#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double k, double m, double x)
{
    return pow(x, m)*sin(k*x);
}

double silnia(int n)
{
    if (n == 0 || n == 1) 
    	return 1;
    else  
    	return (silnia(n-1)*n);
}

double suma(double i, double m, double k, double x)
{
    double licznik=pow(-1, i)*pow(k*x, 2*i+m+2);
    double mianownik=pow(k, m+1)*silnia(2*i+1)*(2*i+m+2);
    return licznik/mianownik;
}

void simspson(double n, double m, double k,  double I, FILE *p)
{

    double h=M_PI/(n-1);
    double f0=f(k, m, 0);
    double fn=f(k, m, M_PI);

    int i=0;
    double suma=0.0;
    double C=0.0;

    for (double j=h; j<M_PI; j+=h)
    {
        if (i%2==1)
            suma+=2*f(k, m, j);
        else
            suma+=4*f(k, m, j);  
        i++;
    }

    suma=suma+f0+fn;
    C=(h/3)*suma;
    fprintf(p, "%f \t %f \t %f \n", n, C, fabs(C - I));
}

void policz(double m, double k, double I,  FILE *file)
{
    double a=0.0;
    double b=M_PI;
    double temp=0.0;
    double wynik, modul;
    int max=30;

    for (int i=0; i<max; i++)
    {
        temp+=suma(i, m, k, b)-suma(i, m, k, a);
        wynik=temp;
        modul=fabs(wynik-I);
        fprintf(file, "%d \t %f \t %f \n", i+1, wynik, modul);
    }
}

void policzS(double m, double k, double I,  FILE *file)
{
    double N[5] = {11, 21, 51, 101, 201};
    for (int j=0; j<5; j++)
    {
        simspson(N[j], m, k, I, file);
    }
}

int main()
{

    FILE *d1=fopen("d1.txt", "w");
	policz(0.0, 1.0, 2.0, d1);

    FILE *d2=fopen("d2.txt", "w");
    policz(1.0, 1.0, M_PI, d2);

    double I = 56.363569;
    FILE *d3 = fopen("d3.txt", "w");
    policz(5.0, 5.0, I, d3);

    FILE *dS1=fopen("dS1.txt", "w");
    FILE *dS2=fopen("dS2.txt", "w");
    FILE *dS3=fopen("dS3.txt", "w");
    
    policzS(0.0, 1.0, 2.0, dS1);
    policzS(1.0, 1.0, M_PI, dS2);
    policzS(5.0, 5.0, I, dS3);

    fclose(d1);
    fclose(d2);
    fclose(d3);
    fclose(dS1);
    fclose(dS2);
    fclose(dS3);

    return 0;
}