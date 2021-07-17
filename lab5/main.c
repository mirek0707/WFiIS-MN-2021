#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define n 7
#define IT_MAX 12

void iloczynMacWek(double macierz[n][n], double wektor[n], double wynik[n]) 
{
    for(int i=0;i<n;i++)
    {
        wynik[i]=0.0;
        for(int j=0;j<n;j++) 
        {
            wynik[i]+=macierz[i][j]*wektor[j]; 
        }
    }
}

double iloczynSkalarny(double wektor1[n], double wektor2[n]) 
{
    double wynik=0.0;
    for(int i=0;i<n;i++) 
    {
        wynik+=wektor1[i]*wektor2[i];
    }
    return wynik;
}

void iloczynMacMac(double wynik[n][n], double macierz1[n][n], double macierz2[n][n]) {
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            wynik[i][j] = 0.0; 
            for (int k = 0; k < n; k++)
                wynik[i][j]+=macierz1[i][k] * macierz2[k][j];
        }
    }
}

void iloczynTensorowy(double result[n][n], double vector[n], double lambda) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
            result[i][j] = lambda * vector[i] * vector[j];
        }
    }
}

void roznicaMacierz(double result[n][n], double macierz1[n][n], double macierz2[n][n]) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
            result[i][j] = macierz1[i][j] - macierz2[i][j];
        }
    }   
}

void copyMatrix(double macierz1[n][n], double macierz2[n][n]) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
            macierz1[i][j] = macierz2[i][j];
        }
    }
}

void transponuj(double macierz1[n][n], double macierz2[n][n]) 
{
    for(int i = 0; i < n; i++) 
    {
        for(int j = 0; j < n; j++)
        {
            macierz1[j][i] = macierz2[i][j];
        }
    }
}

int main()
{
	double A[n][n];
    double W[n][n];
    double Wn[n][n];
    double X[n][n];
    double X_T[n][n];
    double D[n][n];
    double x_old[n];
    double x_new[n];
    double lambda;
    double temp[n][n];
     for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = (1.0) / sqrt(2.0 + fabs(i - j));
            W[i][j] = A[i][j];
            Wn[i][j]=A[i][j];
        }
    }

    FILE *w_wlasne;
    FILE *macierz_D;
    w_wlasne = fopen("wynik.dat", "w");
    macierz_D = fopen("macierz_D.dat", "w");

    for (int k = 0; k < n; k++)
    {
    	for (int i=0;i<n;i++)
    		x_old[i]=1.0;
    	for(int i = 1; i <= IT_MAX; i++)
    	{
    		iloczynMacWek(W, x_old, x_new);
    		lambda=iloczynSkalarny(x_new,x_old)/iloczynSkalarny(x_old,x_old);
    		for (int j=0;j<n;j++)
            {
                x_old[j] = x_new[j]/sqrt(iloczynSkalarny(x_new, x_new));
            }

    		fprintf(w_wlasne, "%d  %f \n", i, lambda);
    	}
    	fprintf(w_wlasne, "\n");

    	iloczynTensorowy(temp, x_old, lambda);
        roznicaMacierz(Wn, W, temp);
        copyMatrix(W, Wn);
        for (int i = 0; i < n; i++)
        {
            X[i][k] = x_old[i];
        }
    }

    iloczynMacMac(temp,A,X);
    transponuj(X_T,X);
    iloczynMacMac(D,X_T,temp);

    for(int i = 0; i < n; i++) 
    {
        for(int j = 0; j < n; j ++) 
        {
            fprintf(macierz_D, "%e ", D[i][j]);
        }
        fprintf(macierz_D, "\n");
    }
    fclose(w_wlasne);
    fclose(macierz_D);

	return 0;
}