#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define n 1000

void iloczyn (double a[n][n], double x[n], double temp[n])
{
	for (int i=0;i<n;i++)
	{
		temp[i]=0.0;
		for (int j=fmax(0,i-10); j<=fmin(i+10,n-1);j++)
		{
			temp[i] += a[i][j] * x[j];
		}
	}
}
double iloczynsk(double w1[n], double w2[n])
{
	double wynik = 0.0;
	for (int i=0;i<n;i++)
	{
		wynik+=w1[i]*w2[i];
	}
	return wynik;
}
void sumaw(double w1[n], double w2[n])
{
	for (int i=0;i<n;i++)
	{
		w2[i]=w1[i]+w2[i];
	}
}
int main ()
{
	int m=10;
	double A[n][n];
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			if (abs(i-j)<=m)
				A[i][j]=1.0/(1.0+abs(i-j));
			else
				A[i][j]=0.0;
		}
	}
	double b[n];
	for (int i=0;i<n;i++)
	{
		b[i]=i;
	}
	int k=0;
	double x[n];
	for (int i=0;i<n;i++)
	{
		x[i]=0.0;
	}
	double temp[n], r[n], alfa;
	FILE *p1 = fopen("zad a.txt", "w");
	do
	{
		k++;
		iloczyn(A, x, temp);
		for(int i=0;i<n;i++)
		{
			r[i]=b[i]-temp[i];
		}
		iloczyn(A, r, temp);
		alfa=iloczynsk(r,r)/iloczynsk(r,temp);
		for (int i=0;i<n;i++)
		{
			temp[i]=r[i]*alfa;
		}
		sumaw(temp, x);
		printf("%d %.6f %g %.6f \n", k, sqrt(iloczynsk(r, r)), alfa, sqrt(iloczynsk(x, x)));
		fprintf(p1, "%d \t %f \t %f \t %f\n", k, sqrt(iloczynsk(r, r)), alfa, sqrt(iloczynsk(x, x)));
	}
	while (sqrt(iloczynsk(r, r)) > pow(10, -6) && k < 500);

	k=0;
	
	for (int i=0;i<n;i++)
	{
		x[i]=1.0;
	}
	FILE *p2 = fopen("zad b.txt", "w");
	do
	{
		k++;
		iloczyn(A, x, temp);
		for(int i=0;i<n;i++)
		{
			r[i]=b[i]-temp[i];
		}
		iloczyn(A, r, temp);
		alfa=iloczynsk(r,r)/iloczynsk(r,temp);
		for (int i=0;i<n;i++)
		{
			temp[i]=r[i]*alfa;
		}
		sumaw(temp, x);
		printf("%d %.6f %g %.6f \n", k, sqrt(iloczynsk(r, r)), alfa, sqrt(iloczynsk(x, x)));
		fprintf(p2, "%d \t %f \t %f \t %f\n", k, sqrt(iloczynsk(r, r)), alfa, sqrt(iloczynsk(x, x)));
	}
	while (sqrt(iloczynsk(r, r)) > pow(10, -6) && k < 500);
	fclose(p1);
	fclose(p2);
	return 0;
}