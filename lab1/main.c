#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
int main ()
{
	int n=4, i, j, k, signum;
	double value;
	gsl_matrix *a=gsl_matrix_calloc(n,n);
	gsl_matrix *acopy=gsl_matrix_calloc(n,n);
	gsl_matrix *ainv=gsl_matrix_calloc(n,n);
	gsl_matrix *c=gsl_matrix_calloc(n,n);

	gsl_permutation *p=gsl_permutation_calloc(n);

	gsl_vector *b=gsl_vector_calloc(n);
	gsl_vector *x=gsl_vector_calloc(n);

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			value=1/(i+j+2.0);
			gsl_matrix_set(a,i,j,value);
			gsl_matrix_set(acopy,i,j,value);
		}
	}
	printf("\n\n macierz A\n");
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("  %10.5g ", gsl_matrix_get(a,i,j));
		}
		printf("\n");
	}

	gsl_linalg_LU_decomp(a,p,&signum);
	printf("\n\n macierz LU\n");
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("  %10.5g ", gsl_matrix_get(a,i,j));
		}
		printf("\n");
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			gsl_vector_set(b,j,0.0);
		}
		gsl_vector_set(b,i,1.0);
		gsl_linalg_LU_solve(a,p,b,x);
		for (j=0;j<n;j++)
		{
			gsl_matrix_set(ainv, j,i, gsl_vector_get(x,j));
		}
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			double suma=0.0;
			for(k=0;k<n;k++)
			{
				double aik=gsl_matrix_get(acopy,i,k);
				double bkj=gsl_matrix_get(ainv,k,j);
				suma=suma+aik*bkj;
			}
			gsl_matrix_set(c,i,j,suma);
		}
	}


	printf("\n\n macierz jednostkowa C=A*A^(-1)\n");
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("  %10.5g ", gsl_matrix_get(c,i,j));
		}
		printf("\n");
	}


	double amax=0.0, ainv_max=0.0, z, cond;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			z=fabs(gsl_matrix_get(acopy,i,j));
			if(z>amax) 
				amax=z;
			z=fabs(gsl_matrix_get(ainv,i,j));
			if(z>ainv_max) 
				ainv_max=z;
		}
	}
	cond=amax*ainv_max;

	printf("\n\n wskaznik uwarunkowania = %f \n", cond);

	double detA = 1.0;
	for (i=0;i<n;i++)
	{
		detA=detA*gsl_matrix_get(a,i,i);
	}
	printf("\n\n wyznacznik = %e \n", detA);
	return 0;
}