#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <stdio.h>
#include <math.h>
double delta(int i, int j)
{
	double d;
	if (i==j)
	{
		d=1;
	}
	else
	{
		d=0;
	}
	return d;
}
int main ()
{
	FILE *a0, *a100, *eigen;
    a0 = fopen("a0.txt", "w");
    a100 = fopen("a100.txt", "w");
    eigen = fopen("eigen.txt", "w");

	double L=10, n=200, N=1;
	int alfa=0;
	double deltaX=(double)L/(n+1);

	gsl_vector *x=gsl_vector_calloc(n);
	gsl_vector *ro=gsl_vector_calloc(n);
	gsl_matrix *A=gsl_matrix_calloc(n, n);
	gsl_matrix *B=gsl_matrix_calloc(n, n);
	gsl_vector *eval = gsl_vector_calloc(n);
    gsl_matrix *evec = gsl_matrix_calloc(n, n);
    gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);

	for (alfa=0;alfa<101;alfa+=2)
	{
		fprintf(eigen, "%d ", alfa);
		
		for(int i=0;i<n;i++)//x
		{
			gsl_vector_set(x, i, (double)(-L/2+deltaX*(i+1)));
		}

		
		for(int i=0;i<n;i++)//ro
		{
			gsl_vector_set(ro, i, (double)(1+4*alfa*pow(gsl_vector_get(x,i), 2)));
		}

		
	    for (int i=0;i<n;i++)//A,B
	    {
	    	for(int j=0;j<n;j++)
	    	{
	    		gsl_matrix_set(A,i,j,(double)((-delta(i,j+1)+2*delta(i,j)-delta(i,j-1))/pow(deltaX,2)));
	    		gsl_matrix_set(B,i,j, (gsl_vector_get(ro,i)/N*delta(i,j)));
	    	}
	    }

	    gsl_eigen_gensymmv(A, B, eval, evec, w);
        gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

        for (int i=0;i<6;i++)
        {
            fprintf(eigen, "%f ", sqrt(gsl_vector_get(eval, i)));
        }
        fprintf(eigen, "\n");

        if (alfa==0)
        {
        	for (int j=0;j<n;j++)
        	{
        		fprintf(a0, "%f ", gsl_vector_get(x,j));
        		for (int i=0;i<6;i++)
        		{
        			fprintf(a0, "%f ", gsl_matrix_get(evec, j, i));
        		}
        		fprintf(a0, "\n");
        	}
        }
        if (alfa==100)
        {
        	for (int j=0;j<n;j++)
        	{
        		fprintf(a100, "%f ", gsl_vector_get(x,j));
        		for (int i=0;i<6;i++)
        		{
        			fprintf(a100, "%f ", gsl_matrix_get(evec, j, i));
        		}
        		fprintf(a100, "\n");
        	}
        }

	}

	gsl_vector_free (x);
	gsl_vector_free(ro);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_gensymmv_free (w);

	fclose(a0);
    fclose(a100);
    fclose(eigen);

    return 0;
		
}