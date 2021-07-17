#include <stdio.h>
#include <math.h>
#include <time.h>
#include </usr/include/gsl/gsl_math.h>
#include </usr/include/gsl/gsl_linalg.h>

const double x0 = 2.0;
const double sigma = 4.0;

double g(double x)
{
	double a0=-pow(x0,2)/(2*pow(sigma,2));
	double a1=x0/(pow(sigma,2));
	double a2=-1/(2*pow(sigma,2));
	return exp(a0+a1*x+a2*pow(x,2));
}

double g2(double x, double alfa)
{
	double U=rand()/(RAND_MAX+1.0);
	return g(x)*(1.0+alfa*(U - 0.5));
}

double G(double x, gsl_vector *vector)
{
    double b0=gsl_vector_get(vector,0);
    double b1=gsl_vector_get(vector,1);
    double b2=gsl_vector_get(vector,2);
    double b3=gsl_vector_get(vector,3);
    return exp(b0 + b1*x + b2*pow(x,2) + b3*pow(x,3));
}

double d(double N)
{
    return ((3*sigma+x0)-(-3*sigma+x0))/(N-1);
}

void zapis(int N, gsl_vector *w, FILE *p)
{
    for(int i = 0; i <= N; i++)
    {
        double x=(-3*sigma+x0)+0.1*i;
        double lG=G(x, w);
        fprintf(p, "%g %g \n",x, lG);
    }
    fprintf(p, "\n \n");
}

void zapisWezly(int N, gsl_vector *wezly, gsl_vector *wartosc, FILE *p)
{
    for(int i = 0; i < N; i++)
    {
        fprintf(p, "%g %g \n",gsl_vector_get(wezly, i),gsl_vector_get(wartosc, i));
    }
    fprintf(p, "\n \n");
}

int main()
{
	srand(time(NULL));
	const int N1=11;
    const int m=4;
    const int N2=101;

    FILE *p1=fopen("wezly.txt","w");
    FILE *p2=fopen("G.txt","w");

    gsl_vector *wezly=gsl_vector_calloc(N1);
    for(int i=0;i<N1;i++)
    {
        gsl_vector_set(wezly,i,-3*sigma+x0+d(N1)*i);
    }

    gsl_vector *wezly2=gsl_vector_calloc(N2);
    for(int i=0;i<N2;i++)
    {
        gsl_vector_set(wezly2,i,-3*sigma+x0+d(N2)*i); 
    }


    gsl_vector *wartosc_g=gsl_vector_calloc(N1);
    gsl_vector *wartosc_g2=gsl_vector_calloc(N1);
    for(int i=0;i<N1;i++)
    {
        gsl_vector_set(wartosc_g,i,g(gsl_vector_get(wezly,i)));
        gsl_vector_set(wartosc_g2,i,g2(gsl_vector_get(wezly,i),0.5));
    }

    gsl_vector *wartosc_g2_1=gsl_vector_calloc(N2);
    for(int i=0;i<N2;i++)
    {
        gsl_vector_set(wartosc_g2_1,i,g2(gsl_vector_get(wezly2,i),0.5));
    }

    zapisWezly(N1,wezly,wartosc_g,p1);
    zapisWezly(N1,wezly,wartosc_g2,p1);
    zapisWezly(N2,wezly2,wartosc_g2_1,p1);

    gsl_vector *wartosc_f=gsl_vector_calloc(N1);
    gsl_vector *wartosc_f2=gsl_vector_calloc(N1);

    for(int i=0;i<N1;i++)
    {
        gsl_vector_set(wartosc_f,i,log(gsl_vector_get(wartosc_g,i)));
        gsl_vector_set(wartosc_f2,i,log(gsl_vector_get(wartosc_g2,i)));
    }

    gsl_vector *wartosc_f2_1=gsl_vector_calloc(N2);
    for(int i=0;i<N2;i++)
    {
        gsl_vector_set(wartosc_f2_1,i,log(gsl_vector_get(wartosc_g2_1, i)));   
    }

    gsl_vector *r=gsl_vector_calloc(m);
    gsl_vector *r2=gsl_vector_calloc(m);
    gsl_vector *r2_1=gsl_vector_calloc(m);

    gsl_matrix *mG=gsl_matrix_calloc(m, m);
    gsl_matrix *mG2=gsl_matrix_calloc(m, m);
    gsl_matrix *mG2_1=gsl_matrix_calloc(m, m);


    

    for(int i=0;i<m;i++)
    {
        double temp=0.0;
        double temp2=0.0; 
        for(int j=0;j<N1;j++)
        {
            temp+=gsl_vector_get(wartosc_f,j)*pow(gsl_vector_get(wezly,j),i);
            temp2+=gsl_vector_get(wartosc_f2, j)*pow(gsl_vector_get(wezly,j),i);
        }
        gsl_vector_set(r,i,temp);
        gsl_vector_set(r2,i,temp2);
    }

    for(int i=0;i<m;i++)
    {
        double temp=0.0;
        for(int j=0;j<N2;j++)
        {
            temp+=gsl_vector_get(wartosc_f2_1,j)*pow(gsl_vector_get(wezly2,j),i);
        }
        gsl_vector_set(r2_1,i,temp);
    }

     for(int i=0;i<m;i++)
     {
		for(int k=0;k<m;k++) 
		{
			double temp=0;
			for(int j=0;j<N1;j++) 
			{
				temp+=pow(gsl_vector_get(wezly,j),i+k);
			}
		gsl_matrix_set(mG,i,k,temp);
        gsl_matrix_set(mG2,i,k,temp);
		}
	}

	for(int i=0;i<m;i++)
	{
		for(int k=0;k<m;k++)
		{
			double temp=0;
			for(int j=0;j<N2;j++)
			{
				temp+=pow(gsl_vector_get(wezly2,j),i+k);
			}
		gsl_matrix_set(mG2_1,i,k,temp);

		}
	}

	gsl_linalg_HH_svx(mG,r);
    gsl_linalg_HH_svx(mG2,r2);
    gsl_linalg_HH_svx(mG2_1,r2_1);

    double a0 = -pow(x0, 2) / (2 * pow(sigma, 2));
    double a1 = x0 / pow(sigma, 2);
    double a2 = -1 / (2 * pow(sigma, 2));

    printf("Analityczne(g1 oraz n = 11) \n");
    printf("%g \n", a0);
    printf("%g \n", a1);
    printf("%g \n", a2);
    printf("\n");

    printf("Numeryczne(g1 oraz n = 11) \n");
    for(int i = 0; i < m; i++){
        printf("%g \n",gsl_vector_get(r, i));
    }
    printf("\n");

    printf("Numeryczne(g2 oraz n = 11) \n");
    for(int i = 0; i < m; i++){
        printf("%g \n",gsl_vector_get(r2, i));
    }
    printf("\n");

    printf("Numeryczne(g2 oraz n = 101) \n");
    for(int i = 0; i < m; i++){
        printf("%g \n",gsl_vector_get(r2_1, i));
    }
    printf("\n");

    int k=((3*sigma+x0)-(-3*sigma+x0))/0.1;
    zapis(k,r,p2);
    zapis(k,r2,p2);
    zapis(k,r2_1,p2);

    gsl_vector_free(wezly);
    gsl_vector_free(wezly2);
    gsl_vector_free(wartosc_g);
    gsl_vector_free(wartosc_g2);
    gsl_vector_free(wartosc_g2_1);
    gsl_vector_free(wartosc_f);
    gsl_vector_free(wartosc_f2);
    gsl_vector_free(wartosc_f2_1);
    gsl_vector_free(r);
    gsl_vector_free(r2);
    gsl_vector_free(r2_1);
    gsl_matrix_free(mG);
    gsl_matrix_free(mG2);
    gsl_matrix_free(mG2_1);

    fclose(p1);
    fclose(p2);

    return 0;
}