#include<stdio.h>
#include<math.h>
#define h 0.1
#define it 1000

double f(double x, double y)
{
	return 5.0*pow((pow(x,2)-y),2)/2.0 +pow(1-x,2);
}
double dfdx(double x, double y, double deltax)
{
	return (f(x+deltax, y)-f(x-deltax, y))/(2.0*deltax);
}
double dfdy(double x, double y, double deltay)
{
	return (f(x, y+deltay)-f(x, y-deltay))/(2.0*deltay);
}

int main()
{
	FILE* p1 = fopen("eps1.dat", "w");
    FILE* p2 = fopen("eps2.dat", "w");

    double delta = 0.0001;;
	double e1 = 0.01;
	double e2 = 0.001;
	double x0 = -0.75;
    double y0 = 1.75;

    fprintf(p1, "%f %f\n", x0, y0);
	fprintf(p2, "%f %f\n", x0, y0);

	for(int i=0; i<it; i++)
	{
		double x1 = x0-h*dfdx(x0, y0, delta);
		double y1 = y0-h*dfdy(x0, y0, delta);

		fprintf(p1, "%f %f\n", x1, y1);
		if (sqrt(pow((x1-x0), 2) + pow((y1-y0), 2)) < e1) 
		{
	        break;
	    }
	    x0=x1;
		y0=y1;
	}

	x0 = -0.75;
    y0 = 1.75;
    for(int i=0; i<it; i++)
	{
		double x1 = x0-h*dfdx(x0, y0, delta);
		double y1 = y0-h*dfdy(x0, y0, delta);

		fprintf(p2, "%f %f\n", x1, y1);
		if (sqrt(pow((x1-x0), 2) + pow((y1-y0), 2)) < e2) 
		{
	        break;
	    }
	    x0=x1;
		y0=y1;
	}

	fclose(p1);
    fclose(p2);
	return 0;
}