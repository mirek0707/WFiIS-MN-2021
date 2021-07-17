#include<stdio.h>
int main()
{
	int n=500;
	double Xa=0.5, Xb=2;
	double h=2*Xb/(n-1);
	double u[n+1], l[n+1];
	double a[n+1], c[n+1], d[n+1];
	double y[n+1], x[n+1], b[n+1], v[n+1], vt[n+1];

	for (int i=1;i<=n;i++)
		x[i]=-Xb+h*(i-1);

	for (int i=1;i<=n;i++)
		d[i]=-2/(h*h);

	for (int i=1;i<=n;i++)
		a[i]=1/(h*h);

	for (int i=1;i<=n;i++)
		c[i]=1/(h*h);

	for(int i=1; i<=n; i++)
	{
    	if(x[i]>=-Xb && x[i]<-Xa || x[i]==0 || x[i]>Xa   && x[i]<=Xb)
        	b[i]=0;
        if(x[i]>=-Xa && x[i]<0) 
        	b[i]=-1;
     	if(x[i]>0 && x[i]<=Xa) 
     		b[i]=1;
    }
    d[1]=1;
    c[1]=0;
    b[1]=0;

    d[n]=1;
    a[n]=0;
    b[n]=0;

    u[1]=d[1];
    for (int i=2;i<=n;i++)
    {
    	l[i]=a[i]/u[i-1];
    	u[i]=d[i]-l[i]*c[i-1];
    }

    y[1]=b[1];
    for(int i=2; i<=n; i++)
        y[i]=b[i]-l[i]*y[i-1];

    v[n]=y[n]/u[n];
    for(int i=n-1;i>0;i--)
    	v[i]=(y[i]-c[i]*v[i+1])/u[i];

	for (int i=0;i<=n;i++)
	{
		if (x[i]>=-Xb && x[i]<=-Xa)
			vt[i]=x[i]/16+1.0/8.0;
		if (x[i]>=-Xa && x[i]<=0)
			vt[i]=-(x[i]*x[i])/2-(7*x[i])/16;
		if (x[i]>=0 && x[i]<=Xa)
			vt[i]=(x[i]*x[i])/2-(7*x[i])/16;
		if (x[i]>=Xa && x[i]<=Xb)
			vt[i]=x[i]/16-1.0/8.0;
	}

	FILE *p;
	p=fopen("N500.txt", "w");
	for(int i=1;i<=n;i++)
	{
		fprintf(p, "%10.4g %10.4g %10.4g \n", x[i], v[i], vt[i]);
	}
	fclose(p);
	return 0;

}