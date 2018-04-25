#include <config.h>

#include <stdio.h>
#include <math.h>

#include "ranlib.h"

#define A0 2.50662823884
#define A1 -18.61500062529
#define A2  41.39119773534
#define A3 -25.44106049637
#define B1  -8.47351093090
#define B2  23.08336743743
#define B3 -21.06224101826
#define B4   3.13082909833
#define C0  -2.78718931138
#define C1  -2.29796479134
#define C2   4.85014127135
#define C3   2.32121276858
#define D1   3.54388924762
#define D2   1.63706781897
#define SPLIT 0.42

double ppnd(double p,int *fault)
{
	double q,r;
	
	*fault=0;
	q=p-0.5;
	if(fabs(q)<SPLIT)
	{
		r=q*q;
		return q*(((A3*r+A2)*r+A1)*r+A0)/((((B4*r+B3)*r+B2)*r+B1)*r+1.0);
	}
	r=(q>0.0)?1.0-p:p;
	if(r<=0.0)
	{
		*fault=1;
		return 0.0;
	}
	r=sqrt(-log(r));
	p=(((C3*r+C2)*r+C1)*r+C0)/((D2*r+D1)*r+1.0);
	return (q<0.0)?-p:p;
}

#undef A0 
#undef A1
#undef A2  
#undef A3
#undef B1
#undef B2  
#undef B3
#undef B4   
#undef C0
#undef C1
#undef C2   
#undef C3   
#undef D1   
#undef D2   
#undef SPLIT 

double trunc_normal(const double a,const double b,const int flag,int *err)
{
	double u,t,t1,p;
	
	u=(double)ranf();
	t=flag==1?0.0:.5*(1.0+erf(a/sqrt(2.0)));
	t1=flag==2?1.0:.5*(1.0+erf(b/sqrt(2.0)));
	p=ppnd(t+u*(t1-t),err);
	return p;
}

void dirichlet(const int n,const double *ct,double *fq)
{
	int i;
	double z=1.0,z1=0.0;
	
	for(i=0;i<n;i++) z1+=ct[i];
	for(i=0;i<n-1;i++) {
		z1-=ct[i];
		fq[i]=z*genbet(ct[i],z1);
		z-=fq[i];
	}
	fq[i]=z;
}

