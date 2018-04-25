#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void) 
{
	int i,j,k,k1,n=24,s,s1;
	double *p,*p1,*p2,theta,x,z,z1,z2,a,b,c,d,x2;
	
	s=1<<n;
	theta=0.1;
	x=(1.0-theta)/theta;
	x2=x*x;
	if(!(p=malloc(sizeof(double)*s))) {
		printf("Can't allocate memory\n");
		exit(-1);
	}
	p1=p;
	for(i=0;i<s;i++) *p1++=1.0;
/*	p[0]=.1;
	p[1]=.2;
	p[2]=.3;
	p[3]=.4;
	p[4]=.1;
	p[5]=.2;
	p[6]=.3;
	p[7]=.4;
	p[8]=.1;
	p[9]=.2;
	p[10]=.3;
	p[11]=.4;
	p[12]=.1;
	p[13]=.2;
	p[14]=.3;
	p[15]=.4; */
/*	for(z=0.0,i=0;i<s;i++) z+=p[i];
	for(i=0;i<s;i++) p[i]/=z; */
	printf("Starting convolution\n");
	k1=s>>2;
	s1=s>>1;
	p1=p;
	for(i=0;i<k1;i++) {
		a=p1[0];
		b=p1[1];
		c=p1[2];
		d=p1[3];
		z=x*(b+c);
		z1=x*(a+d);
		*p1++=a*x2+z+d;
		*p1++=b*x2+z1+c;
		*p1++=c*x2+z1+b;
		*p1++=d*x2+z+a;
	}
	/*
	p1=p1;
	p2=p+1;
	for(i=0;i<k1;i++,p1+=2,p2+=2) {
		z1=*p1;
		z2=*p2;
		*p2=x*z2+z1;
		*p1=x*z1+z2;
	}
	k1>>=1;
	p1=p;
	p2=p+2;
	for(i=0;i<k1;i++,p1+=3,p2+=3) {
		z1=*p1;
		z2=*p2;
		*p2++=x*z2+z1;
		*p1++=x*z1+z2;
		z1=*p1;
		z2=*p2;
		*p2=x*z2+z1;
		*p1=x*z1+z2;
	} */ 
	k1>>=1; 
	p1=p;
	p2=p+4;
	for(i=0;i<k1;i++,p1+=5,p2+=5) {
		z1=*p1;
		z2=*p2;
		*p2++=x*z2+z1;
		*p1++=x*z1+z2;
		z1=*p1;
		z2=*p2;
		*p2++=x*z2+z1;
		*p1++=x*z1+z2;
		z1=*p1;
		z2=*p2;
		*p2++=x*z2+z1;
		*p1++=x*z1+z2;
		z1=*p1;
		z2=*p2;
		*p2=x*z2+z1;
		*p1=x*z1+z2;
	} 
	k1>>=1;
	j=8; 
	while(j<s1) {
		p1=p;
		p2=p1+j;
		for(i=0;i<k1;i++) {
		  for(k=0;k<j;k++) {
		    z1=*p1;
		    z2=*p2;
		    *p2++=x*z2+z1;
		    *p1++=x*z1+z2;
		  }
		  p1+=j;
		  p2+=j;
		}
		j<<=1;
		k1>>=1;
	}
	p1=p;
	p2=p1+j;
	z=exp(log(theta)*(double)n);
	for(i=0;i<k1;i++) {
	  for(k=0;k<j;k++) {
	    z1=*p1*z;
	    z2=*p2*z;
	    *p2++=x*z2+z1;
	    *p1++=x*z1+z2;
	  }
	  p1+=j;
	  p2+=j;
	}
	for(i=0;i<16;i++) printf("%d %10g %g\n",i,p[i],z); 
}
