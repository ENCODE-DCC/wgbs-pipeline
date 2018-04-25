#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>

int main(void) 
{
	int i,j,k,k1,n=23,s,s1;
	double *p,*p1,*p2,*pp,theta,x,z,z1,z2;
	int fd,fd1;
	
	s=1<<n;
	fd=open("junk",O_RDONLY,0);
	fd1=open("junk1",O_RDWR,0600);
	if(fd<0 || fd1<0) exit(0);
	theta=0.1;
	x=(1.0-theta)/theta;
	pp=mmap(0,sizeof(double)*s,PROT_READ,0,fd,0);
	p=mmap(0,sizeof(double)*s,PROT_READ|PROT_WRITE,0,fd1,0);
	if(pp==MAP_FAILED) {
		printf("couldn't mmap input file\n");
		perror(0);
		exit(-1);
	}
	if(p==MAP_FAILED) {
		printf("couldn't mmap output file\n");
		perror(0);
		exit(-1);
	}
	printf("Starting convolution\n");
	k1=s1=s>>1;
	p1=p;
	p2=p+1;
	for(k=i=0;i<k1;i++,p1+=2,p2+=2) {
		z1=pp[k++];
		z2=pp[k++];
		*p2=x*z2+z1;
		*p1=x*z1+z2;
	}
	munmap(pp,sizeof(double)*s);
	close(fd);
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
	} 
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
	for(i=0;i<10;i++) printf("%d %10g %g\n",i,p[i],z); 
	munmap(p,sizeof(double)*s);
	close(fd1);
}
