/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - Rockefeller University                         *
 *                                                                          *
 *                       October 1997                                       *
 *                                                                          *
 * loki_identity.c:                                                         *
 *                                                                          *
 * Calculate the 9 condensed identity coefficients for relative pairs       *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_ibd.h"

#define SWAP (a,b) {swap_temp=(a);(a)=(b);(b)=swap_temp;}

static double coeffmat[9][9]={
	 {0,0,0,.25,-.25,-.25,.25,0,0},
	 {1,-1,-1,-.25,.25,.25,-.25,1,0},
	 {0,0,0,-1,1,.5,-.5,0,0},
	 {-2,2,1,1,-1,-.5,.5,-1,0},
	 {0,0,0,-1,.5,1,-.5,0,0},
	 {-2,1,2,1,-.5,-1,.5,-1,0},
	 {0,0,0,0,0,0,-.5,0,.5},
	 {0,0,0,4,-2,-2,2,0,-1},
	 {4,-2,-2,-4,2,2,-1.5,1,.5}};

static const struct Id_Record *id_array;

static double phi2(int a, int b)
{
	double x;
	int swap_temp;
	
	if(!(a&&b)) x=0.0;
	else if (a==b) x=0.5+0.5*phi2(id_array[a-1].sire,id_array[a-1].dam);
	else {
		if(a<b) {
			swap_temp=a;
			a=b;
			b=swap_temp;
		}
		x=.5*(phi2(id_array[a-1].sire,b)+phi2(id_array[a-1].dam,b));
	}
	return x;
}
static int intcomp(const void *i,const void *j)
{
	int x,y;
	
	x=*(int *)i;
	y=*(int *)j;
	if(x>y) return -1;
	if(x<y) return 1;
	return 0;
}

static double phi3(int a, int b, int c)
{
	double x;
	int p[3],ids,idd;
	
	if(!(a&&b&&c)) x=0.0;
	else {
		p[0]=a;
		p[1]=b;
		p[2]=c;
		gnu_qsort((void *)p,3,sizeof(int),intcomp);
		ids=id_array[p[0]-1].sire;
		idd=id_array[p[0]-1].dam;
		if(p[0]==p[1]) {
			if(p[0]==p[2]) x=.25*(1.0+3.0*phi2(ids,idd));
			else x=.5*(phi2(p[0],p[2])+phi3(ids,idd,p[2]));
		} else x=.5*(phi3(ids,p[1],p[2])+phi3(idd,p[1],p[2]));
	}
	return x;
}

static double phi4(int a, int b, int c, int d) 
{
	double x;
	int p[4],ids,idd;
	
	if(!(a&&b&&c&&d)) x=0.0;
	else {
		p[0]=a;
		p[1]=b;
		p[2]=c;
		p[3]=d;
		gnu_qsort((void *)p,4,sizeof(int),intcomp);
		ids=id_array[p[0]-1].sire;
		idd=id_array[p[0]-1].dam;
		if(p[0]==p[1]) {
			if(p[0]==p[2]) {
				if(p[0]==p[3]) x=.125*(1.0+7.0*phi2(ids,idd));
				else x=.25*(phi2(p[0],p[3])+3*phi3(ids,idd,p[3]));
			} else x=.5*(phi3(p[0],p[2],p[3])+phi4(ids,idd,p[2],p[3]));
		} else x=.5*(phi4(ids,p[1],p[2],p[3])+phi4(idd,p[1],p[2],p[3]));
	}
	return x;
}

static double phi22(int a, int b, int c, int d) 
{
	double x;
	int p,ids,idd;
	
	if(!(a&&b&&c&&d)) x=0.0;
	else if(a==b && a==c && a==d) x=.25*(1.0+3.0*phi2(id_array[a-1].sire,id_array[a-1].dam));
	else {
		if(a<b) {p=a;a=b;b=p;}
		if(c<d) {p=c;c=d;d=p;}
		if(a==c && b<d) {p=b;b=d;d=p;}
		if(a==b && a==c) x=.5*(phi2(a,d)+phi3(id_array[a-1].sire,id_array[a-1].dam,d));
		else {
			if(a<c) {p=a;a=c;c=p;p=b;b=d;d=p;}
			ids=id_array[a-1].sire;
			idd=id_array[a-1].dam;
			if(a==b) x=.5*(phi2(c,d)+phi22(ids,idd,c,d));
			else if(a==c) x=.25*(2.0*phi3(a,b,d)+phi22(ids,b,idd,d)+phi22(idd,b,ids,d));
			else x=.5*(phi22(ids,b,c,d)+phi22(idd,b,c,d));
		}
	}
	return x;
}

void loki_identity(double *k,int a,int b,const struct Id_Record *idr)
{
	int i,j;
	double p2,y[9];
 
	id_array=idr;
	for(i=0;i<9;i++) k[i]=0.0;
	p2=phi2(a,b);
	if(p2>0.0) {
		y[0]=1.0;
		y[1]=2.0*phi2(a,a);
		y[2]=2.0*phi2(b,b);
		y[3]=4.0*p2;
		y[4]=8.0*phi3(a,a,b);
		y[5]=8.0*phi3(a,b,b);
		y[6]=16.0*phi4(a,a,b,b);
		y[7]=4.0*phi22(a,a,b,b);
		y[8]=16.0*phi22(a,b,a,b);
		for(i=0;i<9;i++) for(j=0;j<9;j++) k[i]+=coeffmat[i][j]*y[j];
	} else k[8]=1.0;
}
