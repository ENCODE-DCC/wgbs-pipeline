/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * calc_var_locus.c:                                                        *
 *                                                                          *
 * Calculate variance contributed by a trait locus                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_utils.h"
#include "mat_utils.h"
#include "calc_var_locus.h"

static void free_calc_var_locus(void) 
{
	calc_var_locus(0,0);
}

void calc_var_locus(struct Locus *loc,struct loki *loki)
{
	int i,j,k,mod,mod1,n_models;
	double z,z1;
	int *gt;
	static double *x,*mu,*n;
	unsigned long mflag=0,a,b;
	struct Id_Record *id_array;
	struct Marker *marker;
	
	if(!loc) {
		if(x) {
			free(x);
			x=0;
		}
		return;
	}
	n_models=loki->models->n_models;
	id_array=loki->pedigree->id_array;
	marker=loki->markers->marker;
	k=n_models*(n_models+1)/2;
	if(!x) {
		if(!(x=malloc(sizeof(double)*(n_models+2*k)))) ABT_FUNC(MMsg);
		n=x+k;
		mu=n+k;
		if(atexit(free_calc_var_locus)) message(WARN_MSG,"Unable to register exit function free_calc_var_locus()\n");
	}
	gt=loc->gt;
	mflag=loc->model_flag;
	for(i=0;i<k;i++) x[i]=n[i]=0.0;
	for(i=0;i<n_models;i++) mu[i]=0.0;
	for(i=0;i<loki->pedigree->ped_size;i++) {
		k=gt[i]-1;
		for(a=1,mod=0;mod<n_models;mod++,a<<=1) if((mflag&a)&&id_array[i].res[mod]) {
			if(k) mu[mod]+=loc->eff[mod][k-1];
			n[mod]++;
		}
	}
	for(mod=0;mod<n_models;mod++) {
		mu[mod]=n[mod]>0.0?mu[mod]/n[mod]:0.0;
		n[mod]=0.0;
	}
	for(i=0;i<loki->pedigree->ped_size;i++) {
		k=gt[i]-1;
		for(a=1,mod=0;mod<n_models;mod++,a<<=1) if((mflag&a)&&id_array[i].res[mod]) {
			z=k?loc->eff[mod][k-1]:0.0;
			z-=mu[mod];
			for(b=1,mod1=0;mod1<=mod;mod1++,b<<=1) if((mflag&b)&&id_array[i].res[mod1]) {
				z1=k?loc->eff[mod1][k-1]:0.0;
				z1-=mu[mod];
				j=IDX(mod,mod1);
				n[j]++;
				x[j]+=z*z1;
			}
		}
	}
	for(j=mod=0;mod<n_models;mod++) {
		for(mod1=0;mod1<=mod;mod1++) {
			z=n[j];
			if(z>0.0) x[j]/=z;
			if(x[j]<1.0e-16) x[j]=0.0;
			j++;
		}
	}
	k=n_models*(n_models+1)/2;
	for(i=0;i<k;i++) loc->variance[i]=x[i];
}

