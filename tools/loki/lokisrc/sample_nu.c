/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * sample_nu.c:                                                             *
 *                                                                          *
 * Sampling routines connected to t-distribution of residuals               *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <errno.h>

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "sample_nu.h"

static double nu_list[]={1,1.5,2,5,10,25,50,100,-1};
static double *nu_prob;
static int n_nu;

void init_sample_nu(struct loki *loki)
{
	int i,type,nrec;
	double n,v;
	struct Id_Record *idr;
	
	n_nu=0;
	while(nu_list[n_nu]>=0.0) n_nu++;	
	if(n_nu) {
		loki->models->res_nu=nu_list[n_nu-1];
		if(!(nu_prob=malloc(sizeof(double)*2*n_nu))) ABT_FUNC(MMsg);
	}
	n=0.0;
	type=loki->models->models[0].var.type;
	idr=loki->pedigree->id_array;
	for(i=0;i<loki->pedigree->ped_size;i++,idr++) if(idr->res) {
		if(type&ST_CONSTANT) nrec=1;
		else nrec=idr->n_rec;
		n+=(double)nrec;
	}
	for(i=0;i<n_nu;i++) {
		v=nu_list[i];
		nu_prob[i]=lgamma((v+n)*.5)-lgamma(.5*v)-log(v*M_PI*loki->models->residual_var[0])*.5*n;
	}
	loki->models->res_nu=.25;
}

void free_sample_nu(void)
{
	if(nu_prob) free(nu_prob);
	nu_prob=0;
}

static void sample_weights(struct loki *loki)
{
	int i,j,type,nrec;
	double y,n,s,kk,res;
	struct Id_Record *idr;
	
	type=loki->models->models[0].var.type;
	if(loki->models->res_nu<1.0e-6) {
		idr=loki->pedigree->id_array;
		for(i=0;i<loki->pedigree->ped_size;i++,idr++) if(idr->res[0]) {
			if(type&ST_CONSTANT) nrec=1;
			else nrec=idr->n_rec;
			for(j=0;j<nrec;j++) idr->vv[0][j]=1.0;
		}
	} else {
		res=loki->models->residual_var[0];
		kk=res/loki->models->res_nu;
		n=(1.0/loki->models->res_nu)+1.0;
		idr=loki->pedigree->id_array;
		for(i=0;i<loki->pedigree->ped_size;i++,idr++) if(idr->res[0]) {
			if(type&ST_CONSTANT) nrec=1;
			else nrec=idr->n_rec;
			for(j=0;j<nrec;j++) {
				y=idr->res[0][j];
				s=y*y;
				y=(s+kk)/(sgamma(n*0.5)*2.0);
				idr->vv[0][j]=y/res;
			}
		}
	}
}

void sample_nu(struct loki *loki)
{
	loki->models->res_nu=.25;
	sample_weights(loki);
	return;
}

