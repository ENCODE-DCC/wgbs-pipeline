/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * sample_rand.c:                                                           *
 *                                                                          *
 * Sampling routine for uncorrelated random variance components             *
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
#include "mat_utils.h"
#include "sparse.h"
#include "loki.h"
#include "loki_utils.h"
#include "sample_rand.h"

static double *work;
static struct loki *lk;

static void free_rand(void)
{
	if(lk->models->c_var) {
		if(lk->models->c_var[0]) free(lk->models->c_var[0]);
		free(lk->models->c_var);
	}
	if(lk->models->rand_list) free(lk->models->rand_list);
	if(lk->models->rand_flag) free(lk->models->rand_flag);
	if(work) free(work);
	lk->models->rand_flag=0;
	lk->models->rand_list=0;
	lk->models->c_var=0;
	work=0;
}

void init_rand(struct loki *loki)
{
	int k,i,j,type,mod;
	struct Variable *var;
	double *tmp;
	struct Model *models;
	
	if(!loki->models->n_models) return;
	lk=loki;
	if(atexit(free_rand)) message(WARN_MSG,"Unable to register exit function free_rand()\n");
	for(loki->models->n_random=i=0;i<loki->data->n_id_records;i++) if(loki->data->id_variable[i].type&ST_RANDOM) loki->models->n_random++;
	for(i=0;i<loki->data->n_nonid_records;i++) if(loki->data->nonid_variable[i].type&ST_RANDOM) loki->models->n_random++;
	if(!(loki->models->rand_flag=malloc(sizeof(unsigned long)*loki->models->n_models))) ABT_FUNC(MMsg);
	if(loki->models->n_random) {
		models=loki->models->models;
		if(!(loki->models->rand_list=malloc(sizeof(struct Variable *)*loki->models->n_random))) ABT_FUNC(MMsg);
		for(j=i=0;i<loki->data->n_id_records;i++) if(loki->data->id_variable[i].type&ST_RANDOM) loki->models->rand_list[j++]=loki->data->id_variable+i;
		for(i=0;i<loki->data->n_nonid_records;i++) if(loki->data->nonid_variable[i].type&ST_RANDOM) loki->models->rand_list[j++]=loki->data->nonid_variable+i;
		for(mod=0;mod<loki->models->n_models;mod++) {
			loki->models->rand_flag[mod]=0;
			for(k=0;k<models[mod].n_terms;k++) {
				type=models[mod].term[k].vars[0].type;
				if(type&ST_RANDOM) {
					i=models[mod].term[k].vars[0].var_index;
					var=(type&ST_CONSTANT)?loki->data->id_variable+i:loki->data->nonid_variable+i;
					for(i=0;i<loki->models->n_random;i++) if(loki->models->rand_list[i]==var) break;
					if(i==loki->models->n_random) ABT_FUNC("Internal error - random variable not found\n");
					loki->models->rand_flag[mod]|=1L<<i;
				}
			}
		}
		k=loki->models->n_models*(loki->models->n_models+1)/2;
		if(!(loki->models->c_var=malloc(sizeof(double *)*loki->models->n_random))) ABT_FUNC(MMsg);
		if(!(tmp=malloc(sizeof(double)*k*loki->models->n_random))) ABT_FUNC(MMsg);
		for(i=0;i<k*loki->models->n_random;i++) tmp[i]=0.0;
		for(i=0;i<loki->models->n_random;i++) {
			loki->models->c_var[i]=tmp;
			for(j=0;j<loki->models->n_models;j++) if(loki->models->rand_flag[j]&(1L<<i)) BB(tmp,j,j)=BB(loki->models->residual_var,j,j);
			tmp+=k;
		}
	}
	for(j=i=0;i<loki->models->n_models;i++) if(loki->models->models[i].polygenic_flag) j++;
	if(j) if(!(work=malloc(sizeof(double)*loki->pedigree->ped_size))) ABT_FUNC(MMsg);
}

void sample_rand(struct loki *loki)
{
	int i,k,r,type,df;
	double *eff,y,ss,n;
	struct Model *models;
	
	models=loki->models->models;
	for(r=k=0;k<models[0].n_terms && r<loki->models->n_random;k++) {
		type=models[0].term[k].vars[0].type;
		if(type&ST_RANDOM) {
			if(!(type&ST_FACTOR)) ABT_FUNC("Can't fit random effect with continuous covariates\n");
			df=models[0].term[k].df;
			eff=models[0].term[k].eff;
			ss=0.0;
			for(i=0;i<df;i++) {
				y=eff[i];
				ss+=y*y;
			}
			n=(double)df+RES_PRIOR_VC0;
			ss+=RES_PRIOR_VC0*RES_PRIOR_SC0;
			loki->models->c_var[r++][0]=ss/(sgamma(n*.5)*2.0);
		}
	}
}

void sample_additive_var(struct loki *loki)
{
	int i,id,j,k,comp;
	double y,n,z,ss,cs;
	struct SparseMatRec *AI;
	struct Id_Record *id_array;
	
	ss=0.0;
	id_array=loki->pedigree->id_array;
	for(id=comp=0;comp<loki->pedigree->n_comp;comp++) {
		AI=loki->models->AIMatrix[comp];
		cs=loki->pedigree->comp_size[comp];
		for(i=0;i<cs;i++) {
			y=id_array[i+id].bv[0];
			work[i]=y*AI[i].val;
			for(j=AI[i].x;j<AI[i+1].x;j++) {
				k=AI[j].x;
				z=AI[j].val;
				work[k]+=y*z;
				work[i]+=id_array[k+id].bv[0]*z;
			}
		}
		for(i=0;i<cs;i++) ss+=id_array[i+id].bv[0]*work[i];
		id+=i;
	}
	n=(double)loki->pedigree->ped_size+RES_PRIOR_VA0;
	ss+=RES_PRIOR_VA0*RES_PRIOR_SA0;
	loki->models->additive_var[0]=ss/(sgamma(n*.5)*2.0);
}
