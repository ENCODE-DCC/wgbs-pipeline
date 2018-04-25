/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*             Simon Heath - University of Washington                       *
*                                                                          *
*                       July 1997                                          *
*                                                                          *
* alloc_loki.c:                                                            *
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
#include "lk_malloc.h"
#include "loki_utils.h"

void init_traitlocus(struct loki *loki,int i)
{
	struct Locus *loc;
	
	loc=loki->models->tlocus+i;
	loc->freq=0;
	loc->flag=0;
	loc->seg[0]=0;
	loc->lk_store=0;
	loc->variance=0;
	loc->gt=0;
	loc->index=i;
	loc->type=ST_TRAITLOCUS;
	loc->eff=0;
	loc->model_flag=0;
	loc->pruned_flag=loki->models->pruned_flag;
	loc->founder_flag=loki->models->founder_flag;
	loc->name=loki->models->qtl_name;
}

void init_trait_loci(struct loki *loki)
{
	int i;
	
	loki->params.n_tloci=(loki->models->qtl_name&&loki->models->models)?DEFAULT_N_TLOCI:0;
	if(loki->params.n_tloci>loki->params.max_tloci) loki->params.n_tloci=loki->params.max_tloci;
	if(loki->params.n_tloci) {
		/* This is just enough for starters.  If we need more then we
		* can realloc this no? */
		if(!(loki->models->tlocus=malloc(sizeof(struct Locus)*loki->params.n_tloci))) ABT_FUNC(MMsg);
		if(!(loki->models->pruned_flag=malloc(sizeof(int)*2*loki->pedigree->ped_size))) ABT_FUNC(MMsg);
		loki->models->founder_flag=loki->models->pruned_flag+loki->pedigree->ped_size;
		for(i=0;i<loki->params.n_tloci;i++) init_traitlocus(loki,i);
	} else {
		loki->models->tlocus=0;
		loki->models->pruned_flag=0;
	}
}

int get_new_traitlocus(const int n_all,struct loki *loki)
{
	int i,j,k,*ff;
	double *tp,**tpp;
	struct Locus *loc;
	
	if(!loki->models->tlocus || n_all<2) return -1;
	for(i=0;i<loki->params.n_tloci;i++) if(!loki->models->tlocus[i].flag) break;
	if(i==loki->params.n_tloci)	{
		if(i<loki->params.max_tloci) {
			loki->params.n_tloci*=1.5;
			if(loki->params.n_tloci>loki->params.max_tloci) loki->params.n_tloci=loki->params.max_tloci;
			if(!(loki->models->tlocus=realloc(loki->models->tlocus,sizeof(struct Locus)*loki->params.n_tloci))) ABT_FUNC(MMsg);
			k=loki->params.n_tloci-i;
			for(j=i;j<loki->params.n_tloci;j++) {
				loc=loki->models->tlocus+j;
				init_traitlocus(loki,j);
			}
		} else return -1;
	}
	j=n_all*(n_all+1)/2-1;
	k=loki->models->n_models;
	k=k*(k+1)/2;
	if(!(tpp=malloc(sizeof(void *)*(loki->models->n_models+loki->pedigree->n_genetic_groups)))) ABT_FUNC(MMsg);
	if(!(tp=malloc(sizeof(double)*(j*loki->models->n_models+k+loki->pedigree->n_genetic_groups*n_all)))) ABT_FUNC(MMsg);
	loc=loki->models->tlocus+i;
	loc->freq=tpp;
	loc->eff=loc->freq+loki->pedigree->n_genetic_groups;
	loc->variance=tp;
	tp+=k;
	for(k=0;k<loki->models->n_models;k++) {
		loc->eff[k]=tp;
		tp+=j;
	}
	for(k=0;k<loki->pedigree->n_genetic_groups;k++) loc->freq[k]=tp+n_all*k;
	loc->n_alleles=n_all;
	if(!(loc->seg[0])) {
		k=loki->pedigree->ped_size;
		if(!(loc->seg[0]=calloc((size_t)4*k,sizeof(int)))) ABT_FUNC(MMsg);
		loc->seg[1]=loc->seg[0]+k;
		loc->genes[0]=loc->seg[1]+k;
		loc->genes[1]=loc->genes[0]+k;
	}
	ff=loc->founder_flag;
	for(k=0;k<loki->pedigree->ped_size;k++) loc->seg[0][k]=loc->seg[1][k]=ff[k]?-1:-2;
	if(!loc->gt) {
		if(!(loc->gt=calloc((size_t)loki->pedigree->ped_size,sizeof(int)))) ABT_FUNC(MMsg);
	}
	if(!loc->lk_store) {
		if(!(loc->lk_store=malloc(sizeof(double)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
	}
	return i;
}

void delete_traitlocus(struct Locus *loc)
{
	loc->flag=0;
	if(loc->variance) free(loc->variance);
	if(loc->freq) free(loc->freq);
	loc->freq=0;
	loc->eff=0;
	loc->variance=0;
}

void delete_all_traitloci(struct loki *loki)
{
	int i;
	struct Locus *loc;
	
	if(loki->models->tlocus) {
		for(i=0;i<loki->params.n_tloci;i++) {
			loc=loki->models->tlocus+i;
			if(loki->models->tlocus[i].flag) delete_traitlocus(loc);
			if(loc->seg[0]) free(loc->seg[0]);
			if(loc->gt) free(loc->gt);
			if(loc->lk_store) free(loc->lk_store);
			loc->seg[0]=0;
			loc->gt=0;
			loc->lk_store=0;
		}
	}
}

void AllocEffects(struct loki *loki)
{
	int i,j,k,k1,type,mod,n_all;
	struct Variable *var;
	double *tp;
	struct Model *model;
	struct Id_Record *id_array;
	
	message(INFO_MSG,"Allocating storage for effects\n");
	id_array=loki->pedigree->id_array;
	/* Count amount of storage required for residuals, BV etc. */
	for(mod=0;mod<loki->models->n_models;mod++) {
		model=loki->models->models+mod;
		for(k=i=0;i<model->n_terms;i++) {
			model->term[i].df=1;
			if(model->term[i].n_vars>1) ABT_FUNC("Sorry - interaction terms not currently supported\n");
			type=model->term[i].vars[0].type;
			j=model->term[i].vars[0].var_index;
			if(type&ST_TRAITLOCUS) model->term[i].df=2;
			else if(type&ST_ID) model->term[i].df=loki->pedigree->ped_size;
			else if(type&ST_MARKER) {
				k1=loki->markers->marker[j].locus.n_alleles;
				model->term[i].df=k1*(k1+1)/2-1;
				loki->markers->marker[j].mterm[mod]=model->term+i;
				if(!(loki->markers->marker[j].locus.variance)) if(!(loki->markers->marker[j].locus.variance=malloc(sizeof(double)*loki->models->n_models*(loki->models->n_models+1)/2))) ABT_FUNC(MMsg);
			} else {
				if(type&ST_CONSTANT) var=loki->data->id_variable+j;
				else var=loki->data->nonid_variable+j;
				if(var->type&ST_FACTOR) {
					model->term[i].vars[0].type|=ST_FACTOR;
					if(type&ST_RANDOM) model->term[i].df=var->n_levels;
					else model->term[i].df=var->n_levels-1;
				}
			}
			if(!(type&ST_TRAITLOCUS)) k+=model->term[i].df;
		}
		if(!k) continue;
		if(!(tp=malloc(sizeof(double)*k))) ABT_FUNC(MMsg);
		loki->sys.RemBlock=AddRemem(tp,loki->sys.RemBlock);
		for(i=0;i<k;i++) tp[i]=0.0;
		for(i=0;i<model->n_terms;i++) if(!(model->term[i].vars[0].type&ST_TRAITLOCUS)) {
			model->term[i].eff=tp;
			if(model->term[i].vars[0].type&ST_MARKER) {
				j=model->term[i].vars[0].var_index;
				if(!(loki->markers->marker[j].locus.eff)) {
					loki->markers->marker[j].locus.eff=lk_malloc(sizeof(void *)*loki->models->n_models);
					for(k=0;k<loki->models->n_models;k++) loki->markers->marker[j].locus.eff[k]=0;
				}
				loki->markers->marker[j].locus.eff[mod]=tp;
			}
			tp+=model->term[i].df;
		}
	}
		if(loki->params.est_aff_freq) {
			for(j=i=0;i<loki->markers->n_markers;i++) j+=loki->markers->marker[i].locus.n_alleles;
			if(j) {
				if(!(tp=malloc(sizeof(double)*j*2))) ABT_FUNC(MMsg);
				loki->sys.RemBlock=AddRemem(tp,loki->sys.RemBlock);
				for(i=0;i<j*2;i++) tp[i]=0.0;
				for(i=0;i<loki->markers->n_markers;i++) {
					n_all=loki->markers->marker[i].locus.n_alleles;
					loki->markers->marker[i].locus.aff_freq=tp;
					loki->markers->marker[i].locus.diff_freq=tp+n_all;
					tp+=2*n_all;
				}
			}
		}
}
