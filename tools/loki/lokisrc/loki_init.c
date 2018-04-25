/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - Rockefeller University                         *
 *                                                                          *
 *                       October 1997                                       *
 *                                                                          *
 * loki_init.c:                                                             *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include "ranlib.h"

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "loki_output.h"
#include "loki_ibd.h"
#include "mat_utils.h"
#include "sample_rand.h"

void InitValues(struct loki *loki)
{
	int i,j,k,k1,type,n,mod,n_rec,*tmp;
	double y,s,mu;
	struct id_data *data;
	struct Id_Record *id;
	
	message(INFO_MSG,"Initializing parameter values\n");
	if(loki->sys.syst_var[SYST_NO_OVERDOMINANT].flag && loki->sys.syst_var[SYST_NO_OVERDOMINANT].data.value) loki->models->no_overdominant=1;
	if(loki->sys.syst_var[SYST_TAU_MODE].flag) {
		if(loki->sys.syst_var[SYST_TAU_MODE].flag==ST_REAL) loki->models->tau_mode=(int)loki->sys.syst_var[SYST_TAU_MODE].data.rvalue;
		else loki->models->tau_mode=loki->sys.syst_var[SYST_TAU_MODE].data.value;
	}
	if(loki->sys.syst_var[SYST_CENSOR_MODE].flag && loki->sys.syst_var[SYST_CENSOR_MODE].data.value) loki->models->censor_mode=1;
	for(mod=0;mod<loki->models->n_models;mod++) {
		type=loki->models->models[mod].var.type;
		j=loki->models->models[mod].var.var_index;
		n=0;
		s=mu=0.0;
		if(!(type&ST_FACTOR)) {
			id=loki->pedigree->id_array;
			for(i=0;i<loki->pedigree->ped_size;i++,id++) {
				if(loki->models->models[mod].polygenic_flag) id->bv[mod]=id->bvsum[mod]=id->bvsum2[mod]=0.0;
				if(type&ST_CONSTANT) id->n_rec=id->data?1:0;
				n_rec=id->n_rec;
				for(k1=k=0;k<n_rec;k++) {
					data=(type&ST_CONSTANT)?id->data+j:id->data1[k]+j;
					if(data && data->flag) {
						if(data->flag&ST_INTTYPE) y=(double)data->data.value;
						else y=data->data.rvalue;
						mu+=y;
						n++;
						k1++;
					}
				}
				if(!k1) id->res[mod]=0;
			}
			if(n) {
				if(loki->models->grand_mean_set[mod]) mu=loki->models->grand_mean[mod];
				else mu/=(double)n;
				id=loki->pedigree->id_array;
				for(i=0;i<loki->pedigree->ped_size;i++,id++) if(id->res[mod]) {
					if(type&ST_CONSTANT) n_rec=id->data?1:0;
					else n_rec=id->n_rec;
					for(k=0;k<n_rec;k++) {
						data=(type&ST_CONSTANT)?id->data+j:id->data1[k]+j;
						if(data && data->flag&ST_INTTYPE) y=(double)data->data.value;
						else y=data->data.rvalue;
						y-=mu;
						s+=y*y;
						id->res[mod][k]=y;
						if(loki->models->use_student_t) id->vv[mod][k]=1.0;
						if((type&ST_CENSORED)&&(data->flag&2)) {
							id->pseudo_qt[mod][k]=0.0;
							loki->models->censored_flag=1;
						}
					}
				}
			}
			if(!loki->models->res_var_set[mod]) {
				s=n?s/(double)n:1.0;
				if(s<loki->models->residual_var_limit[mod]) s=loki->models->residual_var_limit[mod];
				BB(loki->models->residual_var,mod,mod)=s;
			} else if(BB(loki->models->residual_var,mod,mod)<loki->models->residual_var_limit[mod]) {
				BB(loki->models->residual_var,mod,mod)=loki->models->residual_var_limit[mod];
				(void)fprintf(stderr,"Warning - residual variance for model %d reset to limit (%g)\n",mod+1,loki->models->residual_var_limit[mod]);
			}
			if(loki->models->add_var_set) {
				if(!loki->models->add_var_set[mod]) {
					BB(loki->models->additive_var,mod,mod)=BB(loki->models->residual_var,mod,mod);
					if(BB(loki->models->additive_var,mod,mod)<loki->models->additive_var_limit[mod]) BB(loki->models->additive_var,mod,mod)=loki->models->additive_var_limit[mod];
				} else if(loki->models->additive_var[mod]<loki->models->additive_var_limit[mod]) {
					loki->models->additive_var[mod]=loki->models->additive_var_limit[mod];
					(void)fprintf(stderr,"Warning - additive variance for model %d reset to limit (%g)\n",mod+1,loki->models->additive_var_limit[mod]);
				}
			}
		}
		loki->models->tau_beta[mod]=2.0;
		if(!loki->models->grand_mean_set[mod]) loki->models->grand_mean[mod]=mu;
		if(loki->sys.syst_var[SYST_TAU].flag) {
			if(loki->sys.syst_var[SYST_TAU].flag==ST_REAL) loki->models->tau_beta[0]=loki->sys.syst_var[SYST_TAU].data.rvalue;
			else loki->models->tau_beta[mod]=(double)loki->sys.syst_var[SYST_TAU].data.value;
		} else if(loki->models->tau_mode==1) loki->models->tau_beta[mod]=loki->models->residual_var[mod];
		switch(loki->models->tau_mode) {
		 case 0:
			loki->models->tau_beta[mod]*=loki->models->residual_var[mod];
		 case 1:
			loki->models->tau[mod]=loki->models->tau_beta[mod];
			break;
		 case 2:
			loki->models->tau[mod]=loki->models->residual_var[mod]*loki->models->tau_beta[mod];
			break;
		}
	}
	if(loki->sys.syst_var[SYST_IBD_OUTPUT].flag) {
		if(loki->sys.syst_var[SYST_IBD_OUTPUT].flag==ST_REAL) loki->params.ibd_mode=(int)loki->sys.syst_var[SYST_IBD_OUTPUT].data.rvalue;
		else loki->params.ibd_mode=loki->sys.syst_var[SYST_IBD_OUTPUT].data.value;
	}
	if(loki->markers->n_markers>1 && ((loki->params.ibd_mode)&IBD_SINGLE_POINT)) { /* Set up singlepoint IBD analysis */
		if(loki->params.analysis&ESTIMATE_IBD) { /* Set up singlepoint IBD analysis */
			if(!(tmp=malloc(sizeof(int)*loki->markers->n_markers))) ABT_FUNC(MMsg);
			for(i=0;i<loki->markers->n_markers;i++) {
				j=loki->markers->marker[i].locus.link_group;
				tmp[i]=loki->markers->linkage[j].type;
			}
			if(loki->markers->linkage) {
				for(i=0;i<loki->markers->n_links;i++) {
					if(loki->markers->linkage[i].ibd_list) {
						free(loki->markers->linkage[i].ibd_list->pos);
						free(loki->markers->linkage[i].ibd_list);
					}
				}
				free(loki->markers->linkage);
			}
			if(!(loki->markers->linkage=malloc(sizeof(struct Link)*loki->markers->n_markers))) ABT_FUNC(MMsg);
			loki->markers->n_links=loki->markers->n_markers;
			for(i=0;i<loki->markers->n_markers;i++) {
				loki->markers->linkage[i].name=loki->markers->marker[i].name;
				for(j=0;j<2;j++) loki->markers->linkage[i].r1[j]=loki->markers->linkage[i].r2[j]=0.0;
				loki->markers->linkage[i].n_markers=1;
				loki->markers->linkage[i].type=tmp[i];
				loki->markers->linkage[i].ibd_list=0;
				loki->markers->linkage[i].sample_pos=0;
				loki->markers->linkage[i].range_set[0]=loki->markers->linkage[i].range_set[1]=0;
				loki->markers->linkage[i].ibd_est_type=IBD_EST_MARKERS;
				loki->markers->marker[i].locus.link_group=i;
			}
			free(tmp);
		}
	}
	if(loki->sys.syst_var[SYST_SI_MODE].flag) {
		if(loki->sys.syst_var[SYST_SI_MODE].flag==ST_REAL) loki->params.si_mode=(int)loki->sys.syst_var[SYST_SI_MODE].data.rvalue;
		else loki->params.si_mode=loki->sys.syst_var[SYST_SI_MODE].data.value;
	}
	if(loki->sys.syst_var[SYST_RNG].flag) {
		if(loki->sys.syst_var[SYST_RNG].flag==ST_REAL) i=(int)loki->sys.syst_var[SYST_RNG].data.rvalue;
		else i=loki->sys.syst_var[SYST_RNG].data.value;
		if(set_mt_idx(i)<0) {
			fprintf(stderr,"Invalid RNG (%d) selected - using default generator\n",i);
		}
	}
	if(loki->sys.syst_var[SYST_GENV_OUT].flag && loki->sys.syst_var[SYST_GENV_OUT].data.value) loki->params.genv_out=1; /* EWD */
#ifdef DEBUG
	if(loki->sys.syst_var[SYST_DEBUG_LEVEL].flag)	{
		if(loki->sys.syst_var[SYST_DEBUG_LEVEL].flag==ST_REAL) loki->params.debug_level=(int)loki->sys.syst_var[SYST_DEBUG_LEVEL].data.rvalue;
		else loki->params.debug_level=loki->sys.syst_var[SYST_DEBUG_LEVEL].data.value;
	}
#endif
#ifdef TRACE_PEEL
	if(loki->sys.syst_var[SYST_PEEL_TRACE].flag) {
		if(loki->sys.syst_var[SYST_PEEL_TRACE].flag==ST_REAL) *peel_trace=(int)loki->sys.syst_var[SYST_PEEL_TRACE].data.rvalue;
		else *peel_trace=loki->sys.syst_var[SYST_PEEL_TRACE].data.value;
	}
#endif
	if(loki->sys.syst_var[SYST_LM_RATIO].flag)	{
		if(loki->sys.syst_var[SYST_LM_RATIO].flag==ST_REAL) loki->params.lm_ratio=loki->sys.syst_var[SYST_LM_RATIO].data.rvalue;
		else loki->params.lm_ratio=(double)loki->sys.syst_var[SYST_LM_RATIO].data.value;
	}
	for(i=0;i<N_MOVE_STATS;i++) loki->sys.move_stats[i].success=loki->sys.move_stats[i].n=0;
	init_rand(loki);
}
