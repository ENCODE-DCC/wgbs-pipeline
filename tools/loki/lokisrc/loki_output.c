/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - MSKCC                                 *
 *                                                                          *
 *                          August 2000                                     *
 *                                                                          *
 * loki_output.c:                                                           *
 *                                                                          *
 * Routines for sample output                                               *
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
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "version.h"
#include "loki_ibd.h"
#include "loki_utils.h"
#include "calc_var_locus.h"
#include "sample_rand.h"
#include "mat_utils.h"
#include "loki_output.h"

static double *tot_gen_var;

void Output_BV(FILE *fptr,struct loki *loki)
{
	int i,j,mx;
	double mu,s,z,z1;
	struct Id_Record *idr;

	mx=(int)get_max_idlen();
	z=loki->params.bv_iter?(double)loki->params.bv_iter:1.0;
	z1=z>1.0?z-1.0:1.0;
	idr=loki->pedigree->id_array;
	fputs("ID ",fptr);
	for(i=3;i<mx;i++) fputc(' ',fptr);
	fputs("         BV           SD(BV)\n",fptr);
	for(i=0;i<loki->pedigree->ped_size;i++,idr++) {
		j=(int)print_orig_id(fptr,i+1);
		for(;j<=mx;j++) fputc(' ',fptr);
		mu=idr->bvsum[0]/z;
		s=(idr->bvsum2[0]-z*mu*mu)/z1;
		fprintf(fptr," %13.8f %13.8f\n",mu,sqrt(s));
	}
}

/* Version 2.1 output routine */
static void OutputSample_1(FILE *fptr,int lp,struct loki *loki)
{
	int i,j,k,l,type,nq,nq1,grp;
	struct Model *mod;
	struct Locus *loc;
	
#ifdef DEBUG
	int mn,mn1;
	double z;
#endif
	
	nq=nq1=0;
	mod=loki->models->models;
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag)	{
		nq++;
		if(loki->models->tlocus[l].flag&TL_LINKED) nq1++;
	}
	(void)fprintf(fptr,"%d %d %d %g %g %g",lp,nq,nq1,loki->models->grand_mean[0],loki->models->residual_var[0],loki->models->tau[0]);
	if(loki->models->use_student_t) (void)fprintf(fptr," %g",loki->models->res_nu);
	if(loki->models->models[0].polygenic_flag) (void)fprintf(fptr," %g",loki->models->additive_var[0]);
	for(i=0;i<loki->models->n_random;i++) (void)fprintf(fptr," %g",loki->models->c_var[i][0]);
	if(loki->params.output_type==OUTPUT_TYPE_ORIGINAL) (void)fprintf(fptr," %d %d",loki->sys.n_cov_columns,loki->pedigree->n_genetic_groups);
	for(i=0;i<mod->n_terms;i++) if(mod->term[i].out_flag) {
		type=mod->term[i].vars[0].type;
		if(type&ST_ID) continue;
		k=mod->term[i].vars[0].var_index;
		if(type&ST_MARKER) {
			for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
				for(j=0;j<loki->markers->marker[k].locus.n_alleles-1;j++) (void)fprintf(fptr," %g",loki->markers->marker[k].locus.freq[grp][j]);
		}
		if(!(type&ST_TRAITLOCUS))
		  for(j=0;j<mod->term[i].df;j++) (void)fprintf(fptr," %g",mod->term[i].eff[j]);
		if(type&ST_MARKER) {
			loc=&loki->markers->marker[k].locus;
			calc_var_locus(loc,loki);
			(void)fprintf(fptr," %g",sqrt(loc->variance[0]));
		}
	}
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag) {
		if(loki->models->tlocus[l].flag&TL_LINKED) {
			(void)fprintf(fptr," %d",loki->models->tlocus[l].link_group+1);
			for(k=0;k<=loki->markers->sex_map;k++) (void)fprintf(fptr," %g",loki->models->tlocus[l].pos[1-k]);
		} else {
			(void)fputs(" 0",fptr);
			if(loki->params.output_type==OUTPUT_TYPE_ORIGINAL) for(k=0;k<=loki->markers->sex_map;k++) (void)fputs(" -1",fptr);
		}
		k=loki->models->tlocus[l].n_alleles;
		for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
		  for(j=0;j<k-1;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].freq[grp][j]);
		k=k*(k+1)/2-1;
		for(j=0;j<k;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].eff[0][j]);
		(void)fprintf(fptr," %g",sqrt(loki->models->tlocus[l].variance[0]));
	}
	(void)fputc('\n',fptr);
#ifdef DEBUG
	if((loki->params.debug_level)&2) {
		mn=mn1=0;
		for(i=0;i<N_MOVE_STATS;i++) {
			z=loki->sys.move_stats[i].n?(double)loki->sys.move_stats[i].success/(double)loki->sys.move_stats[i].n:0.0;
			(void)fputc(i?' ':'[',fptr);
			(void)fprintf(fptr," %g",z);
			if(i<4) {
				mn+=loki->sys.move_stats[i].success;
				mn1+=loki->sys.move_stats[i].n;
			}
		}
		z=mn1?(double)mn/(double)mn1:0.0;
		(void)fprintf(fptr," %g]\n",z);
	}
#endif
	(void)fflush(fptr);
}

/* Version 2.2 output routine */
static void OutputSample_2(FILE *fptr,int lp,struct loki *loki)
{
	int i,j,k,l,type,grp;
	struct Model *mod;
	struct Locus *loc;
#ifdef DEBUG
	int mn,mn1;
	double z;
#endif

	tot_gen_var[0]=0.0;
	mod=loki->models->models;
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag) {
		tot_gen_var[0]+=loki->models->tlocus[l].variance[0];
	}
	if(loki->models->models[0].polygenic_flag) tot_gen_var[0]+=loki->models->additive_var[0];
	for(i=0;i<mod->n_terms;i++) {
		k=mod->term[i].vars[0].var_index;
		type=mod->term[i].vars[0].type;
		if(type&ST_MARKER) {
			loc=&loki->markers->marker[k].locus;
			calc_var_locus(loc,loki);
			tot_gen_var[0]+=loc->variance[0];
		}
	}
	(void)fprintf(fptr,"%d %g %g",lp,loki->models->grand_mean[0],loki->models->residual_var[0]);
	if(loki->models->use_student_t) (void)fprintf(fptr," %g",loki->models->res_nu);
	if(loki->models->models[0].polygenic_flag) (void)fprintf(fptr," %g",loki->models->additive_var[0]);
	for(i=0;i<loki->models->n_random;i++) (void)fprintf(fptr," %g",loki->models->c_var[i][0]);
	for(i=0;i<mod->n_terms;i++) if(mod->term[i].out_flag) {
		type=mod->term[i].vars[0].type;
		if(type&ST_ID) continue;
		k=mod->term[i].vars[0].var_index;
		if(type&ST_MARKER) {
			for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
				for(j=0;j<loki->markers->marker[k].locus.n_alleles-1;j++) (void)fprintf(fptr," %g",loki->markers->marker[k].locus.freq[grp][j]);
		}
		if(!(type&ST_TRAITLOCUS))
		  for(j=0;j<mod->term[i].df;j++) (void)fprintf(fptr," %g",mod->term[i].eff[j]);
		if(type&ST_MARKER) {
			(void)fprintf(fptr," %g",sqrt(loki->markers->marker[k].locus.variance[0]));
		}
	}
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag) {
		if(loki->models->tlocus[l].flag&TL_LINKED) {
			(void)fprintf(fptr," %d",loki->models->tlocus[l].link_group+1);
			for(k=0;k<=loki->markers->sex_map;k++) (void)fprintf(fptr," %g",loki->models->tlocus[l].pos[1-k]);
		} else (void)fputs(" 0",fptr);
		k=loki->models->tlocus[l].n_alleles;
		for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
		  for(j=0;j<k-1;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].freq[grp][j]);
		k=k*(k+1)/2-1;
		for(j=0;j<k;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].eff[0][j]);
		(void)fprintf(fptr," %g",sqrt(loki->models->tlocus[l].variance[0]));
	}
	(void)fputc('\n',fptr);
#ifdef DEBUG
	if((loki->params.debug_level)&2) {
		mn=mn1=0;
		for(i=0;i<N_MOVE_STATS;i++) {
			z=loki->sys.move_stats[i].n?(double)loki->sys.move_stats[i].success/(double)loki->sys.move_stats[i].n:0.0;
			(void)fputc(i?' ':'[',fptr);
			(void)fprintf(fptr,"%g",z);
			if(i<4) {
				mn+=loki->sys.move_stats[i].success;
				mn1+=loki->sys.move_stats[i].n;
			}
		}
		z=mn1?(double)mn/(double)mn1:0.0;
		(void)fprintf(fptr," %g]\n",z);
	}
#endif
	(void)fflush(fptr);
}

/* Version 2.3 output routine */
static void OutputSample_3(FILE *fptr,int lp,struct loki *loki)
{
	int i,j,k,l,type,grp,mod,mod1;
	int n_models;
#ifdef DEBUG
	int mn,mn1;
	double z;
#endif

	(void)fprintf(fptr,"%d",lp);
	n_models=loki->models->n_models;
	for(j=mod=0;mod<n_models;mod++) {
		for(k=0;k<mod;k++) (void)fprintf(fptr," %g",loki->models->residual_var[j++]);
		if(loki->models->res_var_set[mod]!=1) (void)fprintf(fptr," %g",loki->models->residual_var[j]);
		j++;
		for(i=0;i<loki->models->models[mod].n_terms;i++) {
			k=loki->models->models[mod].term[i].vars[0].var_index;
			type=loki->models->models[mod].term[i].vars[0].type;
			if(type&ST_MARKER) calc_var_locus(&loki->markers->marker[k].locus,loki);
		}
	}
	k=n_models*(n_models+1)/2;
	for(i=0;i<k;i++) tot_gen_var[i]=0.0;
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag) {
		for(i=0;i<k;i++) tot_gen_var[i]+=loki->models->tlocus[l].variance[i];
	}
	for(l=0;l<n_models;l++) if(loki->models->models[l].polygenic_flag) break;
	if(l<n_models) for(i=0;i<k;i++) tot_gen_var[i]+=loki->models->additive_var[i];
	for(j=0;j<loki->markers->n_markers;j++) if(loki->markers->marker[j].locus.variance) for(i=0;i<k;i++) tot_gen_var[i]+=loki->markers->marker[j].locus.variance[i];
	if(loki->models->use_student_t) (void)fprintf(fptr," %g",loki->models->res_nu);
	for(i=0;i<k;i++) (void)fprintf(fptr," %g",tot_gen_var[i]);
	if(l<n_models) {
		for(mod=0;mod<n_models;mod++) if(loki->models->models[mod].polygenic_flag) {
			for(mod1=0;mod1<mod;mod1++) if(loki->models->models[mod1].polygenic_flag) (void)fprintf(fptr," %g",BB(loki->models->additive_var,mod,mod1));
			if(loki->models->add_var_set[mod]!=1) (void)fprintf(fptr," %g",BB(loki->models->additive_var,mod,mod));
		}
	}
	for(i=0;i<loki->models->n_random;i++) {
		for(mod=0;mod<n_models;mod++) if(loki->models->rand_flag[mod]&(1<<i)) {
			for(mod1=0;mod1<=mod;mod1++) if(loki->models->rand_flag[mod1]&(1<<i)) (void)fprintf(fptr," %g",BB(loki->models->c_var[i],mod,mod1));
		}
	}
	for(mod=0;mod<n_models;mod++) {
		if(loki->models->grand_mean_set[mod]!=1) (void)fprintf(fptr," %g",loki->models->grand_mean[mod]);
		for(i=0;i<loki->models->models[mod].n_terms;i++) if(loki->models->models[mod].term[i].out_flag) {
			type=loki->models->models[mod].term[i].vars[0].type;
			if(type&ST_ID) continue;
			k=loki->models->models[mod].term[i].vars[0].var_index;
			if(!(type&(ST_MARKER|ST_TRAITLOCUS))) {
				for(j=0;j<loki->models->models[mod].term[i].df;j++) (void)fprintf(fptr," %g",loki->models->models[mod].term[i].eff[j]);
			}
		}
	}
	for(k=0;k<loki->markers->n_markers;k++) if(loki->markers->marker[k].locus.variance) {
		for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
		  for(j=0;j<loki->markers->marker[k].locus.n_alleles-1;j++) (void)fprintf(fptr," %g",loki->markers->marker[k].locus.freq[grp][j]);
		for(mod=0;mod<n_models;mod++) if(loki->markers->marker[k].mterm[mod]) {
			for(j=0;j<loki->markers->marker[k].mterm[mod]->df;j++) (void)fprintf(fptr," %g",loki->markers->marker[k].mterm[mod]->eff[j]);
		}
		for(mod=0;mod<n_models;mod++) if(loki->markers->marker[k].mterm[mod]) 
		  for(mod1=0;mod1<=mod;mod1++) if(loki->markers->marker[k].mterm[mod1]) (void)fprintf(fptr," %g",sqrt(BB(loki->markers->marker[k].locus.variance,mod,mod1)));
	}
	for(l=0;l<loki->params.n_tloci;l++) if(loki->models->tlocus[l].flag) {
		if(loki->models->tlocus[l].flag&TL_LINKED) {
			(void)fprintf(fptr," %d",loki->models->tlocus[l].link_group+1);
			for(k=0;k<=loki->markers->sex_map;k++) (void)fprintf(fptr," %g",loki->models->tlocus[l].pos[1-k]);
		} else (void)fputs(" 0",fptr);
		k=loki->models->tlocus[l].n_alleles;
		for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++)
		  for(j=0;j<k-1;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].freq[grp][j]);
		k=k*(k+1)/2-1;
		for(mod=0;mod<n_models;mod++) if(loki->models->tlocus[l].model_flag&(1<<mod))
		  for(j=0;j<k;j++) (void)fprintf(fptr," %g",loki->models->tlocus[l].eff[mod][j]);
		for(mod=0;mod<n_models;mod++) if(loki->models->tlocus[l].model_flag&(1<<mod))
		  for(mod1=0;mod1<=mod;mod1++) if(loki->models->tlocus[l].model_flag&(1<<mod1)) (void)fprintf(fptr," %g",sqrt(BB(loki->models->tlocus[l].variance,mod,mod1)));
	}
	(void)fputc('\n',fptr);
#ifdef DEBUG
	if((loki->params.debug_level)&2) {
		mn=mn1=0;
		for(i=0;i<N_MOVE_STATS;i++) {
			z=loki->sys.move_stats[i].n?(double)loki->sys.move_stats[i].success/(double)loki->sys.move_stats[i].n:0.0;
			(void)fputc(i?' ':'[',fptr);
			(void)fprintf(fptr,"%g",z);
			if(i<4) {
				mn+=loki->sys.move_stats[i].success;
				mn1+=loki->sys.move_stats[i].n;
			}
		}
		z=mn1?(double)mn/(double)mn1:0.0;
		(void)fprintf(fptr," %g]\n",z);
	}
#endif
	(void)fflush(fptr);
}

static void free_outputsample(void) 
{
	OutputSample(0,0,0);
}

void OutputSample(FILE *fptr,int lp,struct loki *loki)
{
	if(!fptr) {
		if(tot_gen_var) {
			free(tot_gen_var);
			tot_gen_var=0;
		}
		return;
	}
	if(!tot_gen_var && loki->models->n_models) {
		if(!(tot_gen_var=malloc(sizeof(double)*loki->models->n_models*(loki->models->n_models+1)/2))) ABT_FUNC(MMsg);
		if(atexit(free_outputsample)) message(WARN_MSG,"Unable to register exit function free_outputsample()\n");
	}
	switch(loki->params.output_type) {
	 case OUTPUT_TYPE_ORIGINAL:
	 case OUTPUT_VERSION_2_1: 
		OutputSample_1(fptr,lp,loki);
		break;
	 case OUTPUT_VERSION_2_2:
		OutputSample_2(fptr,lp,loki);
	 case OUTPUT_VERSION_2_3:
		OutputSample_3(fptr,lp,loki);
		break;
	}
}

void OutputFreqHeader(FILE *fptr,struct loki *loki)
{
	int i,j;
	
	(void)fprintf(fptr,"Created by %s: %s",LOKI_NAME,ctime(&loki->sys.lktime.start_time));
	for(i=0;i<loki->markers->n_markers;i++) {
		j=loki->markers->marker[i].locus.link_group;
		if(loki->markers->linkage[j].type&LINK_PSEUDO) continue;
		(void)fprintf(fptr,"%d ",i);
		print_marker_name(fptr,i);
		(void)fputc(':',fptr);
		for(j=0;j<loki->markers->marker[i].locus.n_alleles-1;j++) {
			(void)fputc(' ',fptr);
			print_allele_type1(fptr,i,j);
		}
		(void)fputs(" (",fptr);
		print_allele_type1(fptr,i,j);
		(void)fputs(")\n",fptr);
	}
	(void)fprintf(fptr,"No. genetic groups: %d\n",loki->pedigree->n_genetic_groups);
	if(loki->params.est_aff_freq) (void)fputs("Estimating allele frequencies amongst affecteds\n",fptr);
	(void)fputs("--------------------\n",fptr);
	(void)fflush(fptr);
}

void OutputFreq(FILE *fptr,int lp,struct loki *loki)
{
	int i,j,k;
	
	(void)fprintf(fptr,"%d ",lp);
	for(i=0;i<loki->markers->n_markers;i++) {
		for(j=0;j<loki->markers->marker[i].locus.n_alleles-1;j++) {
			for(k=0;k<loki->pedigree->n_genetic_groups;k++) (void)fprintf(fptr,"%g ",loki->markers->marker[i].locus.freq[k][j]);
			if(loki->params.est_aff_freq) {
				(void)fprintf(fptr,"%g ",loki->markers->marker[i].locus.aff_freq[j]);
			}
		}
	}
	(void)fputc('\n',fptr);
	(void)fflush(fptr);
}

void OutputHeader(FILE *fptr,struct loki *loki)
{
	int i,j,k,k1,k2,k3,type,n_all,grp,mod,mod1;
	struct Variable *var;
	struct Locus **locilist=0;
	char *s;
	
	(void)fprintf(fptr,"Created by %s: %s",LOKI_NAME,ctime(&loki->sys.lktime.start_time));
	if(!loki->params.analysis) {
		(void)fprintf(fptr,"Output format: %d\n",loki->params.output_type);
		for(mod=0;mod<loki->models->n_models;mod++) {
			type=loki->models->models[mod].var.type;
			if(type) {
				if(loki->models->n_models>1) {
					if(!mod) (void)fputs("Models: \n  ",fptr);
					else (void)fputs("  ",fptr);
				} else (void)fputs("Model: ",fptr);
				i=loki->models->models[mod].var.var_index;
				var=(type&ST_CONSTANT)?loki->data->id_variable+i:loki->data->nonid_variable+i;
				if(type&ST_CENSORED) fputc('(',fptr);
				(void)fputs(var->name,fptr);
				if(type&ST_CENSORED) fputc(')',fptr);
				(void)fputs(" = ",fptr);
				for(i=0;i<loki->models->models[mod].n_terms;i++) {
					type=loki->models->models[mod].term[i].vars[0].type;
					if(i) (void)fputs(" + ",fptr);
					if(type&ST_TRAITLOCUS) (void)fputs("QTL",fptr);
					else if(type&ST_ID) (void)fputs("ID\'",fptr);
					else {
						j=loki->models->models[mod].term[i].vars[0].var_index;
						if(type&ST_MARKER) {
							print_marker_name(fptr,j);
						} else {
							var=(type&ST_CONSTANT)?loki->data->id_variable+j:loki->data->nonid_variable+j;
							(void)fputs(var->name,fptr);
							if(type&ST_RANDOM) (void)fputc('\'',fptr);
						}
					}
				}
				(void)fputc('\n',fptr);
			}
		}
	}
	if(!(loki->params.analysis&NULL_ANALYSIS)) {
		if(loki->markers->n_links) {
			(void)fprintf(fptr,"Input map function: %s\nOutput map function: Haldane\n",loki->params.map_function?"Kosambi":"Haldane");
			if((loki->params.n_tloci+loki->markers->n_markers) && !(locilist=malloc(sizeof(void *)*(loki->params.n_tloci+loki->markers->n_markers)))) ABT_FUNC(MMsg);
			(void)fputs("Linkage groups:\n",fptr);
			for(i=0;i<loki->markers->n_links;i++) {
				s=loki->markers->linkage[i].name;
				if(!s) s="<UN-NAMED>";
				(void)fprintf(fptr,"  %d: %s  Map range: ",i+1,s);
				for(k3=0;k3<=loki->markers->sex_map;k3++) (void)fprintf(fptr," (%gcM to %gcM) ",loki->markers->linkage[i].r1[1-k3],loki->markers->linkage[i].r2[1-k3]);
				(void)fputc('\n',fptr);
				if(locilist) {
					get_locuslist(locilist,i,&k2,1);
					gnu_qsort(locilist,(size_t)k2,sizeof(void *),cmp_loci);
					for(k1=0;k1<k2;k1++) {
						(void)fprintf(fptr,"    %s -",loki->markers->marker[locilist[k1]->index].name);
						for(k3=0;k3<=loki->markers->sex_map;k3++) (void)fprintf(fptr," %g",locilist[k1]->pos[1-k3]);
						(void)fputc('\n',fptr);
					}
				}
				switch(loki->markers->linkage[i].ibd_est_type) {
				 case IBD_EST_DISCRETE:
					if(!loki->markers->linkage[i].ibd_list) ABT_FUNC("Internal error - no ibd list\n");
					(void)fprintf(fptr,"IBD Matrix estimated at position%s",loki->markers->linkage[i].ibd_list->idx==1?":":"s:");
					for(j=0;j<loki->markers->linkage[i].ibd_list->idx;j++) {
						(void)fputc(j?',':' ',fptr);
						(void)fprintf(fptr,"%g",loki->markers->linkage[i].ibd_list->pos[j]);
					}
					(void)fputc('\n',fptr);
					break;
				 case IBD_EST_GRID:
					if(!loki->markers->linkage[i].ibd_list) ABT_FUNC("Internal error - no ibd list\n");
					(void)fprintf(fptr,"IBD Matrix estimated at grid of positions from %g to %g step %g\n",
									  loki->markers->linkage[i].ibd_list->pos[0],loki->markers->linkage[i].ibd_list->pos[1],loki->markers->linkage[i].ibd_list->pos[2]);
					break;
				 case IBD_EST_MARKERS:
					(void)fputs("IBD Matrix estimated at all marker locations\n",fptr);
					break;
				}
			}
			if(locilist) free(locilist);
			(void)fputs("Total Map Length:",fptr);
			for(k3=0;k3<=loki->markers->sex_map;k3++) (void)fprintf(fptr," %gcM",loki->markers->total_maplength[1-k3]);
			(void)fputc('\n',fptr);		}
	}
	if(!loki->params.analysis && loki->models->models) {
		(void)fputs("Output columns:\n",fptr);
		k1=0;
		(void)fprintf(fptr,"  %d: Iteration count\n",++k1);
		if(loki->params.output_type<=OUTPUT_VERSION_2_1) {
			(void)fprintf(fptr,"  %d: No. QTL's in models[0]\n",++k1);
			(void)fprintf(fptr,"  %d: No. linked QTL's\n",++k1);
		}
		if(loki->params.output_type<3) {
			if(loki->models->n_models>1) {
				for(mod=0;mod<loki->models->n_models;mod++) (void)fprintf(fptr,"  %d: Grand mean %d\n",++k1,mod+1);
			} else (void)fprintf(fptr,"  %d: Grand mean\n",++k1);
		}
		if(loki->models->n_models>1) {
			for(k=mod=0;mod<loki->models->n_models;mod++) {
				for(mod1=0;mod1<mod;mod1++) {
					if(loki->models->res_var_set[k++]!=1) (void)fprintf(fptr,"  %d: residual covariance %d %d\n",++k1,mod1+1,mod+1);
				}
				if(loki->models->res_var_set[k++]!=1) (void)fprintf(fptr,"  %d: residual variance   %d\n",++k1,mod+1);
			}
		} else {
			if(loki->models->res_var_set[0]!=1) (void)fprintf(fptr,"  %d: residual variance\n",++k1);
		}
		if(loki->models->use_student_t) (void)fprintf(fptr,"  %d: d.f. for student t distribution of loki->models->residuals\n",++k1);
		if(loki->params.output_type<=OUTPUT_VERSION_2_1) {
			if(loki->models->n_models>1) {
				for(mod=0;mod<loki->models->n_models;mod++) (void)fprintf(fptr,"  %d: tau %d\n",++k1,mod+1);
			} else (void)fprintf(fptr,"  %d: loki->models->tau\n",++k1);
		}
		if(loki->params.output_type>=OUTPUT_VERSION_2_3) {
			if(loki->models->n_models>1) {
				for(mod=0;mod<loki->models->n_models;mod++) {
					for(mod1=0;mod1<mod;mod1++) {
						(void)fprintf(fptr,"  %d: Total genetic covariance %d %d\n",++k1,mod1+1,mod+1);
					}
					(void)fprintf(fptr,"  %d: Total genetic variance   %d\n",++k1,mod+1);
				}
			} else {
				(void)fprintf(fptr,"  %d: Total genetic variance\n",++k1);
			}
		} 
		if(loki->models->n_models>1) {
			for(k=mod=0;mod<loki->models->n_models;mod++) if(loki->models->models[mod].polygenic_flag) {
				for(mod1=0;mod1<mod;mod1++) if(loki->models->models[mod1].polygenic_flag) {
					if(loki->models->add_var_set[k++]!=1) (void)fprintf(fptr,"  %d: additive covariance %d %d\n",++k1,mod1+1,mod+1);
				}
				if(loki->models->add_var_set[k++]!=1) (void)fprintf(fptr,"  %d: additive variance   %d\n",++k1,mod+1);
			}
		} else if(loki->models->models[0].polygenic_flag) {
			if(loki->models->add_var_set[0]!=1) (void)fprintf(fptr,"  %d: additive variance\n",++k1);
		}
		for(j=0;j<loki->models->n_random;j++) {
			var=loki->models->rand_list[j];
			for(mod=0;mod<loki->models->n_models;mod++) if(loki->models->rand_flag[mod]&(1<<j)) {
				for(mod1=0;mod1<mod;mod1++) if(loki->models->rand_flag[mod]&(1<<j)) {
					(void)fprintf(fptr,"  %d: Additional random covariance for %s",++k1,var->name);
					if(var->index) (void)fprintf(fptr,"(%d)",var->index);
					(void)fprintf(fptr," %d %d\n",mod1+1,mod+1);
				}
				(void)fprintf(fptr,"  %d: Additional random variance for %s",++k1,var->name);
				if(var->index) (void)fprintf(fptr,"(%d)",var->index);
				if(loki->models->n_models>1) (void)fprintf(fptr,"   %d\n",mod+1);
				else (void)fputc('\n',fptr);
			}
		}
		if(loki->params.output_type==OUTPUT_TYPE_ORIGINAL) {
			(void)fprintf(fptr,"  %d: No. covariate columns\n",++k1);
			(void)fprintf(fptr,"  %d: No. genetic groups\n",++k1);
		}
		k3=k1;
		for(mod=0;mod<loki->models->n_models;mod++) {
			(void)fputs(" covariate data",fptr);
			if(loki->models->n_models>1) (void)fprintf(fptr," - model %d:\n",mod+1);
			else (void)fputs(":\n",fptr);
			if(loki->params.output_type==DEFAULT_OUTPUT_TYPE) {
				if(loki->models->grand_mean_set[mod]!=1) (void)fprintf(fptr,"  %d: Grand mean\n",++k1);
			}
			for(i=0;i<loki->models->models[mod].n_terms;i++) if(loki->models->models[mod].term[i].out_flag) {
				type=loki->models->models[mod].term[i].vars[0].type;
				if(type&ST_ID) continue;
				k=loki->models->models[mod].term[i].vars[0].var_index;
				if(type&(ST_TRAITLOCUS|ST_MARKER)) continue;
				if(type&ST_CONSTANT) var=loki->data->id_variable+k;
				else var=loki->data->nonid_variable+k;
				k=(type&ST_RANDOM)?0:1;
				for(j=0;j<loki->models->models[mod].term[i].df;j++) {
					(void)fprintf(fptr,"  %d: %s",++k1,var->name);
					(void)fputs(" effect ",fptr);
					if(var->type&ST_FACTOR)	{
						if(var->rec_flag==ST_STRING) (void)fputs(var->recode[j+k].string,fptr);
						else (void)fprintf(fptr,"%d",var->recode[j+k].value);
					}
					(void)fputc('\n',fptr);
				}
			}
		}
		for(k=0;k<loki->markers->n_markers;k++) if(loki->markers->marker[k].locus.variance) {
			n_all=loki->markers->marker[k].locus.n_alleles;
			for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++) {
				for(j=0;j<n_all-1;j++) {
					(void)fprintf(fptr,"  %d: %s",++k1,loki->markers->marker[k].name);
					if(loki->pedigree->group_recode.recode) {
						(void)fputs(" [Genetic group ",fptr);
						if(loki->pedigree->group_recode.flag==ST_STRING) (void)fputs(loki->pedigree->group_recode.recode[grp].string,fptr);
						else (void)fprintf(fptr,"%d",loki->pedigree->group_recode.recode[grp].value);
						(void)fputc(']',fptr);
					}
					(void)fputs(" freq. ",fptr);
					print_allele_type1(fptr,k,j);
					(void)fputc('\n',fptr);
				}
			}
			for(mod=0;mod<loki->models->n_models;mod++) if(loki->markers->marker[k].mterm[mod]) {
				for(j=1;j<n_all;j++) for(k2=0;k2<=j;k2++)	{
					(void)fprintf(fptr,"  %d: %s",++k1,loki->markers->marker[k].name);
					(void)fputs(" effect ",fptr);
					print_allele_type1(fptr,k,k2);
					(void)fputc(',',fptr);
					print_allele_type1(fptr,k,j);
					if(loki->models->n_models>1) (void)fprintf(fptr," for model %d ",mod+1);
					(void)fputc('\n',fptr);
				}
			}
			for(mod=0;mod<loki->models->n_models;mod++) if(loki->markers->marker[k].mterm[mod]) 
			  for(mod1=0;mod1<=mod;mod1++) if(loki->markers->marker[k].mterm[mod1]) {
				  (void)fprintf(fptr,"  %d: %s",++k1,loki->markers->marker[k].name);
				  (void)fputs(" size ",fptr);
				  if(loki->models->n_models>1) {
					  if(mod==mod1) (void)fprintf(fptr,"model %d",mod+1);
					  else (void)fprintf(fptr,"models %d,%d\n",mod+1,mod1+1);
				  }
				  (void)fputc('\n',fptr);
			  }
		}
		if(loki->models->tlocus && loki->params.max_tloci) {
			(void)fputs(" QTL data blocks:\n",fptr);
			if(loki->params.output_type==OUTPUT_TYPE_ORIGINAL) {
 				if(loki->markers->sex_map) (void)fputs("  male position\n  female position\n",fptr);
				else (void)fputs("  position\n",fptr);
			} else {
 				if(loki->markers->sex_map) (void)fputs("  [male position if linked]\n  [female position if linked]\n",fptr);
				else (void)fputs("  [position if linked]\n",fptr);
			}
			for(grp=0;grp<loki->pedigree->n_genetic_groups;grp++) {
				if(loki->pedigree->group_recode.recode) {
					(void)fputs("  [Genetic group ",fptr);
					if(loki->pedigree->group_recode.flag==ST_STRING) (void)fputs(loki->pedigree->group_recode.recode[grp].string,fptr);
					else (void)fprintf(fptr,"%d",loki->pedigree->group_recode.recode[grp].value);
					(void)fputc(']',fptr);
				}
				(void)fputs("  freq. 1\n",fptr);
			}
			if(loki->models->n_models>1) {
				(void)fputs("  Effect block:\n   Model no.\n   effect 1,2\n   effect 2,2\n",fptr);
				(void)fputs("  Size block:\n",fptr);
			} else (void)fputs("  effect 1,2\n  effect 2,2\n  size\n",fptr);
			if(loki->params.min_tloci==loki->params.max_tloci) (void)fprintf(fptr,"No. QTL: %d\n",loki->params.max_tloci);
			else {
				(void)fprintf(fptr,"Number of QTL: %d to %d\n",loki->params.min_tloci,loki->params.max_tloci);
				if(loki->models->tloci_mean_set) (void)fprintf(fptr,"Mean of poisson prior on QTL number: %g\n",loki->models->tloci_mean);
			}
		}
		if(loki->models->res_var_set[0]==1) (void)fprintf(fptr,"residual variance: %g\n",loki->models->residual_var[0]);
		if(loki->models->grand_mean_set[0]==1) (void)fprintf(fptr,"Grand mean: %g\n",loki->models->grand_mean[0]);
		(void)fprintf(fptr,"tau Mode: %d\ntau Beta: %g\n",loki->models->tau_mode,loki->models->tau_beta[0]);
		(void)fprintf(fptr,"No. fixed output columns: %d\n",k1);
		if(loki->models->no_overdominant) (void)fputs("Over-dominant QTLs not allowed\n",fptr);
		k1-=k3;
		loki->sys.n_cov_columns=k1;
	} else if(!(loki->params.analysis&ESTIMATE_IBD)) {
		(void)fputs("Affected only ",fptr);
		if(loki->params.analysis&IBD_ANALYSIS) (void)fputs("IBD analysis\n",fptr);
		else (void)fputc('\n',fptr);
	}
	if(loki->params.lm_ratio>0.0) (void)fprintf(fptr,"LM ratio: %g\n",loki->params.lm_ratio);
	(void)fprintf(fptr,"SI_mode: %d\n",loki->params.si_mode);
	(void)fprintf(fptr,"No. genetic groups: %d\n",loki->pedigree->n_genetic_groups);
	if(loki->markers->sex_map) (void)fputs("Sex specific map\n",fptr);
	(void)fputs("--------------------\n",fptr);
	(void)fflush(fptr);
}

