/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                            May 2003                                      *
 *                                                                          *
 * set_recode.c:                                                            *
 *                                                                          *
 * Perform set recoding                                                     *  
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
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
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "gen_elim.h"

#define TMEM_SIZE 2048

struct temp_mem {
	struct temp_mem *next;
	int p[TMEM_SIZE];
	int i;
};

/* Set recoding from gtype sets.  Only works for autosomal or X-linked
 * markers (not required for Y or mitochondrial loci) */
void set_recode(struct Marker *mark,struct loki *loki,int comp)
{
	int i,j,k,k1,k2,k3,base,psize,n_all,g[4],mask,ltype,diploid=0,comp1,gt_known;
	int *allele_set,*allele_mask,*haplo,*ngens,**req_set[2],nk,sex,*tp,*tp1,cand_gene;
	struct temp_mem *first,*current;
	gtype **gtypes;
	struct Id_Record *id_array,**kids,*kid;
	lk_ulong a;
	
	/* General strategy - we've already performed genotype elimination, so
	 * we can just go through the pedigree in reverse order (i.e., kids before
	 * parents), passing required alleles upwards. */
	
	cand_gene=mark->mterm?1:0;
	ltype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
	if(ltype!=LINK_AUTO && ltype!=LINK_X) ABT_FUNC("Illegal link type\n");
	if(comp<0) { /* Do all components */
		message(DEBUG_MSG,"Set recoding for locus %s\n",mark->name);
		comp=base=0;
		comp1=loki->pedigree->n_comp;
	} else { /* Do a particular component */
		message(DEBUG_MSG,"Set recoding for locus %s component %d\n",mark->name,comp+1);
		base=loki->pedigree->comp_start[comp];
		comp1=comp+1;
	}
	psize=loki->pedigree->ped_size;
	id_array=loki->pedigree->id_array;
	n_all=mark->locus.n_alleles;
	gtypes=mark->gtypes;
	haplo=mark->haplo;
	ngens=mark->ngens;
	if(!(allele_set=malloc(sizeof(int)*n_all*2))) ABT_FUNC(MMsg);
	allele_mask=allele_set+n_all;
	if(!(req_set[0]=malloc(sizeof(void *)*psize*2))) ABT_FUNC(MMsg);
	if(!(first=malloc(sizeof(struct temp_mem)))) ABT_FUNC(MMsg);
	first->next=0;
	first->i=0;
	current=first;
	req_set[1]=req_set[0]+psize;
	for(;comp<comp1;comp++) {
		n_all=mark->n_all1[comp];
		psize=loki->pedigree->comp_size[comp];
		for(i=psize-1;i>=0;i--) {
			j=i+base;
			memset(allele_set,0,sizeof(int)*n_all);
			sex=id_array[j].sex;
			if(ltype==LINK_AUTO || sex==2) diploid=1;
			else if(sex==1) diploid=0;
			else ABT_FUNC("Unknown sex\n");
			/* Check if unordered genotype is known */
			k=ngens[j];
			gt_known=0;
			if(diploid) {
				if(k==1) {
					if(((g[0]=gtypes[j]->mat))&&((g[1]=gtypes[j]->pat))) {
						gt_known=1;
						allele_set[g[0]-1]=1;
						allele_set[g[1]-1]|=2;
					}
				} else if(k==2) {
					for(k3=k2=0;k2<2;k2++) {
						if(!(g[k3++]=gtypes[j][k2].mat)) break;
						if(!(g[k3++]=gtypes[j][k2].pat)) break;
					}
					if(k2==2 && g[0]==g[3] && g[1]==g[2]) {
						gt_known=1;
						allele_set[g[0]-1]=3;
						allele_set[g[1]-1]=3;
					}
				}
			} else {
				if(k==1) {
					if((g[0]=gtypes[j]->mat)) {
						gt_known=1;
						allele_set[g[0]-1]=1;
					}
				}
			}
			if(!gt_known) { /* Genotype unknown, go to kids */
				/* First we sort out which alleles can come from which parent */
				if(k) {
					memset(allele_mask,0,sizeof(int)*n_all);
					for(mask=k1=0;k1<k && mask!=3;k1++) {
						if(!(mask&1)) {
							if((k3=gtypes[j][k1].mat)) allele_mask[k3-1]|=1;
							else mask|=1;
						}
						if(diploid && !(mask&2)) {
							if((k3=gtypes[j][k1].pat)) allele_mask[k3-1]|=2;
							else mask|=2;
						}
					}
					if(mask) for(k1=0;k1<n_all;k1++) allele_mask[k1]|=mask;
				} else {
					k2=diploid?3:1;
					for(k1=0;k1<n_all;k1++) allele_mask[k1]=k2;
				}
				/* First get observed alleles (if any) */
				if((k1=haplo[j])) {
					g[0]=(int)(sqrt(.25+2.0*(double)k1)-.49999);
					g[1]=k1-(g[0]*(g[0]+1)/2);
					if(diploid) {
						if(g[0]) allele_set[g[0]-1]=3;
						if(g[1]) allele_set[g[1]-1]=3;
					} else {
						if(g[0]) allele_set[g[0]-1]=1;
						if(g[1] && g[1]!=g[0]) ABT_FUNC("Heterozygous male for X-linked marker - this should not have got this far\n");
					}
				}
				/* Now do kids */
				nk=id_array[j].nkids;
				kids=id_array[j].kids;
				for(k1=0;k1<nk;k1++) {
					kid=kids[k1];
					if(!diploid && kid->sex==1) continue; /* No transmission from father to male offspring */
					k2=kids[k1]->idx;
					tp=req_set[2-sex][k2];
					assert(tp);
					while((k3=(*tp++))) allele_set[k3-1]=3;
				}
				/* Now get parental origin */
				for(k1=0;k1<n_all;k1++) allele_set[k1]&=allele_mask[k1];
			}
			/* Store required sets */
			for(k1=k2=1,k=0;k<n_all;k++) {
				if(cand_gene) allele_set[k]=3;
				k3=allele_set[k];
				if(k3&1) k1++;
				if(k3&2) k2++;
			}
			if(current->i+k1+k2>TMEM_SIZE) {
				/* This is unlikely to happen, but still... */
				if(k1+k2>TMEM_SIZE) abt(__FILE__,__LINE__,"%s(): TMEM_SIZE too small.  Set to at least %d\n",__func__,k1+k2);
				if(!(current->next=malloc(sizeof(struct temp_mem)))) ABT_FUNC(MMsg);
				current=current->next;
				current->i=0;
				current->next=0;
			}
			k=current->i;
			req_set[0][j]=tp=current->p+k;
			req_set[1][j]=tp1=current->p+k+k1;
			current->i+=k1+k2;
			k=diploid?3:1;
			mask=(gt_known?0:k);
			if(k1==n_all) mask&=~1;
			if(k2==n_all) mask&=~2;
			for(k2=k1=k=0;k<n_all;k++) {
				k3=allele_set[k];
				if(k3&1) tp[k1++]=k+1;
				if(k3&2) tp1[k2++]=k+1;
			}
			tp[k1]=0;
			tp1[k2]=0;
			/* Compat routine to generate lumped allele structures for oldstyle peel routines. To be replaced */
			if(n_all>(int)(sizeof(lk_ulong)<<3)) ABT_FUNC("Too many alleles for current version\n");
			for(k=0;k<2;k++) {
				a=0L;
				k2=0;
				if(mask&(k+1)) {
					for(k1=0;k1<n_all;k1++) if(!(allele_set[k1]&(k+1))) {
						if(allele_mask[k1]&(k+1)) {
							a|=(1L<<k1);
							k2++;
						}
					}
				}
				if(a) a|=(1<<(n_all-1));
				mark->req_set[k][j]=k2>1?a:0L;
			}
		}
		base+=psize;
	}
	/* Clean up */
	while(first) {
		current=first->next;
		free(first);
		first=current;
	}
	free(req_set[0]);
	free(allele_set);
}
