/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*                     Simon Heath - CNG, Evry                              *
*                                                                          *
*                         February 2003                                    *
*                                                                          *
* recode_alleles.c:                                                        *
*                                                                          *
* Copyright (C) Simon C. Heath 2003                                        *
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
#include <math.h>
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "gen_elim.h"

void recode_alleles(int i,struct loki *loki) 
{
	int j,k,k1,a1,a2,id;
	int comp,n_comp,cs,n_all,n_all1,*tab1,*trans;
	struct Marker *mark;
	struct Locus *loc;
	struct Id_Record *id_array;
	
	/* Count how many alleles are found in each component */
	mark=loki->markers->marker+i;
	loc=&mark->locus;
	id_array=loki->pedigree->id_array;
	message(DEBUG_MSG,"Recoding alleles for marker %s\n",mark->name);
	n_all=loc->n_alleles;
	if(n_all<2) return;
	tab1=lk_malloc(sizeof(int)*n_all);
	n_comp=loki->pedigree->n_comp;
	for(id=comp=0;comp<n_comp;comp++) {
		trans=mark->allele_trans[comp];
		for(j=0;j<n_all-1;j++) tab1[j]=0;
		tab1[j]=1;
		cs=loki->pedigree->comp_size[comp];
		for(j=0;j<cs;j++) {
			for(k=0;k<id_array[id+j].n_gt_sets;k++) {
				if((k1=mark->orig_gt[id_array[id+j].gt_idx[k]])) {
					a1=(int)(sqrt(.25+2.0*(double)k1)-.49999);
					a2=k1-(a1*(a1+1)/2);
					assert(a1);
					tab1[a1-1]=1;
					if(a2) tab1[a2-1]=1;
				}
			}
		}
		for(j=0;j<n_all;j++) trans[j]=-1;
		for(n_all1=j=0;j<n_all;j++) {
			if(tab1[j]) {
				trans[n_all1++]=j;
				tab1[j]=n_all1;
			}
		}
		if(n_all1<1) { 
			/* No observed alleles in this component */
			mark->n_all1[comp]=0;
		} else {
			if(n_all1<n_all) n_all1++;
			mark->n_all1[comp]=n_all1;
			if(n_all1<n_all) { 
				/* Recode alleles */
				message(DEBUG_MSG,"No. alleles reduced from %d to %d in component %d\n",n_all,n_all1,comp+1);
				for(j=0;j<cs;j++) {
					for(k=0;k<id_array[id+j].n_gt_sets;k++) {
						if((k1=mark->orig_gt[id_array[id+j].gt_idx[k]])) {
							a1=(int)(sqrt(.25+2.0*(double)k1)-.49999);
							a2=k1-(a1*(a1+1)/2);
							a1=tab1[a1-1];
							if(a2) a2=tab1[a2-1];
							assert(a1>=a2);
							k1=a1*(a1+1)/2+a2;
							mark->orig_gt[id_array[id+j].gt_idx[k]]=k1;
							if(mark->haplo[id+j]>0) mark->haplo[id+j]=k1;
						}
					}
				}
			}
		}
		id+=cs;
	}
	free(tab1);
}

