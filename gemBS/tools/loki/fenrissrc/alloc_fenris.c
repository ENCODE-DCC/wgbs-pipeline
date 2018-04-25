/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * alloc_fenris.c:                                                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris.h"

int n_genetic_groups=1,n_tloci=0;
struct TraitLocus *tlocus;

void AllocFenrisStruct(void)
{
	int i,j,k,n_all,*ff,*tp1;
#ifdef DEBUG
	int k2;
#endif
	double **tpp,*tp;
	signed char *tcp,**tcpp;
	struct Locus *loc;
	
	if(n_markers) {
		for(i=0;i<n_markers;i++) {
			ff=founder_flag[i];
			for(j=0;j<ped_size;j++) {
				if(id_array[j].pruned_flag[i]) ff[j]=2;
				else {
					k=id_array[j].sire;
					if(k) {
						if(id_array[k-1].pruned_flag[i]) ff[j]=1;
						else ff[j]=0;
					} else ff[j]=1;
#ifdef DEBUG
					k=id_array[j].dam;
					if(k) {
						if(id_array[k-1].pruned_flag[i]) k2=1;
						else k2=0;
					} else k2=1;
					if(k2!=ff[j]) ABT_FUNC("Bad pruning - half pruned family\n");
#endif
				}
			}
		}
		for(j=i=0;i<n_markers;i++) j+=marker[i].locus.n_alleles;
		if(j) {
			if(!(tpp=malloc(sizeof(void *)*n_genetic_groups*n_markers))) ABT_FUNC(MMsg);
			if(!(tcpp=malloc(sizeof(void *)*n_genetic_groups*n_markers))) ABT_FUNC(MMsg);
			if(!(tcp=malloc((size_t)j*n_genetic_groups*n_markers))) ABT_FUNC(MMsg);
			if(!(tp=malloc(sizeof(double)*j*n_genetic_groups*n_markers))) ABT_FUNC(MMsg);
			if(!(tp1=malloc(sizeof(int)*n_genetic_groups*n_markers))) ABT_FUNC(MMsg);
			RemBlock=AddRemem(tpp,RemBlock);
			RemBlock=AddRemem(tcpp,RemBlock);
			RemBlock=AddRemem(tcp,RemBlock);
			RemBlock=AddRemem(tp,RemBlock);
			RemBlock=AddRemem(tp1,RemBlock);
			for(i=0;i<n_markers;i++) {
				loc=&marker[i].locus;
				n_all=loc->n_alleles;
				loc->aff_freq=loc->diff_freq=0;
				if(!n_all) {
					loc->freq=0;
					marker[i].freq_set=0;
					marker[i].count_flag=0;
					continue;
				}
				loc->freq=tpp;
				marker[i].freq_set=tcpp;
				marker[i].count_flag=tp1;
				tcpp+=n_genetic_groups;
				tpp+=n_genetic_groups;
				tp1+=n_genetic_groups;
				for(j=0;j<n_genetic_groups;j++) {
					marker[i].count_flag[j]=0;
					loc->freq[j]=tp;
					marker[i].freq_set[j]=tcp;
					tp+=n_all;
					tcp+=n_all;
					for(k=0;k<n_all;k++) {
						loc->freq[j][k]=0.0;
						marker[i].freq_set[j][k]=0;
					}
				}
			}
		} else for(i=0;i<n_markers;i++) {
			marker[i].locus.freq=0;
			marker[i].freq_set=0;
		}
		for(i=0;i<n_markers;i++) {
			marker[i].pos_set=0;
			marker[i].locus.pos[0]=marker[i].locus.pos[1]=0.0;
		}
	}
}

