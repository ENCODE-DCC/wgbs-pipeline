/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                       CNG  - Evry                                        *
 *                                                                          *
 *                       July 2002                                          *
 *                                                                          *
 * update_segs.c:                                                           *
 *                                                                          *
 * Interface routines between the L and M samplers                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "seg_pen.h"
#include "gen_pen.h"
#include "update_segs.h"

/* 
 * fg controls whether we are sampling, and whether we will update frequency estimates.
 * If bit 1 is set, we sample marker genotypes.  If bit 2 is set, we sample QTL genotypes.
 * If bit 4 is set, we sample allele frequencies (markers + QTL).  If bit 8 is set
 * we update allele frequencies in affecteds (markers only).
 *
 * fg1 controls whether we need to update founder genes for markers
 * and QTLs.  If bit 1 is set, founder genes for markers need updating,
 * and if bit 2 is set, founder genes for QTLs need updating
 */
void update_seg_probs(int fg,int fg1,struct loki *loki)
{
	int i,k,comp;
	struct Locus *loc;
	double z;
	
	for(k=0;k<loki->markers->n_markers;k++) {
		loc=&loki->markers->marker[k].locus;
		if(fg&12) seg_init_freq(loc,loki);
		if(fg1&1) pass_founder_genes(loc,loki);
		for(comp=0;comp<loki->pedigree->n_comp;comp++) {
			loc->lk_store[comp]=seg_pen(loc,comp,&i,fg&~2,loki);
			if(i) {
				(void)fprintf(stderr,"seg_pen returned error code %d for marker %s",(int)loc->lk_store[comp],loki->markers->marker[k].name);
				(void)fprintf(stderr," (comp=%d, fg=%d, fg1=%d)\n",comp,fg,fg1);
				ABT_FUNC("Illegal segregation pattern\n");
			}
#ifndef NDEBUG
			z=loc->lk_store[comp];
			if(isinf(z) || isnan(z)) {
				fprintf(stderr,"seg_pen returned illegal value %g (loc=%s, comp=%d, fg=%d, fg1=%d)\n",z,loc->name,comp,fg,fg1);
			}
#endif
		}
		if(fg&4) {
			seg_sample_freq(loc,loki);
		}
		if(fg&8) seg_update_aff_freq(loc,loki);
		if(fg) loc->flag|=LOCUS_SAMPLED;
	}
	for(k=0;k<loki->params.n_tloci;k++) {
		loc=loki->models->tlocus+k;
		if(loc->flag) {
			if(fg&4) seg_init_freq(loc,loki);
			if(fg1&2) pass_founder_genes(loc,loki);
			for(comp=0;comp<loki->pedigree->n_comp;comp++) {
				loc->lk_store[comp]=gen_pen(loc,comp,&i,fg&~1,loki);
				if(i) {
					(void)fprintf(stderr,"%d %d %d %d\n",k,comp,i,(int)loc->lk_store[comp]);
					ABT_FUNC("Illegal segregation pattern\n");
				}
			}
			if(fg&~1) loc->flag|=LOCUS_SAMPLED;
			if(fg&4) seg_sample_freq(loc,loki);
		}
	}
}

void reprune_segs(struct loki *loki)
{
	int i,j,k,**seg,*ff;
	
	for(j=0;j<loki->markers->n_markers;j++) {
		seg=loki->markers->marker[j].locus.seg;
		ff=loki->markers->marker[j].locus.founder_flag;
		for(i=0;i<loki->pedigree->ped_size;i++) if(loki->pedigree->id_array[i].sire) {
			k=ff[i];
			if(k==2) seg[0][i]=seg[1][i]=-2;
			else if(k==1) seg[0][i]=seg[1][i]=-1;
		}
	}
}

