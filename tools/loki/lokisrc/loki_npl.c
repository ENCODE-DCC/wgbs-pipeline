/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - CNG, Evry                                *
 *                                                                          *
 *                        July 2002                                         *
 *                                                                          *
 * loki_npl.c:                                                              *
 *                                                                          *
 * Routines for calculating NPL scores & associated stats                   *
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
#include "seg_pen.h"
#include "loki_ibd.h"
#include "loki_npl.h"

static int *naff,**affs;

void SetupNPL(struct loki *loki)
{
	int i,i1,i2,j,k,k1,k2,comp,*tp,n_long,cs,*n_longs;
	int ids,idd,ids1,idd1,sh,Spairs;
	unsigned long *tpl,*founders;
	struct Id_Record *id_array;
	
	get_founder_params(&founders,&n_longs,0,0,loki);
	id_array=loki->pedigree->id_array;
	if(!(affs=malloc(sizeof(void *)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
	if(!(naff=malloc(sizeof(int)*loki->pedigree->n_comp))) ABT_FUNC(MMsg);
	/* Count the number of affecteds in each component */
	for(k=i=comp=0;comp<loki->pedigree->n_comp;comp++) {
		naff[comp]=0;
		affs[comp]=0;
		for(j=0;j<loki->pedigree->comp_size[comp];j++) if(id_array[i+j].affected==2) naff[comp]++;
		k+=naff[comp];
	}
	/* List the affecteds in each component and calculate the maximum Spairs value */
	if(k) {
		if(!(tp=malloc(sizeof(int)*k))) ABT_FUNC(MMsg);
		tpl=founders;
		for(i=comp=0;comp<loki->pedigree->n_comp;comp++) {
			n_long=n_longs[comp];
			cs=loki->pedigree->comp_size[comp];
			if(naff[comp]) {
				affs[comp]=tp;
				/* Make list of affecteds in this component */
				for(k=j=0;j<cs;j++) if(id_array[i+j].affected==2) {
					tp[k++]=i+j;
				}
				Spairs=0;
				/* Check which affected pairs are related */
				for(j=0;j<k;j++) {
					i1=tp[j];
					ids=id_array[i1].sire;
					idd=id_array[i1].dam;
					for(k1=0;k1<j;k1++) {
						i2=tp[k1];
						for(k2=0;k2<n_long;k2++) if(tpl[i1*n_long+k2]&tpl[i2*n_long+k2]) break;
						if(k2<n_long) {
							/* Pair is related.  Check if possible to share 2 genes IBD */
							ids1=id_array[i2].sire;
							idd1=id_array[i2].dam;
							sh=1;
							if(ids&&ids1&&idd&&idd1) {
								for(k2=0;k2<n_long;k2++) if(tpl[(ids-1)*n_long+k2]&tpl[(ids1-1)*n_long+k2]) break;
								if(k2<n_long) for(k2=0;k2<n_long;k2++) if(tpl[(idd-1)*n_long+k2]&tpl[(idd1-1)*n_long+k2]) break;
								if(k2<n_long) sh=2;
								else {
									for(k2=0;k2<n_long;k2++) if(tpl[(ids-1)*n_long+k2]&tpl[(idd1-1)*n_long+k2]) break;
									if(k2<n_long) for(k2=0;k2<n_long;k2++) if(tpl[(idd-1)*n_long+k2]&tpl[(ids1-1)*n_long+k2]) break;
									if(k2<n_long) sh=2;
								}
								Spairs+=sh;
							}
						}
					}
				}
				printf("Comp %d, Sp %d\n",comp,Spairs);
				tp+=k;
			} 
			tpl+=n_long*cs;
		}
	}
}

