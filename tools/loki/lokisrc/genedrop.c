/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * genedrop.c:                                                              *
 *                                                                          *
 * Drop genes down pedigree consistent with seg pattern                     *
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
#include "ranlib.h"
#include "loki_peel.h"
#include "genedrop.h"

void genedrop(int locus,const struct loki *loki)
{
	int i,ped_size,**seg;
	struct Id_Record *id_array;
	
	id_array=loki->pedigree->id_array;
	ped_size=loki->pedigree->ped_size;
	seg=loki->markers->marker[locus].locus.seg;
	for(i=0;i<ped_size;i++) {
		if(id_array[i].sire) seg[X_PAT][i]=ranf()<.5?0:1;
		else seg[X_PAT][i]= -1;
		if(id_array[i].dam) seg[X_MAT][i]=ranf()<.5?0:1;
		else seg[X_MAT][i]= -1;
	}
}

void drop_genotypes(int locus,const struct loki *loki)
{
	int i,j,n_all,ids,idd,ped_size;
	double *freq,*cm,z;
	
	struct Id_Record *id_array;
	
	freq=loki->markers->marker[locus].locus.freq[0];
	n_all=loki->markers->marker[locus].locus.n_alleles;
	if(n_all<1) return;
	id_array=loki->pedigree->id_array;
	ped_size=loki->pedigree->ped_size;
	if(!(cm=malloc(sizeof(double)*n_all))) ABT_FUNC(MMsg);
	cm[0]=freq[0];
	for(i=1;i<n_all;i++) cm[i]=cm[i-1]+freq[i];
	for(i=0;i<ped_size;i++)	{
		if((ids=id_array[i].sire)) id_array[i].allele[X_PAT]=id_array[ids-1].allele[ranf()<.5?0:1];
		else {
			z=ranf();
			for(j=0;j<n_all;j++) if(z<=cm[j]) break;
			id_array[i].allele[X_PAT]=j+1;
		}
		if((idd=id_array[i].dam)) id_array[i].allele[X_MAT]=id_array[idd-1].allele[ranf()<.5?0:1];
		else {
			z=ranf();
			for(j=0;j<n_all;j++) if(z<=cm[j]) break;
			id_array[i].allele[X_MAT]=j+1;
		}
	}
	free(cm);
}

