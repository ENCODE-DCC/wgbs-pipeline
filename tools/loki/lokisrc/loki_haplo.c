/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - CNG                                   *
 *                                                                          *
 *                          December 2002                                   *
 *                                                                          *
 * loki_haplo.c:                                                            *
 *                                                                          *
 * Routines for haplotype estimation                                        *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
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
#include "loki_peel.h"
#include "version.h"
#include "seg_pen.h"
#include "loki_haplo.h"

void sample_haplo(FILE *fptr,int *naff,int **affs,int ***haplo_store,int si_mode)
{
	int i,j,k,k1,k2,id;
	double z;
	
	for(j=0;j<n_markers;j++) {
		pass_founder_genes(&marker[j].locus);
		for(k2=i=0;i<n_comp;i++) if(naff[i]) { /* Check components with 'affected' individuals */
			z=seg_pen(j,i,&k,1,si_mode);
			if(k) {
				fprintf(stderr,"seg_pen() returned an error (%g)\n",z);
				ABT_FUNC("Aborting\n");
			}
			for(k=0;k<naff[i];k++,k2++) {
				id=affs[i][k];
				for(k1=0;k1<2;k1++) haplo_store[k1][k2][j]=id_array[id].allele[k1];
			}
		}
	}
	for(k2=i=0;i<n_comp;i++) {
		for(k=0;k<naff[i];k++,k2++) {
			id=affs[i][k];
			print_orig_id(fptr,id);
			for(j=0;j<n_markers;j++) {
				printf("
		}
	}
}
