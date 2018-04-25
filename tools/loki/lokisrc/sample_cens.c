/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                      Simon Heath - MSKCC                                 *
 *                                                                          *
 *                          August 2000                                     *
 *                                                                          *
 * sample_cens.c:                                                           *
 *                                                                          *
 * Routines for sampling censored data                                      *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "ranlib.h"
#include "sample_cens.h"

/* Sample truncated (censored) data points */
void Sample_Censored(struct loki *loki)
{
	int i,j,type,idx,er;
	double sd,y,z;
	struct Model *mod;
	struct Id_Record *id_array;

	sd=sqrt(loki->models->residual_var[0]);
	mod=loki->models->models;
	type=mod->var.type;
	idx=mod->var.var_index;
	id_array=loki->pedigree->id_array;
	if(type&ST_CONSTANT)	{
		for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0] && (id_array[i].data[idx].flag&2)) {
			y=id_array[i].res[0][0]-id_array[i].pseudo_qt[0][0];
			z=sd*trunc_normal(y/sd,0.0,2,&er); /* Sample from normal left truncated at y */
			if(!er) {
				id_array[i].pseudo_qt[0][0]=z-y;
				id_array[i].res[0][0]=z;
			}
		}
	} else for(i=0;i<loki->pedigree->ped_size;i++) {
		for(j=0;j<id_array[i].n_rec;j++) if(id_array[i].res[j] && (id_array[i].data1[j][idx].flag&2)) {
			y=id_array[i].res[0][j]-id_array[i].pseudo_qt[0][j];
			z=sd*trunc_normal(y/sd,0.0,2,&er);
			if(!er) {
				id_array[i].pseudo_qt[0][j]=z-y;
				id_array[i].res[0][j]=z;
			}
		}
	}
}

