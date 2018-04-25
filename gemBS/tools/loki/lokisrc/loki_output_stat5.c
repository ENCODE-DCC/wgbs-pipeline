/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Paris                             *
 *                                                                          *
 *                          January 2004                                    *
 *                                                                          *
 * loki_output_stat5.c:                                                     *
 *                                                                          *
 * Routines for sample output for EWD's stat5 code                          *
 * Modified by SCH to merge with current loki code - Jan 2004               *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2004                                        *
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
#include "loki_utils.h"
#include "loki_output_stat5.h"

/* June 3, 1998 EWD - output a QTL vector to file */
/* Jan 12, 2004 EWD - update? for new structure */
void OutputQTLvect(FILE *qptr,int lp,struct loki *loki)
{
   int l,n;
	struct Locus *loc;
	
	loc=loki->models->tlocus;
   for(l=0;l<loki->params.n_tloci;l++) if(loc[l].flag) {
      fprintf(qptr,"%d ",lp);
      if(loc[l].flag&TL_LINKED) {
         fprintf(qptr,"%d ",loc[l].link_group+1);
      } else fputs("0 ",qptr);
      for(n=0;n<loki->pedigree->ped_size;n++) fprintf(qptr,"%d ",loc[l].gt[n]);
      fputc('\n',qptr);
   }
}

void OutputQTLHead(FILE *qptr, struct loki *loki)
/* June 3, 1998 EWD - print header information for QTL genotype vector file */
{
   int i;
   fprintf(qptr," %d %d\n",loki->pedigree->ped_size,loki->pedigree->n_comp);
   for(i=0;i<loki->pedigree->ped_size;i++) {
      print_orig_id(qptr,i+1);
      fputc(' ',qptr);
   }
   fputc('\n',qptr);
   for(i=0;i<loki->pedigree->ped_size;i++) {
      fprintf(qptr,"%d ",loki->pedigree->id_array[i].comp);
   }
   fputc('\n',qptr);
}

