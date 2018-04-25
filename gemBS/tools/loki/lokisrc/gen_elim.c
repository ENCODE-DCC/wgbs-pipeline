/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * gen_elim.c:                                                              *
 *                                                                          *
 * Perform error checking, genotype elimination and set recoding            *  
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
#include "oldpeel_compat.h"
#include "gen_elim.h"

static void free_gtlists(struct loki *loki)
{
  int i,nf;
  nuc_fam *fam;
	
  fam=loki->family->families;
  nf=loki->family->cp_start[loki->pedigree->n_comp];
  for(i=0;i<nf;i++) {
    if(fam[i].gtypes) {
      free(fam[i].gtypes);
      fam[i].gtypes=0;
    }
  }
}

static void free_gtypes(struct Marker *mark,struct loki *loki)
{
  int i,pedsize;
  gtype **gtypes;
	
  gtypes=mark->gtypes;
  if(gtypes) {
    pedsize=loki->pedigree->ped_size;	
    for(i=0;i<pedsize;i++) if(gtypes[i]) {
      free(gtypes[i]);
      gtypes[i]=0;
    }
  }
}

gen_elim_err *gen_elim_marker(int ix,struct loki *loki) 
{
  struct Marker *mark;
  int linktype,i,comp;
  gen_elim_err *gerr=0;
  static char *marker_types[]={"","microsatellite","snp","super-locus"};
  
  mark=loki->markers->marker+ix;
  message(DEBUG_MSG,"Genotype elimination for %s %s\n",marker_types[mark->marker_type],mark->name);
  if(mark->n_children) {
	  for(comp=0;comp<loki->pedigree->n_comp;comp++) {
		  for(i=0;i<mark->n_children;i++) mark->child_flag[i]=0;
		  gerr=process_mark(mark,loki,comp);
	  }
  } else gerr=process_mark(mark,loki,-1);
  free_gtlists(loki);
  if(gerr) message(WARN_MSG,"Mendelian errors detected for marker %s\n",mark->name);
  else {
	  linktype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
	  if(linktype==LINK_AUTO || linktype==LINK_X) {
      alloc_reqset(mark,loki->pedigree->ped_size);
      set_recode(mark,loki,-1);
      setup_allset(mark,loki);
    }
    /* If there are errors we don't clean up gtypes yet, because we use them for error reporting */
    free_gtypes(mark,loki);
  }
  return gerr;
}

int Gen_Elim(struct loki *loki) 
{
  int i,n_markers,err=0;
  gen_elim_err **errors;
	
  n_markers=loki->markers->n_markers;
  if(n_markers) {
    message(INFO_MSG,"Genotype elimination\n");
    if(!(errors=malloc(sizeof(void *)*n_markers))) ABT_FUNC(MMsg);
    for(i=0;i<n_markers;i++) recode_alleles(i,loki);
    for(i=0;i<n_markers;i++) if(!loki->markers->marker[i].parent) {
      errors[i]=gen_elim_marker(i,loki);
      if(errors[i]) err++;
    }
    if(err) handle_errors(loki,errors);
    for(i=0;i<n_markers;i++) if(errors[i]) free_gtypes(loki->markers->marker+i,loki);
    free(errors);
  }
  return err;
}
