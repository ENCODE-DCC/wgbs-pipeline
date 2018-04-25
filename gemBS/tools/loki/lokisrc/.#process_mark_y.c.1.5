/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                            May 2003                                      *
 *                                                                          *
 * process_mark_y.c:                                                        *
 *                                                                          *
 * Check inheritance of y, w or mitochondrial linked markers                *
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

#include "utils.h"
#include "libhdr.h"
#include "version.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "ped_utils.h"
#include "gen_elim.h"

/* Assign initial genotype lists based on Y-linked or mitochondrial microsat marker data */
static gen_elim_err *assign_gtypes_y(int ix,int *hap,gtype **gts,int *ngens,int sex,gen_elim_err *gerr,int ltype,struct loki *loki)
{
  int k,a1,a2,err=0;
  gtype *gt=0;
	
  if((k=hap[ix])) {
    a1=(int)(sqrt(.25+2.0*(double)k)-.49999);
    a2=k-(a1*(a1+1)/2);
    if(ltype==LINK_Y) {
      if(sex==2) {
	gerr=add_gen_elim_err(gerr,0,loki->pedigree->id_array+ix,GEN_ELIM_Y_OBS_FEMALE);
	err=1;
      } else {
	if(a2 && a1!=a2) {
	  gerr=add_gen_elim_err(gerr,0,loki->pedigree->id_array+ix,GEN_ELIM_Y_HET_MALE);
	  err=1;
	} else {
	  if(!(gt=malloc(sizeof(gtype)))) ABT_FUNC(MMsg);
	  gt[0].pat=a1;
	}
      }
    } else { /* Mitochondrial or W markers act the same way */
      if(sex==1) {
	gerr=add_gen_elim_err(gerr,0,loki->pedigree->id_array+ix,GEN_ELIM_MIT_OBS_MALE);
	err=1;
      } else {
	if(a2 && a1!=a2) {
	  gerr=add_gen_elim_err(gerr,0,loki->pedigree->id_array+ix,GEN_ELIM_MIT_HET_FEMALE);
	  err=1;
	} else {
	  if(!(gt=malloc(sizeof(gtype)))) ABT_FUNC(MMsg);
	  gt[0].mat=a1;
	}
      }
    }
    if(err) {
      ngens[ix]=0;
      gts[ix]=0;
    } else {
      gts[ix]=gt;
      ngens[ix]=1;
    }
  } else {
    ngens[ix]=0;
    gts[ix]=0;
  }
  return gerr;
}

gen_elim_err *process_mark_y(struct Marker *mark,struct loki *loki,int ltype,int comp)
{
  int i,j,k,k1,a,par,*prune,*haplo,*ngens,ped_start,ped_stop,sx,*line,line_sz=8;
  struct Id_Record *id;
  gtype **gtypes;
  gen_elim_err *gerr=0;
  nuc_fam **fam_ptr;
	
  prune=mark->locus.pruned_flag;
  haplo=mark->haplo;
  gtypes=mark->gtypes;
  ngens=mark->ngens;
  id=loki->pedigree->id_array;
  if(comp<0) { /* Do all components */
    ped_start=0;
    ped_stop=loki->pedigree->ped_size;
  } else { /* Do a particular component */
    ped_start=loki->pedigree->comp_start[comp];
    ped_stop=ped_start+loki->pedigree->comp_size[comp];
  }
  if(!(line=malloc(sizeof(int)*line_sz))) ABT_FUNC(MMsg);
  sx=(ltype==LINK_Y?1:2);
  message(DEBUG_MSG,"Pass 1\n");
  for(i=ped_start,j=0;i<ped_stop;i++) {
    ngens[i]=0;
    gtypes[i]=0;
    id[i].flag=0;
    if(!prune[i]) {
      gerr=assign_gtypes_y(i,haplo,gtypes,ngens,id[i].sex,gerr,ltype,loki);
      if(id[i].sex!=sx) continue;
      par=(sx==1?id[i].sire:id[i].dam);
      if(par && !prune[par-1]) {
	k=id[par-1].flag;
#ifdef DEBUG
	if(!k || k>j) ABT_FUNC("Internal error\n");
#endif
	id[i].flag=k;
	if(ngens[i]) line[k-1]++;
      } else {
	id[i].flag=++j;
	if(j==line_sz) {
	  line_sz<<=1;
	  if(!(line=realloc(line,sizeof(int)*line_sz))) ABT_FUNC(MMsg);
	}
	line[j-1]=ngens[i]?1:0;
      }
    }
  }
  if(j) {
    message(DEBUG_MSG,"Pass 2\n");
    if(!(fam_ptr=malloc(sizeof(void *)*j))) ABT_FUNC(MMsg);
    for(i=0;i<j;i++) {
      if(line[i]) {
	if(!(fam_ptr[i]=malloc(sizeof(nuc_fam)))) ABT_FUNC(MMsg);
	if(!(fam_ptr[i]->gtypes=malloc(sizeof(struct f_gtype)))) ABT_FUNC(MMsg);
	if(!(fam_ptr[i]->kids=malloc(sizeof(void *)*(line[i]+1)))) ABT_FUNC(MMsg);
	fam_ptr[i]->n_err=fam_ptr[i]->flag=0;
	line[i]=0;
      } else fam_ptr[i]=0;
    }
    for(i=ped_start;i<ped_stop;i++) if((k=id[i].flag)) {
      k--;
      if(!fam_ptr[k]) continue;
      if(!fam_ptr[k]->flag) {
	if(ltype==LINK_Y) fam_ptr[k]->father=id+i;
	else fam_ptr[k]->mother=id+i;
	fam_ptr[k]->flag=1;
      }
      if(ngens[i]) {
	k1=line[k]++;
	fam_ptr[k]->kids[k1]=id+i;
	fam_ptr[k]->kids[k1+1]=0;
	if(fam_ptr[k]->flag<2) {
	  if(ltype==LINK_Y) {
	    a=gtypes[i]->pat;
	    fam_ptr[k]->gtypes->par[1].pat=a;
	  } else {
	    a=gtypes[i]->mat;
	    fam_ptr[k]->gtypes->par[0].mat=a;
	  }
	  fam_ptr[k]->flag=2;
	} else {
	  if(ltype==LINK_Y) {
	    a=gtypes[i]->pat;
	    if(a!=fam_ptr[k]->gtypes->par[1].pat) fam_ptr[k]->n_err++;
	  } else {
	    a=gtypes[i]->mat;
	    if(a!=fam_ptr[k]->gtypes->par[0].mat) fam_ptr[k]->n_err++;
	  }
	}
      }
    }
    for(i=0;i<j;i++) if(fam_ptr[i]) {
      free(fam_ptr[i]->gtypes);
      if(fam_ptr[i]->n_err) {
	gerr=add_gen_elim_err(gerr,fam_ptr[i],ltype==LINK_Y?fam_ptr[i]->father:fam_ptr[i]->mother,GEN_ELIM_PASS2_YM);
      } else free(fam_ptr[i]);
    }
    free(fam_ptr);
  }
  free(line);
  return invert_err_list(gerr);
}
