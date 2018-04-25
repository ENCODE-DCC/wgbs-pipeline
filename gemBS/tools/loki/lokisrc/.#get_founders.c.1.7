/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - Centre National de Genotypage, Evry            *
 *                                                                          *
 *                       July 2002                                          *
 *                                                                          *
 * get_founders.c:                                                          *
 *                                                                          *
 * Make compact list of founders for each individual (used for assessing    *
 * who is related to who)                                                   *
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
#include "loki_utils.h"
#include "loki_ibd.h"

void get_founders(unsigned long **found,int **inb,int **n_l,int **n_p,const struct loki *loki)
{
  int *n_longs,i,j,k,k1,ids,idd,nf,nf1,longbit,nlong,cs,comp,*n_pairs,*inbr;
  unsigned long *founders,*tp,*tp1,*tp2;
	
  message(INFO_MSG,"Setting up founder lists\n");
  if(!(n_longs=malloc(sizeof(int)*(2*loki->pedigree->n_comp+loki->pedigree->ped_size)))) exit(EXIT_FAILURE);
  n_pairs=n_longs+loki->pedigree->n_comp;
  inbr=n_pairs+loki->pedigree->n_comp;
  longbit=sizeof(long)<<3;
  for(i=nf1=comp=0;comp<loki->pedigree->n_comp;comp++) {
    cs=loki->pedigree->comp_size[comp];
    for(nf=j=0;j<cs;j++,i++) {
      if(!(loki->pedigree->id_array[i].sire)) nf++;
      if(!(loki->pedigree->id_array[i].dam)) nf++;
    }
    n_longs[comp]=nf/longbit;
    if(n_longs[comp]*longbit<nf) n_longs[comp]++;
    nf1+=n_longs[comp]*cs;
  }
  if(!(founders=malloc(sizeof(long)*nf1))) exit(EXIT_FAILURE);
  tp=founders;
  for(i=nf1=comp=0;comp<loki->pedigree->n_comp;comp++) {
    n_pairs[comp]=0;
    cs=loki->pedigree->comp_size[comp];
    nlong=n_longs[comp];
    tp2=tp;
    for(nf=j=0;j<cs;j++,tp+=nlong) {
      inbr[i+j]=0;
      ids=loki->pedigree->id_array[i+j].sire;
      if(!ids) {
	for(k=0;k<nlong;k++) tp[k]=0;
	k=nf/longbit;
	k1=nf%longbit;
	tp[k]=1L<<k1;
	nf++;
      } else {
	tp1=tp2+(ids-1-i)*nlong;
	for(k=0;k<nlong;k++) tp[k]=tp1[k];
      }
      idd=loki->pedigree->id_array[i+j].dam;
      if(!idd) {
	k=nf/longbit;
	k1=nf%longbit;
	tp[k]|=1L<<k1;
	nf++;
      } else {
	tp1=tp2+(idd-1-i)*nlong;
	for(k=0;k<nlong;k++) {
	  if(tp[k]&tp1[k]) inbr[i+j]=1;
	  tp[k]|=tp1[k];
	}
      }
      tp1=tp2;
      for(k=0;k<j;k++,tp1+=nlong) {
	for(k1=0;k1<nlong;k1++) if(tp[k1]&tp1[k1]) break;
	if(k1<nlong) if((inbr[i+k]||inbr[i+j])||((ids!=i+k+1)&&(idd!=i+k+1))) n_pairs[comp]++;
      }
      if(inbr[i+j]) n_pairs[comp]++;
    }
    i+=cs;
  }
  *found=founders;
  *n_l=n_longs;
  *n_p=n_pairs;
  *inb=inbr;
}
