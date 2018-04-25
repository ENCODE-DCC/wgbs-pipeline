/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                            April 2004                                    *
 *                                                                          *
 * pseudo_chrom.c:                                                          *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <signal.h>
#include <sys/wait.h>
#include <assert.h>

#include "utils.h"
#include "libhdr.h"
#include "ranlib.h"
#include "loki.h"
#include "sparse.h"
#include "min_deg.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "loki_compress.h"
#include "loki_tlmoves.h"
#include "lk_malloc.h"
#include "locus.h"
#include "seg_pen.h"
#include "gen_elim.h"
#include "get_peelseq.h"
#include "pseudo_chrom.h"

#define link my_link

void handle_pseudochrom(struct loki *loki)
{
  int link,i,n,n1,*mk_ix;
  size_t s;
  struct Link *lk,*lk1;
  struct Marker *mk,*mk1;
	
  for(n=n1=link=0;link<loki->markers->n_links;link++) {
    if(loki->markers->linkage[link].type&LINK_MIRRORED) {
      n+=loki->markers->linkage[link].n_markers;
      n1++;
    }
  }
  if(n) { /* Extra markers from mirroring chromosomes */
    i=n+loki->markers->n_markers;
    loki->markers->marker=lk_realloc(loki->markers->marker,(size_t)i*sizeof(struct Marker));
    memset(loki->markers->marker+loki->markers->n_markers,0,n*sizeof(struct Marker));
    mk_ix=lk_malloc(sizeof(int)*n);
    loki->sys.RemBlock=AddRemem(mk_ix,loki->sys.RemBlock);
    /* Now add extra linkage groups */
    n1+=loki->markers->n_links;
    loki->markers->linkage=lk_realloc(loki->markers->linkage,sizeof(struct Link)*n1);
    n1=loki->markers->n_links;
    for(link=0;link<loki->markers->n_links;link++) if(loki->markers->linkage[link].type&LINK_MIRRORED) {
      lk=loki->markers->linkage+link;
      lk1=loki->markers->linkage+(n1++);
      lk1->real_chr=lk;
      s=strlen(lk->name);
      lk1->name=lk_malloc(s+2);
      memcpy(lk1->name,lk->name,s);
      lk1->name[s]='\'';
      lk1->name[s+1]=0;
      lk1->ibd_list=0;
      for(i=0;i<2;i++) {
	lk1->r1[i]=lk->r1[i];
	lk1->r2[i]=lk->r2[i];
	lk1->range_set[i]=lk->range_set[i];
      }
      lk1->ibd_est_type=0;
      lk1->n_markers=lk->n_markers;
      lk1->type=(lk->type&LINK_TYPES_MASK)|LINK_PSEUDO;
      lk1->sample_pos=lk->sample_pos;
      lk1->mk_index=mk_ix;
      mk_ix+=lk1->n_markers;
      for(i=0;i<lk1->n_markers;i++) {
	lk1->mk_index[i]=loki->markers->n_markers;
	mk=loki->markers->marker+lk->mk_index[i];
	mk1=loki->markers->marker+lk1->mk_index[i];
	s=strlen(mk->name);
	mk1->name=lk_malloc(s+2);
	loki->sys.RemBlock=AddRemem(mk1->name,loki->sys.RemBlock);
	memcpy(mk1->name,mk->name,s);
	mk1->name[s]='\'';
	mk1->name[s+1]=0;
	mk1->locus.n_alleles=mk->locus.n_alleles+1; /* Space for unsampled allele */
	alloc_marker(loki->markers->n_markers++,loki);
	mk1->pos_set=mk->pos_set;
	mk1->locus.pos[X_MAT]=mk->locus.pos[X_MAT];
	mk1->locus.pos[X_PAT]=mk->locus.pos[X_PAT];
	mk1->locus.link_group=n1-1;
      }
    }
    loki->markers->n_links=n1;
  }
}

static void sample_pseudo_segs(int link,struct loki *loki)
{
  int i,j,k,n,par,pedsize,**seg,s,newlink;
  struct Locus **list,**list1;
  double *recom[2],x,x1,*p[2],pp[2],z,z1,newpos[2],oldpos[2];
  struct Id_Record *id_array;
  struct Link *lk,*lk_real;
  struct Marker *mk,*mk1;
	
  lk=loki->markers->linkage+link;
  assert(lk->type&LINK_PSEUDO);
  message(INFO_MSG,"Sampling pseudochromosome segs for chromsome '%s'\n",lk->name);
  lk_real=lk->real_chr;
  id_array=loki->pedigree->id_array;
  pedsize=loki->pedigree->ped_size;
  for(i=0;i<lk->n_markers;i++) {
    mk=loki->markers->marker+lk->mk_index[i];
    mk1=loki->markers->marker+lk_real->mk_index[i];
    memcpy(mk->locus.founder_flag,mk1->locus.founder_flag,sizeof(int)*pedsize);
    memcpy(mk->locus.pruned_flag,mk1->locus.pruned_flag,sizeof(int)*pedsize);
    init_marker_segs(mk,loki);
  }
  list=get_sorted_locuslist(link,&n,0);
  /* Set up recombination frequencies */
  if(!n) return;
  recom[0]=lk_malloc(sizeof(double)*2*n);
  recom[1]=recom[0]+n;
  for(par=0;par<2;par++) {
    x=list[0]->pos[par];
    for(j=1;j<n;j++) {
      x1=list[j]->pos[par];
      recom[par][j-1]=.5*(1.0-exp(.02*(x-x1)));
      x=x1;
    }
  }
  /* Do we have any linked trait loci ? */
  for(i=j=0;i<n;i++) if(list[i]->type&ST_TRAITLOCUS) j++;
  /* Yes - sample segs conditional on trait locus pattern */
  if(j) { /* Yes - sample new position for linked trait loci */
    list1=lk_malloc(sizeof(void *)*j);
    for(i=j=0;i<n;i++) if(list[i]->type&ST_TRAITLOCUS) list1[j++]=list[i];
    for(i=0;i<j;i++) {
      newlink=get_tl_position(newpos,loki);
      if(j==1 && (link==newlink || newlink<0)) {
	for(k=0;k<2;k++) list1[0]->pos[k]=newpos[k];
	if(link!=newlink) {
	  list1[0]->link_group=newlink;
	  list1[0]->flag&=~TL_LINKED;
	  list1[0]->flag|=TL_UNLINKED;
	}
      } else {
	z=peel_locus(list1,i,j,0,loki);
	for(k=0;k<2;k++) {
	  oldpos[k]=list1[i]->pos[k];
	  list1[i]->pos[k]=newpos[k];
	}
	list1[i]->link_group=newlink;
	if(newlink<0) {
	  list1[i]->flag&=~TL_LINKED;
	  list1[i]->flag|=TL_UNLINKED;
	}
	z1=peel_locus(list1,i,j,0,loki);
	z=safe_exp(z1-z);
	if(ranf()>=z) {
	  for(k=0;k<2;k++) list1[i]->pos[k]=oldpos[k];
	  list1[i]->link_group=link;
	  list1[i]->flag&=~TL_UNLINKED;
	  list1[i]->flag|=TL_LINKED;
	  (void)peel_locus(list1,i,j,1,loki);
	} else {
	  (void)calc_tl_like(list1[i],1,loki);
	  if(link!=newlink) list1[i--]=list1[--j];
	}
      }
    }
    free(list1);
    list=get_sorted_locuslist(link,&n,0);
  }
  /* Still have linked trait loci? */
  if(j) {/* Yes - sample segs conditional on trait locus pattern */
    p[0]=lk_malloc(sizeof(double)*2*n);
    p[1]=p[0]+n;
    for(i=0;i<pedsize;i++) {
      if(id_array[i].sire) {
	for(par=0;par<2;par++) {
	  if(list[0]->type&ST_TRAITLOCUS) {
	    s=list[0]->seg[par][i];
	    if(s>=0) {
	      p[s][0]=1.0;
	      p[s^1][0]=0.0;
	    } else p[0][0]=p[1][0]=0.5;
	  } else p[0][0]=p[1][0]=0.5;
	  for(j=1;j<n;j++) {
	    if(list[j]->type&ST_TRAITLOCUS) {
	      s=list[j]->seg[par][i];
	      if(s>=0) {
		pp[s]=1.0;
		pp[s^1]=0.0;
	      } else pp[0]=pp[1]=0.5;
	    } else pp[0]=pp[1]=0.5;
	    z=recom[par][j-1];
	    p[0][j]=pp[0]*((1.0-z)*p[0][j-1]+z*p[1][j-1]);
	    p[1][j]=pp[1]*((1.0-z)*p[1][j-1]+z*p[0][j-1]);
	  }
	  j--;
	  z=p[0][j]+p[1][j];
	  s=(ranf()*z<p[0][j])?0:1;
	  list[j--]->seg[par][i]=s;
	  for(;j>=0;j--) {
	    pp[s]=p[s][j]*((1.0-z)*p[s][j+1]+z*p[s^1][j+1]);
	    pp[s^1]=p[s^1][j]*((1.0-z)*p[s^1][j+1]+z*p[s][j+1]);
	    z=pp[0]+pp[1];
	    s=(ranf()*z<pp[0])?0:1;
	    list[j]->seg[par][i]=s;
	  }
	}
      } else for(j=0;j<n;j++) {
	seg=list[j]->seg;
	seg[X_MAT][i]=seg[X_PAT][i]=-1;
      }
    }
    free(p[0]);
  } else { /* No, assign seg pattern at random */
    for(i=0;i<pedsize;i++) {
      if(id_array[i].sire) {
	for(par=0;par<2;par++) {
	  s=list[0]->seg[par][i]=ranf()<.5?1:0;
	  for(j=1;j<n;j++) {
	    if(ranf()<recom[par][j-1]) s^=1;
	    list[j]->seg[par][i]=s;
	  }
	}
      } else for(j=0;j<n;j++) {
	seg=list[j]->seg;
	seg[X_MAT][i]=seg[X_PAT][i]=-1;
      }
    }
  }
  free(recom[0]);
}

static void sample_pseudo_gens(int link,struct loki *loki)
{
  int i,j,k,k1,par,s,ng,nall,**seg,*ff,pedsize;
  double **fq,**fq1,z,zz;
  struct Link *lk,*lk_real;
  struct Id_Record *id_array;
  struct Marker *mk,*mk1;
	
  lk=loki->markers->linkage+link;
  assert(lk->type&LINK_PSEUDO);
  lk_real=lk->real_chr;
  message(INFO_MSG,"Sampling pseudochromosome genotypes for chromosome '%s'\n",lk->name);
  id_array=loki->pedigree->id_array;
  assert(lk->n_markers==lk_real->n_markers);
  ng=loki->pedigree->n_genetic_groups;
  pedsize=loki->pedigree->ped_size;
  for(i=0;i<lk->n_markers;i++) {
    mk=loki->markers->marker+lk->mk_index[i];
    mk1=loki->markers->marker+lk_real->mk_index[i];
    ff=mk->locus.founder_flag;
    /* Copy frequencies from 'real' chromosome */
    fq=mk->locus.freq;
    fq1=mk1->locus.freq;
    nall=mk1->locus.n_alleles;
    for(j=0;j<nall;j++) for(k=0;k<ng;k++) fq[k][j]=fq1[k][j];
    for(k=0;k<ng;k++) fq[k][j]=0.0;
    /* Sample marker genotypes (conditional on previously sampled marker pattern) */
    seg=mk->locus.seg;
    for(j=0;j<pedsize;j++) {
      if(ff[j]) {
	/* Sample founder alleles */
	for(par=0;par<2;par++) {
	  k=id_array[j].group-1;
	  assert(k>=0);
	  zz=fq[k][nall-1];
	  z=ranf()*(1.0-zz);
	  for(k1=0;k1<nall-1;k1++) {
	    z-=fq[k][k1];
	    if(z<=0.0) break;
	  }
	  assert(k1<nall-1);
	  id_array[j].allele[par]=k1+1;
	}
      } else {
	s=seg[X_MAT][j];
	assert(s>=0 && s<2);
	id_array[j].allele[X_MAT]=id_array[id_array[j].dam-1].allele[s];
	s=seg[X_PAT][j];
	assert(s>=0 && s<2);
	id_array[j].allele[X_PAT]=id_array[id_array[j].sire-1].allele[s];
      }
      if(mk1->haplo[j]) {
	k=id_array[j].allele[X_MAT];
	k1=id_array[j].allele[X_PAT];
	if(k>k1) k=k*(k+1)/2+k1;
	else k=k1*(k1+1)/2+k;
	mk->haplo[j]=k;
      } else mk->haplo[j]=0;
    }
  }
}

void sample_pseudo_chromosomes(struct loki *loki)
{
  int link,nm,i,j,ltype,comp;
  struct Link *lk;
  gen_elim_err *err;
	
  for(link=0;link<loki->markers->n_links;link++) {
    if(loki->markers->linkage[link].type&LINK_PSEUDO) {
      lk=loki->markers->linkage+link;
      nm=lk->n_markers;
      if(nm) {
	ltype=lk->type&LINK_TYPES_MASK;
	sample_pseudo_segs(link,loki);
	sample_pseudo_gens(link,loki);
	message(INFO_MSG,"Genotype elimination & peeling sequence determination for pseudochromosome '%s'\n",lk->name);
	for(i=0;i<nm;i++) {
	  j=lk->mk_index[i];
	  free_marker_data(loki->markers->marker+j);
	  recode_alleles(j,loki); 
	  err=gen_elim_marker(j,loki);
	  assert(!err);
	  if(loki->peel->peelseq_head[j]) {
	    for(comp=0;comp<loki->pedigree->n_comp;comp++) free_peelseq(loki->peel->peelseq_head[j]+comp);
	    free(loki->peel->peelseq_head[j]);
	  }
	  loki->peel->peelseq_head[j]=get_peelseq(&loki->markers->marker[j].locus,loki,ltype);
	  loki->markers->marker[j].locus.flag|=LOCUS_SAMPLED;
	}
	get_peelseq(0,0,0);
	min_deg(0,0,0,0,0);
      }
    }
  }
  setup_obslist(loki);
}
