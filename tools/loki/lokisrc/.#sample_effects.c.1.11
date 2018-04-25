/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                   Simon Heath - MSKCC                                    *
 *                                                                          *
 *                       August 2000                                        *
 *                                                                          *
 * sample_effects.c:                                                        *
 *                                                                          *
 * Sampling routine for all effects                                         *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include <errno.h>

#include "ranlib.h"
#include "utils.h"
#include "sparse.h"
#include "loki.h"
#include "loki_utils.h"
#include "sample_rand.h"
#include "sample_effects.h"
#include "min_deg.h"

static int n_var,n_lev,cov_start,mk_start,rand_start,poly_start;
static int bsize=1024,*order,XX_size,entry_size,*zero,init_flag;

static double *B;
static struct Diag *XX;
static struct Off *freelist;
struct entry {
  double val;
  int pos;
};
static double **full_xx;
static struct entry *entry;
static int full_store=60;
static struct loki *loki;

static struct Off *alloc_new_nodes(int n)
{
  struct Off *p;
  int i;
	
#ifdef DEBUG
  if(!n) ABT_FUNC("Internal error - called with zero argument\n");
#endif
  if(!(p=malloc(sizeof(struct Off)*n))) ABT_FUNC(MMsg);
  loki->sys.RemBlock=AddRemem(p,loki->sys.RemBlock);
  /* Link blocks together */
  for(i=0;i<n-1;i++) p[i].Next=p+i+1;
  p[n-1].Next=0;
  return p;
}

static struct Off *get_node(void)
{
  struct Off *p;
	
  if(!freelist) freelist=alloc_new_nodes(bsize);
  p=freelist;
  freelist=p->Next;
  return p;
}

static void add_XX_entry1(int x,int y,double z)
{
  struct Off *p,**pp;
	
  pp=&XX[x].First;
  p=*pp;
  while(p) {
    if(p->col>=y) break;
    pp=&p->Next;
    p=p->Next;
  }
  if(p && p->col==y) p->val+=z;
  else {
    p=get_node();
    p->col=y;
    p->val=z;
    p->Next=*pp;
    *pp=p;
    XX[x].count++;
  }
}

static void add_XX_entry(int x,int y,double z)
{
  struct Off *p,*p2,**pp;
	
  pp=&XX[x].First;
  p=*pp;
  while(p) {
    if(p->col>=y) break;
    pp=&p->Next;
    p=*pp;
  }
  if(p && p->col==y) p->val+=z;
  else {
    if(!freelist) freelist=alloc_new_nodes(bsize);
    p2=freelist;
    freelist=p2->Next;
    p2->col=y;
    p2->val=z;
    p2->Next=p;
    *pp=p2;
  }
}

static void free_sample_effects(void)
{
  if(XX) free(XX);
  if(zero) free(zero);
  if(B) free(B);
  if(entry) free(entry);
  if(order) free(order);
  if(full_xx) {
    if(full_xx[0]) free(full_xx[0]);
    free(full_xx);
  }
  XX=0;
  B=0;
  entry=0;
  order=0;
  zero=0;
  XX_size=entry_size=0;
  full_xx=0;
  init_flag=1;
}

static void init_sample_effects(void)
{
  int i,j,k,k1,k2,k3,type,rec,nrec,n_lev1,n_var1,comp,*order1,*mat;
  struct id_data *data;
  struct SparseMatRec *AI;
  struct Off *p,*p1;
  struct Model *mod;
  struct Id_Record *id_array;
	
  id_array=loki->pedigree->id_array;
  init_flag=1;
  if(full_store) {
    i=full_store*(full_store+1)/2;
    if(!(full_xx=malloc(sizeof(void *)*full_store))) ABT_FUNC(MMsg);
    if(!(full_xx[0]=malloc(sizeof(double)*i))) ABT_FUNC(MMsg);
    for(i=1;i<full_store;i++) full_xx[i]=full_xx[i-1]+i;
  }
  if(atexit(free_sample_effects)) message(WARN_MSG,"Unable to register exit function free_sample_effects()\n");
  n_lev=n_var=0;
  poly_start=mk_start=rand_start=cov_start=-1;
  mod=loki->models->models;
  /* Count variables and levels (ignoring QTL's and mean) */
  for(i=k=0;k<mod->n_terms;k++) {
    type=mod->term[k].vars[0].type;
    if(type&ST_MARKER) {
      if(mk_start<0) mk_start=n_lev;
      n_lev+=mod->term[k].df;
      n_var++;
    }
  }
  if(mod->polygenic_flag) {
    poly_start=n_lev;
    n_lev+=loki->pedigree->ped_size;
    n_var++;
  }
  if(loki->models->n_random) for(k=0;k<mod->n_terms;k++) {
    type=mod->term[k].vars[0].type;
    if(type&ST_RANDOM) {
      if(rand_start<0) rand_start=n_lev;
      n_lev+=mod->term[k].df;
      n_var++;
    }
  }
  for(k=0;k<mod->n_terms;k++) {
    type=mod->term[k].vars[0].type;
    if(type&(ST_MARKER|ST_TRAITLOCUS|ST_ID|ST_RANDOM)) continue;
    if(cov_start<0) cov_start=n_lev;
    n_lev+=mod->term[k].df;
    n_var++;
  }
  if(!n_lev) return;
  /* Allocate X'X Matrix */
  /* Make it a bit bigger so we can use it later */
  XX_size=(int)(1.1*(double)(n_lev+2+loki->params.max_tloci*2));
  if(!(XX=malloc(sizeof(struct Diag)*XX_size))) ABT_FUNC(MMsg);
  if(!(zero=malloc(sizeof(int)*XX_size))) ABT_FUNC(MMsg);
  if(!(B=malloc(sizeof(double)*XX_size))) ABT_FUNC(MMsg);
  for(i=0;i<n_lev;i++) {
    XX[i].First=0;
    XX[i].count=0;
  }
  entry_size=(int)(1.1*(double)n_var+2+loki->params.max_tloci);
  if(!(entry=malloc(sizeof(struct entry)*entry_size))) ABT_FUNC(MMsg);
  freelist=alloc_new_nodes(bsize);
  /* Assemble the X'X matrix with non-zero elements.  Need to do this so we can
     call min_deg() and get an order for factorizing the matrix */
  for(i=0;i<loki->pedigree->ped_size;i++) if(id_array[i].res[0]) {
    nrec=id_array[i].n_rec;
    if(!nrec) nrec=1;
    for(rec=0;rec<nrec;rec++) {
      /* Set up contribution vector */
      n_lev1=n_var1=0;
      /* Candidate genes */
      for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_MARKER) {
	  k1=mod->term[k].vars[0].var_index;
	  k2=loki->markers->marker[k1].locus.gt[i]-1;
#ifdef DEBUG
	  if(k2<0) ABT_FUNC("Internal error - candidate gene with unsampled genotype\n");
#endif
	  if(k2--)	{
	    entry[n_var1].pos=n_lev1+k2;
	  } else entry[n_var1].pos=-1;
	  n_var1++;
	  n_lev1+=mod->term[k].df;
	}
      }
      /* Polygenic effect */
      if(mod->polygenic_flag) for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_ID) {
	  entry[n_var1].pos=n_lev1+i; 
	  n_lev1+=loki->pedigree->ped_size;
	  n_var1++;
	}
      }
      /* Additional random effects */
      if(loki->models->n_random) for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_RANDOM) {
	  k1=mod->term[k].vars[0].var_index;
	  data=0;
	  if(type&ST_CONSTANT)	{
	    if(id_array[i].data) data=id_array[i].data+k1;
	  } else if(id_array[i].data1) data=id_array[i].data1[rec]+k1;
	  if(type&ST_FACTOR) { /* Discrete factors */
	    k2=(int)data->data.value-1;
#ifdef DEBUG
	    if(k2<0 || k2>mod->term[k].df) ABT_FUNC("Internal error - illegal factor value\n");
#endif
	    entry[n_var1].pos=n_lev1+k2;
	  } 
#ifdef DEBUG
	  else ABT_FUNC("Illegal Random type\n");
#endif
	  n_lev1+=mod->term[k].df;
	  n_var1++;
	}
      }
      /* Other Covariates */
      for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&(ST_ID|ST_TRAITLOCUS|ST_RANDOM|ST_MARKER)) continue;
	k1=mod->term[k].vars[0].var_index;
	data=0;
	if(type&ST_CONSTANT)	{
	  if(id_array[i].data) data=id_array[i].data+k1;
	} else if(id_array[i].data1) data=id_array[i].data1[rec]+k1;
	if(type&ST_FACTOR) { /* Discrete factors */
	  k2=(int)data->data.value-1;
#ifdef DEBUG
	  if(k2<0 || k2>mod->term[k].df) ABT_FUNC("Internal error - illegal factor value\n");
#endif
	  if(k2--)	entry[n_var1].pos=n_lev1+k2;
	  else entry[n_var1].pos=-1;
	} else { /* Continuous factors */
	  entry[n_var1].pos=n_lev1;
	}
	n_var1++;
	n_lev1+=mod->term[k].df;
      }
#ifdef DEBUG
      if(n_var1!=n_var || n_lev1!=n_lev) {
	(void)fprintf(stderr,"%d %d  %d %d\n",n_var,n_var1,n_lev,n_lev1);
	ABT_FUNC("Internal error - level mismatch\n");
      }
#endif
      /* Add contributions to X'X */
      k=mod->var.var_index;
      for(k=0;k<n_var;k++) {
	k1=entry[k].pos;
	if(k1>=0) {
	  XX[k1].val=1.0;
	  for(k2=k+1;k2<n_var;k2++) if((k3=entry[k2].pos)>=0) {
	    add_XX_entry1(k3,k1,1.0);
	  }
	}
      }
    }
  }
  if(mod->polygenic_flag) {
    /* Add A^-1 Matrix */
    j=poly_start;
    for(comp=0;comp<loki->pedigree->n_comp;comp++) {
      AI=loki->models->AIMatrix[comp];
      for(i=0;i<loki->pedigree->comp_size[comp];i++) {
	k=j+i;
	for(k1=AI[i].x;k1<AI[i+1].x;k1++) {
	  k2=AI[k1].x+j;
	  add_XX_entry1(k,k2,1.0);
	}
      }
      j+=i;
    }
  }
  /* Assemble matrix for min_deg */
  /* Count non-zero off-diagonal elements */
  for(j=i=0;i<n_lev;i++) j+=XX[i].count;
  if(!(order=malloc(sizeof(int)*n_lev))) ABT_FUNC(MMsg);
  if(!(order1=malloc(sizeof(int)*(n_lev*2+1+j)))) ABT_FUNC(MMsg);
  mat=order1+n_lev;
  k3=n_lev+1;
  for(i=0;i<n_lev;i++) {
    mat[i]=k3;
    p=p1=XX[i].First;
    while(p) {
      mat[k3++]=p->col;
      p1=p;
      p=p->Next;
    }
    if(p1) {
      p1->Next=freelist;
      freelist=XX[i].First;
    }
  }
  mat[i]=k3;
  /* Get factorization order */
  min_deg(n_lev,mat,order1,0,0);
  for(i=0;i<n_lev;i++) order[order1[i]]=i+1;
  free(order1);
  min_deg(0,0,0,0,0);
}

void sample_effects(struct loki *lk)
{
  int i,j,k,k1,k2,k3,mtype,type,idx,rec,nrec,n_qt,n_qtlev,n_var1,n_lev1,comp;
  int b_var,b_lev,t_var,t_lev,qt_start;
  double y,z,z1,ss,ssn,wt,wt1,*tdp,*tdp1,*tdp2,*tdp3;
  struct id_data *data;
  struct SparseMatRec *AI;
  struct Id_Record *id;
  struct Off *p,*p1,*p2,*p3,*p4,**pp;
  struct Model *mod;
  struct Id_Record *id_array;
	
  if(!lk->models->models) return;
  loki=lk;
  /* Get Ordering of equations */
  if(!init_flag) init_sample_effects();
  id_array=loki->pedigree->id_array;
  mod=loki->models->models;
  n_qt=n_qtlev=0;
  /* Count QTL's */
  for(j=0;j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag) {
    n_qt++;
    n_qtlev+=2;
  }
  /* Add one for data */
  b_var=n_qt+1;
  b_lev=n_qtlev+1;
  /* Add one for the grand mean if not fixed */
  if(loki->models->grand_mean_set[0]!=1) {
    b_lev++;
    b_var++;
  }
  qt_start=b_lev-n_qtlev;
  t_var=b_var+n_var;
  t_lev=b_lev+n_lev;
  if(!(t_var&&t_lev)) return;
  /* Check size of X'X Matrix */
  if(t_lev>XX_size) {
    XX_size=(int)((double)t_lev*1.1);
    if(!(XX=realloc(XX,sizeof(struct Diag)*XX_size))) ABT_FUNC(MMsg);
    if(!(zero=realloc(zero,sizeof(int)*XX_size))) ABT_FUNC(MMsg);
    if(!(B=realloc(B,sizeof(double)*XX_size))) ABT_FUNC(MMsg);
  }
  for(i=0;i<t_lev;i++) {
    XX[i].First=0;
    XX[i].val=0;
  }
  if(full_store) {
    j=full_store>t_lev?t_lev:full_store;
    tdp=full_xx[0];
    k=j*(j+1)/2;
    while(k--) *(tdp++)=0;
  }
  if(t_var>entry_size) {
    entry_size=(int)(1.1*(double)t_var);
    if(!(entry=realloc(entry,sizeof(struct entry)*entry_size))) ABT_FUNC(MMsg);
  }
  /* Assemble the X'X matrix */
  mtype=mod->var.type;
  idx=mod->var.var_index;
  id=id_array;
  for(i=0;i<loki->pedigree->ped_size;i++,id++) if(id->res[0]) {
    nrec=id->n_rec;
    for(rec=0;rec<nrec;rec++) {
      wt=loki->models->use_student_t?id->vv[0][rec]:1.0;
      /* Load raw data value into y */
      y=0.0;
      if(mtype&ST_CONSTANT) {
	data=id->data+idx;
	if(data->flag&2) y=id->pseudo_qt[0][0];
      } else {
	data=id->data1[rec]+idx;
	if(data->flag&2) y=id->pseudo_qt[0][rec];
      }
      if(data->flag&ST_INTTYPE) y+=(double)data->data.value;
      else y+=data->data.rvalue;
      /* Set up contribution vector */
      n_lev1=n_var1=1;
      /* ...mean */
      if(loki->models->grand_mean_set[0]!=1) {
	entry[n_var1].val=1.0;
	entry[n_var1++].pos=n_lev1++;
      } else y-=loki->models->grand_mean[0];
      /* ...QTL's */
      if(n_qt) {
	for(j=0;j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag)	{
	  k2=loki->models->tlocus[j].gt[i]-1;
	  if(k2--)	{
	    entry[n_var1].val=1.0;
	    entry[n_var1++].pos=n_lev1+k2;
	  }
	  n_lev1+=2;
	}
      }
      b_var=n_var1;
      /* Candidate genes */
      if(mk_start>=0) for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_MARKER) { 
	  k1=mod->term[k].vars[0].var_index;
	  k2=loki->markers->marker[k1].locus.gt[i]-1;
#ifdef DEBUG
	  if(k2<0) ABT_FUNC("Internal error - candidate gene with unsampled genotype\n");
#endif
	  if(k2--)	{
	    entry[n_var1].val=1.0;
	    entry[n_var1++].pos=n_lev1+k2;
	  }
	  n_lev1+=mod->term[k].df;
	}
      }
      /* Polygenic effect */
      if(mod->polygenic_flag) for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_ID) {
	  entry[n_var1].val=1.0;
	  entry[n_var1++].pos=n_lev1+i; /* Polygenic effect */
	  n_lev1+=loki->pedigree->ped_size;
	}
      }
      /* Additional random effects */
      if(loki->models->n_random) for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&ST_RANDOM) {
	  k1=mod->term[k].vars[0].var_index;
	  data=0;
	  if(type&ST_CONSTANT)	{
	    if(id->data) data=id->data+k1;
	  } else if(id->data1) data=id->data1[rec]+k1;
	  if(type&ST_FACTOR) { /* Discrete factors */
	    k2=(int)data->data.value-1;
#ifdef DEBUG
	    if(k2<0 || k2>mod->term[k].df) ABT_FUNC("Internal error - illegal factor value\n");
#endif
	    entry[n_var1].pos=n_lev1+k2;
	    entry[n_var1].val=1.0;
	  }
#ifdef DEBUG
	  else ABT_FUNC("Illegal Random type\n");
#endif
	  n_lev1+=mod->term[k].df;
	  n_var1++;
	}
      }
      /* ...other effects */
      for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&(ST_TRAITLOCUS|ST_ID|ST_MARKER|ST_RANDOM)) continue; /* Already done */
	k1=mod->term[k].vars[0].var_index;
	data=0;
	if(type&ST_CONSTANT)	{
	  if(id->data) data=id->data+k1;
	} else if(id->data1) data=id->data1[rec]+k1;
	if(type&ST_FACTOR) { /* Discrete factors */
	  k2=(int)data->data.value-1;
#ifdef DEBUG
	  if(k2<0 || k2>mod->term[k].df) ABT_FUNC("Internal error - illegal factor value\n");
#endif
	  if(k2--)	{
	    entry[n_var1].pos=n_lev1+k2;
	    entry[n_var1++].val=1.0;
	  }
	} else { /* Continuous factors */
	  if(data->flag&ST_INTTYPE) z=(double)data->data.value;
	  else z=data->data.rvalue;
	  entry[n_var1].pos=n_lev1;
	  entry[n_var1++].val=z;
	}
	n_lev1+=mod->term[k].df;
      }
#ifdef DEBUG
      if(n_lev1!=t_lev) ABT_FUNC("Internal error - level mismatch\n");
#endif
      id->res[0][rec]=y;
      entry[0].pos=0;
      entry[0].val=y;
      /* Add contributions to X'X */
      /* Reorder equations */
      k3=0;
      for(k=b_var;k<n_var1;k++) {
	k1=entry[k].pos;
#ifdef DEBUG
	if(k1>=b_lev) { 
	  k1=n_lev+b_lev-order[(k1-b_lev)];
	} else ABT_FUNC("Illegal equation position\n");
	if(k1<0 || k1>=t_lev) ABT_FUNC("Illegal equation position\n");
#else
	k1=n_lev+b_lev-order[(k1-b_lev)];
#endif
	entry[k].pos=k1;
      }
      wt1=1.0/(wt*loki->models->residual_var[0]);
      if(t_lev<=full_store) {
	for(k=0;k<n_var1;k++) {
	  k1=entry[k].pos;
	  if(k1>=0) {
	    z=entry[k].val;
	    if(fabs(z)<DBL_EPSILON) continue;
	    tdp=full_xx[k1];
	    tdp[k1]+=z*z*wt1;
	    for(k2=k+1;k2<n_var1;k2++) {
	      k3=entry[k2].pos;
	      y=entry[k2].val;
	      if(k3<k1) tdp[k3]+=z*y*wt1;
	      else full_xx[k3][k1]+=z*y*wt1;
	    }
	  }
	}
      } else {
	for(k=0;k<n_var1;k++) {
	  k1=entry[k].pos;
	  if(k1>=0) {
	    z=entry[k].val;
	    if(fabs(z)<DBL_EPSILON) continue;
	    if(k1<full_store) {
	      tdp=full_xx[k1];
	      tdp[k1]+=z*z*wt1;
	      for(k2=k+1;k2<n_var1;k2++) {
		k3=entry[k2].pos;
		y=entry[k2].val;
		if(k3<k1) tdp[k3]+=z*y*wt1;
		else {
		  if(k3>=full_store) add_XX_entry(k3,k1,z*y*wt1);
		  else full_xx[k3][k1]+=z*y*wt1;
		}
	      }
	    } else {
	      XX[k1].val+=z*z*wt1;
	      for(k2=k+1;k2<n_var1;k2++) {
		k3=entry[k2].pos;
		y=entry[k2].val;
		if(k3<k1) add_XX_entry(k1,k3,z*y*wt1);
		else add_XX_entry(k3,k1,z*y*wt1);
	      }
	    }
	  }
	}
      }
    }
  }
  if(mod->polygenic_flag) {
    z=1.0/loki->models->additive_var[0];
    /* Add A^-1 Matrix */
    y=0.0;
    j=poly_start;
    for(comp=0;comp<loki->pedigree->n_comp;comp++) {
      AI=loki->models->AIMatrix[comp];
      for(i=0;i<loki->pedigree->comp_size[comp];i++) {
	k=n_lev+b_lev-order[j+i];
	if(k<full_store) {
	  tdp=full_xx[k];
	  tdp[k]+=z*AI[i].val;
	  for(k1=AI[i].x;k1<AI[i+1].x;k1++) {
	    k2=AI[k1].x+j;
	    k2=n_lev+b_lev-order[k2];
	    if(k2<k) tdp[k2]+=z*AI[k1].val;
	    else {
	      if(k2>=full_store) add_XX_entry(k2,k,z*AI[k1].val);
	      else full_xx[k2][k]+=z*AI[k1].val;
	    }
	  }
	} else {
	  XX[k].val+=z*AI[i].val;
	  for(k1=AI[i].x;k1<AI[i+1].x;k1++) {
	    k2=AI[k1].x+j;
	    k2=n_lev+b_lev-order[k2];
	    if(k2<k) add_XX_entry(k,k2,z*AI[k1].val);
	    else add_XX_entry(k2,k,z*AI[k1].val);
	  }
	}
      }
      j+=i;
    }
  }
  /* Add random terms for additional random effects */
  if(loki->models->n_random) for(j=rand_start,i=k=0;k<mod->n_terms;k++) {
    type=mod->term[k].vars[0].type;
    if(type&ST_RANDOM) {
      z=1.0/loki->models->c_var[i][0];
      for(k1=0;k1<mod->term[k].df;k1++) {
	k2=n_lev+b_lev-order[j++];
	if(k2>=full_store) XX[k2].val+=z;
	else full_xx[k2][k2]+=z;
      }
      i++;
    }
  }
  z=1.0/loki->models->tau[0];
  if(n_qt) {
    /* Add random terms for QTLs */
    if(b_lev<=full_store) {
      for(j=qt_start,i=0;i<loki->params.n_tloci;i++) if(loki->models->tlocus[i].flag)	{
	for(k=0;k<2;k++,j++) full_xx[j][j]+=z;
      }
    } else {
      for(j=qt_start,i=0;i<loki->params.n_tloci;i++) if(loki->models->tlocus[i].flag)	{
	for(k=0;k<2;k++,j++) {
	  if(j>=full_store)	XX[j].val+=z;
	  else full_xx[j][j]+=z;
	}
      }
    }
  }
  /* and for candidate genes */
  if(mk_start>=0) {
    for(j=mk_start,k=0;k<mod->n_terms;k++) {
      type=mod->term[k].vars[0].type;
      if(type&ST_MARKER) {
	for(k1=0;k1<mod->term[k].df;k1++) {
	  k2=n_lev+b_lev-order[j++];
	  if(k2>=full_store) XX[k2].val+=z;
	  else full_xx[k2][k2]+=z;
	}
      }
    }
  }
  /* Gaussian elimination step  - sparse region */
  for(i=t_lev-1;i>=full_store && i>0;i--) {
    y=XX[i].val;
    if(fabs(y)<1.0e-12) {
      zero[i]=1;
      continue;
    } else zero[i]=0;
    if(y<DBL_EPSILON) ABT_FUNC("Effects matrix not positive definite\n");
    y=1.0/y;
    /* Absorb pivot row into remaining rows */
    if((p=XX[i].First)) {
      z=p->val;
      j=p->col;
      if(j<full_store) full_xx[j][j]-=z*z*y;
      else XX[j].val-=z*z*y;
      p2=p;
      while((p=p->Next))	{
	j=p->col;
	if(j<full_store) {
	  p1=p2;
	  z1=p->val;
	  z=-z1*y;
	  tdp=full_xx[j];
	  while(p1!=p) {
	    tdp[p1->col]+=p1->val*z;
	    p1=p1->Next;
	  }
	  tdp[j]+=z1*z;
	} else break;
      }
      while(p)	{
	j=p->col;
	p1=p2;
	z1=p->val;
	z=-z1*y;
	pp=&XX[j].First;
	k1=p1->col;
	while((p3=*pp)) {
	  if(p3->col==k1) {
	    p3->val+=p1->val*z;
	    p1=p1->Next;
	    if(p1==p) break;
	    pp=&p3->Next;
	    k1=p1->col;
	  } else if(p3->col>k1) {
	    if(!freelist) freelist=alloc_new_nodes(bsize);
	    p4=freelist;
	    freelist=p4->Next;
	    p4->col=k1;
	    p4->val=p1->val*z;
	    p4->Next=p3;
	    *pp=p4;
	    p1=p1->Next;
	    if(p1==p) break;
	    pp=&p4->Next;
	    k1=p1->col;
	  } else pp=&p3->Next;
	}
	if(p1!=p) {
	  do {
	    if(!freelist) freelist=alloc_new_nodes(bsize);
	    p3=freelist;
	    freelist=p3->Next;
	    p3->col=p1->col;
	    p3->val=p1->val*z;
	    *pp=p3;
	    pp=&p3->Next;
	    p1=p1->Next;
	  } while(p1!=p);
	  *pp=0;
	}
	XX[j].val-=z1*z1*y;
	p=p->Next;
      }
    }
  }
  /* Gaussian elimination - full_stored region */
  for(;i>0;i--) {
    tdp3=tdp=full_xx[i];
    y=tdp[i];
    if(fabs(y)<1.0e-12) {
      zero[i]=1;
      continue;
    } else zero[i]=0;
    if(y<DBL_EPSILON) 
      ABT_FUNC("Effects matrix not positive definite\n");
    y=1.0/y;
    /* Absorb pivot row into remaining rows */
    tdp1=*full_xx;
    z1=*(tdp++);
    *(tdp1++)-=z1*z1*y;
    for(j=1;j<i;j++) {
      z1=*(tdp++);
      z=-z1*y;
      tdp2=tdp3;
      k=j;
      while(k--) *(tdp1++)+=*(tdp2++)*z;
      *(tdp1++)+=z1*z;
    }
  }
  /* And now go backwards, sampling as we go */
  k=full_store>t_lev?t_lev:full_store;
  for(i=1;i<k;i++) if(!zero[i]) {
    tdp=full_xx[i];
    y=*tdp;
    for(j=1;j<i;j++) if(!zero[j]) y-=B[j]*tdp[j];
    z=1.0/tdp[i];
    B[i]=y*z+snorm()*sqrt(z);
    if(loki->models->no_overdominant && n_qt && i>qt_start && i<qt_start+n_qtlev) {
      for(j=qt_start,k1=0;i>j && k1<loki->params.n_tloci;k1++) if(loki->models->tlocus[k1].flag)	{
	if(i==j+1) {
	  z=B[i]?B[i-1]/B[i]:0.0;
	  if(z<0.0 || z>1.0) {
	    B[i-1]=loki->models->tlocus[k1].eff[0][0];
	    B[i]=loki->models->tlocus[k1].eff[0][1];
	  }
	  break;
	}
	j+=2;
      }
    }
  }
  for(;i<t_lev;i++) if(!zero[i]) {
    if((p=XX[i].First)) {
      p1=p;
      if(!p->col) {
	y=p->val;
	p=p->Next;
      } else y=0.0;
      while(p) {
	j=p->col;
	if(!zero[j]) y-=B[j]*p->val;
	p1=p;
	p=p->Next;
      }
      /* Free up what we don't need */
      if(p1) {
	p1->Next=freelist;
	freelist=XX[i].First;
      }
    } else y=0.0;
    z=1.0/XX[i].val;
    B[i]=y*z+snorm()*sqrt(z);
    if(loki->models->no_overdominant && n_qt && i>qt_start && i<qt_start+n_qtlev) {
      for(j=qt_start,k1=0;i>j && k1<loki->params.n_tloci;k1++) if(loki->models->tlocus[k1].flag)	{
	if(i==j+1) {
	  z=B[i]?B[i-1]/B[i]:0.0;
	  if(z<0.0 || z>1.0) {
	    B[i-1]=loki->models->tlocus[k1].eff[0][0];
	    B[i]=loki->models->tlocus[k1].eff[0][1];
	  }
	  break;
	}
	j+=2;
      }
    }
  }
  /* Store new effect samples */
  /* Grand Mean */
  if(loki->models->grand_mean_set[0]!=1) {
    if(!zero[1]) loki->models->grand_mean[0]=B[1];
  }
  /* QTLs */
  if(n_qt) {
    for(j=qt_start,i=0;i<loki->params.n_tloci;i++) if(loki->models->tlocus[i].flag)	{
      for(k=0;k<2;k++) loki->models->tlocus[i].eff[0][k]=B[j++];
    }
  }
  /* candidate genes */
  if(mk_start>=0) {
    for(j=mk_start,k=0;k<mod->n_terms;k++) {
      type=mod->term[k].vars[0].type;
      if(type&ST_MARKER) {
	for(k1=0;k1<mod->term[k].df;k1++) {
	  k2=n_lev+b_lev-order[j++];
	  mod->term[k].eff[k1]=B[k2];
	}
      }
    }
  }
  /* additional random effects */
  if(loki->models->n_random) for(j=rand_start,i=k=0;k<mod->n_terms;k++) {
    type=mod->term[k].vars[0].type;
    if(type&ST_RANDOM) {
      for(k1=0;k1<mod->term[k].df;k1++) {
	k2=n_lev+b_lev-order[j++];
	mod->term[k].eff[k1]=B[k2];
      }
    }
  }
  /* Polygenic component */
  if(mod->polygenic_flag) {
    j=poly_start;
    id=id_array;
    for(i=0;i<loki->pedigree->ped_size;i++,id++) {
      k=n_lev+b_lev-order[j++];
      z=B[k];
      id->bv[0]=z;
      id->bvsum[0]+=z;
      id->bvsum2[0]+=z*z;
    }
    loki->params.bv_iter++;
  }
  /* and other covariates effects */
  if(cov_start>=0) {
    j=cov_start;
    for(k=0;k<mod->n_terms;k++) {
      type=mod->term[k].vars[0].type;
      if(type&(ST_TRAITLOCUS|ST_ID|ST_MARKER|ST_RANDOM)) continue; /* Already done */
      for(k1=0;k1<mod->term[k].df;k1++) {
	k2=n_lev+b_lev-order[j++];
	mod->term[k].eff[k1]=B[k2];
      }
    }
  }
  ss=ssn=0.0;
  /* Correct loki->models->residuals for new effects */
  id=id_array;
  for(i=0;i<loki->pedigree->ped_size;i++,id++) if(id->res[0]) {
    nrec=id->n_rec;
    for(rec=0;rec<nrec;rec++) {
      if(mtype&ST_CONSTANT) y=id->res[0][0];
      else y=id->res[0][rec];
      if(loki->models->grand_mean_set[0]!=1) y-=loki->models->grand_mean[0];
      if(mod->polygenic_flag) y-=id->bv[0];
      for(j=0;j<loki->params.n_tloci;j++) if(loki->models->tlocus[j].flag)	{
	k2=loki->models->tlocus[j].gt[i]-1;
	if(k2) y-=loki->models->tlocus[j].eff[0][k2-1];
      }
      for(k=0;k<mod->n_terms;k++) {
	type=mod->term[k].vars[0].type;
	if(type&(ST_ID|ST_TRAITLOCUS)) continue;
	k1=mod->term[k].vars[0].var_index;
	if(type&ST_MARKER) {
	  k2=loki->markers->marker[k1].locus.gt[i]-1;
	  if(k2) y-=mod->term[k].eff[k2-1];
	} else {
	  data=0;
	  if(type&ST_CONSTANT)	{
	    if(id->data) data=id->data+k1;
	  } else if(id->data1) data=id->data1[rec]+k1;
#ifdef DEBUG
	  if(!data) ABT_FUNC("Internal error - null data pointer\n");
#endif
	  if(type&ST_FACTOR) {
	    k2=(int)data->data.value-1;
	    if(!(type&ST_RANDOM)) k2--;
	    if(k2>=0) y-=mod->term[k].eff[k2];
	  } else {
	    if(data->flag&ST_INTTYPE) z=(double)data->data.value;
	    else z=data->data.rvalue;
	    y-=mod->term[k].eff[0]*z;
	  }
	}
      }
      id->res[0][rec]=y;
      ss+=y*y;
      ssn+=1.0;
    }
  }
}
