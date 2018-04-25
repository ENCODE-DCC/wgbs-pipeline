/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                        August 1997                                       *
 *                                                                          *
 * loki_pen.c:                                                              *
 *                                                                          *
 * Penetrance routines for peeling calculations                             *
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
#include "loki_peel.h"
#include "loki_utils.h"

void penetrance(double *val,int id,struct Locus *loc,int n_all,int n_bits,const struct loki *loki)
{
  int i,j,k,mtype,rec,nrec,censflag,idx,censor_mode,use_student_t;
  double p,y,m,*eff,kon1,kon2,sd,wt,res;
  struct Id_Record *idr;
  struct Model *mod;

  idr=loki->pedigree->id_array+id;
  if(!idr->res[0]) return;
  eff=loc->eff[0];
#ifdef DEBUG
  if(n_bits) j=1<<(n_bits+n_bits);
  else j=n_all*n_all;
  for(p=0.0,i=0;i<j;i++) p+=val[i];
  if(p<=0.0) {
    (void)fprintf(stderr,"penetrance() called with invalid R-Function (p=%g) - ",p);
    print_orig_id(stderr,id+1);
    (void)fputc('\n',stderr);
    ABT_FUNC("Aborting\n");
  }
#endif
  res=loki->models->residual_var[0];
  kon1=1.0/(2.0*res);
  kon2=sqrt(kon1/M_PI);
  sd=sqrt(res);
  mod=loki->models->models;
  idx=mod->var.var_index;
  mtype=mod->var.type;
  nrec=idr->n_rec;
  if(!nrec) nrec=1;
  censor_mode=loki->models->censor_mode;
  use_student_t=loki->models->use_student_t;
  for(rec=0;rec<nrec;rec++) {
    censflag=0;
    y=idr->res[0][rec];
    wt=use_student_t?1.0/idr->res[0][rec]:1.0;
    if(!censor_mode && (mtype&ST_CENSORED)) {
      if(mtype&ST_CONSTANT) {
	if(idr->data[idx].flag&2) censflag=1;
      } else {
	if(idr->data1[rec][idx].flag&2) censflag=1;
      }
    }
    /* Add on current genotype effect */
    if(loc->flag&LOCUS_SAMPLED) {
      k=loc->gt[id]-1;
      if(k) y+=eff[k-1];
    }
    if(censflag) {
      y-=idr->pseudo_qt[0][rec];
      p=.5*erfc(y*sqrt(wt)/sd);
    } else p=kon2*exp(-y*y*kon1*wt);
    val[0]*=p;
    if(!n_bits) { /* Non-sparse representation */
      for(k=0,i=1;i<n_all;i++) for(j=0;j<=i;j++) {
	m=y-eff[k++];
	if(censflag) p=.5*erfc(m*sqrt(wt)/sd);
	else p=kon2*exp(-m*m*wt*kon1);
	val[i*n_all+j]*=p;
	if(i!=j) val[j*n_all+i]*=p;
      }
    } else {
      for(k=0,i=1;i<n_all;i++) for(j=0;j<=i;j++) { /* Sparse representation */
	m=y-eff[k++];
	if(censflag) p=.5*erfc(m*sqrt(wt)/sd);
	else p=kon2*exp(-m*m*wt*kon1);
	val[(i<<n_bits)|j]*=p;
	if(i!=j) val[(j<<n_bits)|i]*=p;
      }
    }
  }
}

void s_penetrance(double *val,int id,struct Locus *loc,const struct loki *loki)
{
  int k,mtype,rec,nrec,censflag,idx,censor_mode,use_student_t;
  double p,y,m,*eff,kon1,kon2,kon1a,sd,sd1,wt,res;
  struct Id_Record *idr;
	
#ifdef DEBUG
  int i;
	
  for(p=0.0,i=0;i<4;i++) p+=val[i];
  if(p<=0.0) {
    (void)fprintf(stderr,"s_penetrance() called with zero function - ");
    print_orig_id(stderr,id+1);
    (void)fputc('\n',stderr);
    ABT_FUNC("Aborting\n");
  }
#endif
  eff=loc->eff[0];
  idr=loki->pedigree->id_array+id;
  nrec=idr->n_rec;
  censor_mode=loki->models->censor_mode;
  use_student_t=loki->models->use_student_t;
  res=loki->models->residual_var[0];
  if(!nrec) nrec=1;
  mtype=loki->models->models[0].var.type;
  kon1=1.0/(2.0*res);
  kon2=sqrt(kon1/M_PI);
  if(!(mtype&ST_CENSORED) || censor_mode) {
    for(rec=0;rec<nrec;rec++) {
      y=idr->res[0][rec];
      kon1a=use_student_t?kon1/idr->vv[0][rec]:kon1;
      if(loc->flag&LOCUS_SAMPLED) {
	k=loc->gt[id]-1;
	if(k) y+=eff[k-1];
      }
      val[0]*=kon2*exp(-y*y*kon1a);
      m=y-eff[0];
      p=kon2*exp(-m*m*kon1a);
      val[1]*=p;
      val[2]*=p;
      m=y-eff[1];
      val[3]*=kon2*exp(-m*m*kon1a);
    }
  } else {
    sd=sqrt(res*2.0);
    idx=loki->models->models[0].var.var_index;
    for(rec=0;rec<nrec;rec++) {
      censflag=0;
      y=idr->res[0][rec];	
      wt=use_student_t?1.0/idr->vv[0][rec]:1.0;
      sd1=sqrt(wt)/sd;
      kon1a=kon1*wt;
      if(mtype&ST_CONSTANT) {
	if(idr->data[idx].flag&2) censflag=1;
      } else {
	if(idr->data1[rec][idx].flag&2) censflag=1;
      }
      if(loc->flag&LOCUS_SAMPLED) {
	k=loc->gt[id]-1;
	if(k) y+=eff[k-1];
      }
      if(censflag) {
	y-=idr->pseudo_qt[0][rec];
	val[0]*=.5*erfc(y*sd1);
	m=y-eff[0];
	p=.5*erfc(m*sd1);
	val[1]*=p;
	val[2]*=p;
	m=y-eff[1];
	val[3]*=.5*erfc(m*sd1);
      } else {
	val[0]*=kon2*exp(-y*y*kon1a);
	m=y-eff[0];
	p=kon2*exp(-m*m*kon1a);
	val[1]*=p;
	val[2]*=p;
	m=y-eff[1];
	val[3]*=kon2*exp(-m*m*kon1a);
      }
    }
  }
}

void s_penetrance1(double *val,int id,struct Locus *loc,const struct loki *loki)
{
  int k;
  double p,y,m,*eff,kon1,kon2;
  struct Id_Record *idr;
	
#ifdef DEBUG
  int i;
	
  for(p=0.0,i=0;i<4;i++) p+=val[i];
  if(p<=0.0) {
    (void)fprintf(stderr,"penetrance() called with zero function - ");
    print_orig_id(stderr,id+1);
    (void)fputc('\n',stderr);
    ABT_FUNC("Aborting\n");
  }
#endif
  eff=loc->eff[0];
  kon1=-1.0/(2.0*loki->models->residual_var[0]);
  kon2=sqrt(-kon1/M_PI);
  idr=loki->pedigree->id_array+id;
  y=idr->res[0][0];
  if((k=loc->gt[id])>1) y+=eff[k-2];
#ifdef TRACE_PEEL
  if(CHK_PEEL(TRACE_LEVEL_4)) {
    (void)fputs("Penetrance routine\n",stdout);
    printf("y=%g, k=%d, eff[0]=%g, eff[1]=%g, kon1=%g, kon2=%g\n",y,k,eff[0],eff[1],kon1,kon2);
  }
#endif
  val[0]*=kon2*exp(y*y*kon1);
  m=y-eff[0];
  p=kon2*exp(m*m*kon1);
  val[1]*=p;
  val[2]*=p;
  m=y-eff[1];
  val[3]*=kon2*exp(m*m*kon1);
#ifdef DEBUG
  for(p=0.0,i=0;i<4;i++) p+=val[i];
  if(p<=0.0) {
    (void)fprintf(stderr,"penetrance() returning with zero function - ");
    print_orig_id(stderr,id+1);
    (void)fprintf(stderr," y=%g, sd=%g",y,sqrt(loki->models->residual_var[0]));
    (void)fputc('\n',stderr);
    ABT_FUNC("Aborting\n");
  }
#endif
}

double q_penetrance(int id,int gt,struct Locus *locp,const struct loki *loki)
{
  int k,mtype,rec,nrec,censflag,idx,censor_mode,use_student_t;
  double p,y,*eff,kon1,kon2,sd,wt,res;
  struct Id_Record *idr;
	
  eff=locp->eff[0];
  gt--;
  idr=loki->pedigree->id_array+id;
  censor_mode=loki->models->censor_mode;
  use_student_t=loki->models->use_student_t;
  nrec=idr->n_rec;
  if(!nrec) nrec=1;
  mtype=loki->models->models[0].var.type;
  p=0.0;
  res=loki->models->residual_var[0];
  kon1=1.0/(2.0*res);
  kon2=log(sqrt(kon1/M_PI));
  if(!(mtype&ST_CENSORED) || censor_mode) {
    for(rec=0;rec<nrec;rec++) {
      y=idr->res[0][rec];
      wt=use_student_t?1.0/idr->vv[0][rec]:1.0;
      if(locp->flag&LOCUS_SAMPLED) {
	k=locp->gt[id]-1;
	if(k) y+=eff[k-1];
      }
      if(gt) y-=eff[gt-1];
      p+=kon2-y*y*kon1*wt;
    }
  } else {
    sd=sqrt(res*2.0);
    idx=loki->models->models[0].var.var_index;
    for(rec=0;rec<nrec;rec++) {
      censflag=0;
      y=idr->res[0][rec];	
      wt=use_student_t?1.0/idr->vv[0][rec]:1.0;
      if(mtype&ST_CONSTANT) {
	if(idr->data[idx].flag&2) censflag=1;
      } else {
	if(idr->data1[rec][idx].flag&2) censflag=1;
      }
      if(locp->flag&LOCUS_SAMPLED) {
	k=locp->gt[id]-1;
	if(k) y+=eff[k-1];
      }
      if(censflag) {
	y-=idr->pseudo_qt[0][rec];
	if(gt) y-=eff[gt-1];
	p+=log(.5*erfc(y*sqrt(wt)/sd));
      } else {
	if(gt) y-=eff[gt-1];
	p+=kon2-y*y*kon1*wt;
      }
    }
  }
  return p;
}
