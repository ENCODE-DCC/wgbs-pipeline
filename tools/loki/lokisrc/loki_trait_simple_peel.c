#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "ranlib.h"
#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "get_par_probs.h"
#include "loki_trait_simple_peel.h"
#ifndef DBL_MAX
#define DBL_MAX MAXDOUBLE
#endif

/* Similar to loki_simple_peelop, but for a trait locus */
double loki_trait_simple_peelop(const struct Simple_Element *element,const int locus,const int s_flag,double **freq,struct R_Func *rf,trait_pen_func *trait_pen,struct loki *loki)
{
  int ids,idd,i,j,k,m,n,pivot,n_off,*off,kid,*ix,link,n_all,n_idx;
  int ix1[]={0,3,12,15};
  int ix2[]={0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15};
  int ix3[]={5,6,9,10};
  double prob=0.0,*tp,*tpp,p1,z,z1,*tmp,*tmp1,*qval,*mval,*pval,*peel_famval;
  struct peel_mem *work;
  struct Id_Record *id_array;
  struct Locus *loc;
	
  pivot=element->pivot-1;
#ifdef DEBUG
  if(pivot== -3) ABT_FUNC("Internal error - this peeling operation should not occur with a trait locus\n");
#endif
  work=&loki->peel->workspace;
  id_array=loki->pedigree->id_array;
  ids=element->sire-1;
  idd=element->dam-1;
  off=element->off;
  n_off=element->n_off;
  loc=&loki->models->tlocus[-1-locus];
  n_all=loc->n_alleles;
  n_idx=n_all*n_all;
#ifdef TRACE_PEEL
  if(CHK_PEEL(TRACE_LEVEL_1))	(void)printf("%s(%p,%d,%d)\n",__func__,(void *)element,locus,s_flag);
  if(CHK_PEEL(TRACE_LEVEL_2)) {
    if(family_id) {
      print_orig_family(stdout,off[0]+1,0);
      (void)fputc(' ',stdout);
    }
    print_orig_id1(stdout,ids+1);
    (void)fputc(',',stdout);
    print_orig_id1(stdout,idd+1);
    (void)fputc(' ',stdout);
    for(i=0;i<n_off;i++) {
      (void)fputc(i?',':'(',stdout);
      print_orig_id1(stdout,off[i]);
    }
    (void)fputs(") -> ",stdout);
    if(pivot==-2) {
      print_orig_id1(stdout,ids+1);
      (void)fputc(',',stdout);
      print_orig_id1(stdout,idd+1);
    } else if(pivot<-2) printf("[%d]",pivot);
    else print_orig_id1(stdout,pivot+1);
    (void)fputc('\n',stdout);
  }
#endif
  peel_famval=work->s2;
  qval=peel_famval+16; 
  if(ids<0) { /* Peeling singletons */
    for(m=0;m<n_off;m++)	{
      kid=off[m]-1;
      p1=get_trait_par_probs(qval,kid,locus,trait_pen,freq,rf,loki);
      if(p1<=0.0) {
	if(!(s_flag&1)) return -DBL_MAX;
	ABT_FUNC("Zero probability in peeling operation\n");
      }
      prob+=log(p1);
      if(s_flag) {
	do {
	  z=ranf();
	  p1=0.0;
	  for(i=0;i<n_idx;i++) if(qval[i])	{
	    p1+=qval[i];
	    if(z<=p1) break;
	  }
	} while(i==n_idx);
	id_array[kid].allele[X_MAT]=1+(i%n_all);
	id_array[kid].allele[X_PAT]=1+((i/n_all)%n_all);
	id_array[kid].flag|=(SAMPLED_MAT|SAMPLED_PAT);
      }
    }
    return prob;
  }
  if(s_flag && !(element->pivot)) 
    return loki_trait_simple_sample(element,locus,s_flag,freq,rf,trait_pen,loki);
  pval=qval+n_idx;
  mval=pval+n_idx; 
  if(idd!=pivot && pivot!= -2) {
    p1=get_trait_par_probs(mval,idd,locus,trait_pen,freq,rf,loki);
#ifdef DEBUG
    if(isnan(p1)) ABT_FUNC("Floating point error\n");
#endif
    if(p1<=0.0)	{
      if(!(s_flag&1)) {
	return -DBL_MAX;
      }
      ABT_FUNC("Zero probability in peeling operation\n");
    }
    prob+=log(p1);
  } else {
    tmp=mval;
    if((k=id_array[idd].rfp)>=0) { /* Insert Previously computed R_Func */
      tmp1=rf[k].p;
      for(j=0;j<4;j++) *(tmp++)=(*tmp1++);
    } else for(j=0;j<4;j++) *(tmp++)=1.0;
  }
  if(ids!=pivot && pivot!= -2) {
    p1=get_trait_par_probs(pval,ids,locus,trait_pen,freq,rf,loki);
#ifdef DEBUG
    if(isnan(p1)) ABT_FUNC("Floating point error\n");
#endif
    if(p1<=0.0)	{
      if(!(s_flag&1)) {
	return -DBL_MAX;
      }
      ABT_FUNC("Zero probability in peeling operation\n");
    }
    prob+=log(p1);
  } else {
    tmp=pval;
    if((k=id_array[ids].rfp)>=0) { /* Insert Previously computed R_Func */
      tmp1=rf[k].p;
      for(j=0;j<n_idx;j++) *(tmp++)=(*tmp1++);
    } else for(j=0;j<n_idx;j++) *(tmp++)=1.0;
  }
  tmp1=peel_famval;
  for(i=0;i<4;i++) {
    p1=pval[i];
    tmp=mval;
    for(k=0;k<4;k++) *(tmp1++)=p1*(*(tmp++));
  }
  link=loc->link_group;
  for(m=0;m<n_off;m++)	{
    kid=off[m]-1;
    if(kid==pivot) continue;
    tmp=qval;
    if((k=id_array[kid].rfp)>=0) { /* Insert Previously computed R_Func */
      tmp1=rf[k].p;
      for(j=0;j<4;j++) *(tmp++)=(*tmp1++);
    } else for(j=0;j<4;j++) *(tmp++)=1.0;
    if(id_array[kid].res[0]) trait_pen(qval,kid,loc,loki);
    p1=0.0;
    tmp=qval;
    for(k=0;k<4;k++) p1+=*(tmp++);
#ifdef DEBUG
    if(isnan(p1)) ABT_FUNC("Floating point error\n");
#endif
    if(p1<=0.0) {
      if(!(s_flag&1)) {
	return -DBL_MAX;
      }
      ABT_FUNC("Zero probability in peeling operation\n");
    }
    prob+=log(p1);
    z=1.0/p1;
    /* First do double homozygote configs */
    tmp=qval;
    ix=ix1;
    for(i=0;i<4;i++) peel_famval[*(ix++)]*=z*(*tmp++);
    if(link<0) { /* unlinked */
      /* Then do pat_hom / mat_het configs */
      z*=.5;
      z1=z*(qval[0]+qval[1]);
      peel_famval[1]*=z1;
      peel_famval[2]*=z1;
      z1=z*(qval[2]+qval[3]);
      peel_famval[13]*=z1;
      peel_famval[14]*=z1;
      /* Then do pat_het / mat_hom configs */
      z1=z*(qval[0]+qval[2]);
      peel_famval[4]*=z1;
      peel_famval[8]*=z1;
      z1=z*(qval[1]+qval[3]);
      peel_famval[7]*=z1;
      peel_famval[11]*=z1;
      /* Then do double het configs */
      ix=ix3;
      for(i=0;i<4;i++) peel_famval[*(ix++)]*=.25;
    } else { /* linked */
      tp=id_array[kid].tp;
      /* Then do pat_hom / mat_het configs */
      tpp=id_array[kid].tpp[X_MAT];
      peel_famval[1]*=z*(tpp[X_MAT]*qval[1]+tpp[X_PAT]*qval[0]);
      peel_famval[2]*=z*(tpp[X_MAT]*qval[0]+tpp[X_PAT]*qval[1]);
      peel_famval[13]*=z*(tpp[X_MAT]*qval[3]+tpp[X_PAT]*qval[2]);
      peel_famval[14]*=z*(tpp[X_MAT]*qval[2]+tpp[X_PAT]*qval[3]);
      /* Then do pat_het / mat_hom configs */
      tpp=id_array[kid].tpp[X_PAT];
      peel_famval[4]*=z*(tpp[X_MAT]*qval[2]+tpp[X_PAT]*qval[0]);
      peel_famval[8]*=z*(tpp[X_MAT]*qval[0]+tpp[X_PAT]*qval[2]);
      peel_famval[7]*=z*(tpp[X_MAT]*qval[3]+tpp[X_PAT]*qval[1]);
      peel_famval[11]*=z*(tpp[X_MAT]*qval[1]+tpp[X_PAT]*qval[3]);
      /* Then do double het configs */
      peel_famval[5]*=z*(tp[X_MM_PM]*qval[3]+tp[X_MP_PM]*qval[2]+tp[X_MM_PP]*qval[1]+tp[X_MP_PP]*qval[0]);
      peel_famval[6]*=z*(tp[X_MM_PM]*qval[2]+tp[X_MP_PM]*qval[3]+tp[X_MM_PP]*qval[0]+tp[X_MP_PP]*qval[1]);
      peel_famval[9]*=z*(tp[X_MM_PM]*qval[1]+tp[X_MP_PM]*qval[0]+tp[X_MM_PP]*qval[3]+tp[X_MP_PP]*qval[2]);
      peel_famval[10]*=z*(tp[X_MM_PM]*qval[0]+tp[X_MP_PM]*qval[1]+tp[X_MM_PP]*qval[2]+tp[X_MP_PP]*qval[3]);
    }
  }
  if(pivot== -2)	{
    p1=0.0;
    tmp=peel_famval;
    for(n=0;n<16;n++) p1+=*(tmp++);
#ifdef DEBUG
    if(isnan(p1)) ABT_FUNC("Floating point error\n");
#endif
    if(p1<=0.0) {
      if(!(s_flag&1)) {
	return -DBL_MAX;
      }
      ABT_FUNC("Zero probability in peeling operation\n");
    }
    prob+=log(p1);
    k=element->out_index;
    rf[k].n_ind=4;
    get_rf_memory(rf+k,16,TRT_MBLOCK,loki);
    tmp=peel_famval;
    tmp1=rf[k].p;
    z=1.0/p1;
    for(n=0;n<16;n++) *(tmp1++)=*(tmp++)*z;
    return prob;
  }
  p1=0.0;
  tmp=peel_famval;
  if(pivot<0) {
    for(i=0;i<16;i++) p1+=*(tmp++);
#ifdef DEBUG
    if(p1<0.0 || isnan(p1)) ABT_FUNC("Internal error - zero prob\n");
#endif
    prob+=log(p1);
  } else {
    if(ids==pivot) {
      tmp1=qval;
      for(i=0;i<4;i++) {
	z=0.0;
	for(j=0;j<4;j++) z+=*(tmp++);
	*(tmp1++)=z;
	p1+=z;
      }
    } else if(idd==pivot) {
      tmp1=qval;
      ix=ix2;
      for(i=0;i<4;i++) {
	z=0.0;
	for(j=0;j<4;j++) z+=peel_famval[*(ix++)];
	*(tmp1++)=z;
	p1+=z;
      }
    } else {
      if((k=id_array[pivot].rfp)>=0) { /* Insert Previously computed R_Func */
	for(j=0;j<n_idx;j++) pval[j]=rf[k].p[j];
      } else for(j=0;j<n_idx;j++) pval[j]=1.0;
      ix=ix1;
      tmp1=qval;
      for(i=0;i<4;i++) *(tmp1++)=peel_famval[*(ix++)];
      if(link<0) {
	z=.5*(peel_famval[1]+peel_famval[2]);
	qval[0]+=z;
	qval[1]+=z;
	z=.5*(peel_famval[13]+peel_famval[14]);
	qval[2]+=z;
	qval[3]+=z;
	z=.5*(peel_famval[4]+peel_famval[8]);
	qval[0]+=z;
	qval[2]+=z;
	z=.5*(peel_famval[7]+peel_famval[11]);
	qval[1]+=z;
	qval[3]+=z;
	ix=ix3;
	for(i=0;i<4;i++) {
	  z=.25*peel_famval[*(ix++)];
	  tmp=qval;
	  for(j=0;j<4;j++) *(tmp++)+=z;
	}
      } else {
	tp=id_array[pivot].tp;
	tpp=id_array[pivot].tpp[X_MAT];
	z=peel_famval[1];
	qval[0]+=tpp[X_PAT]*z;
	qval[1]+=tpp[X_MAT]*z;
	z=peel_famval[2];
	qval[1]+=tpp[X_PAT]*z;
	qval[0]+=tpp[X_MAT]*z;
	z=peel_famval[13];
	qval[2]+=tpp[X_PAT]*z;
	qval[3]+=tpp[X_MAT]*z;
	z=peel_famval[14];
	qval[3]+=tpp[X_PAT]*z;
	qval[2]+=tpp[X_MAT]*z;
	tpp=id_array[pivot].tpp[X_PAT];
	z=peel_famval[4];
	qval[0]+=tpp[X_PAT]*z;
	qval[2]+=tpp[X_MAT]*z;
	z=peel_famval[8];
	qval[2]+=tpp[X_PAT]*z;
	qval[0]+=tpp[X_MAT]*z;
	z=peel_famval[7];
	qval[1]+=tpp[X_PAT]*z;
	qval[3]+=tpp[X_MAT]*z;
	z=peel_famval[11];
	qval[3]+=tpp[X_PAT]*z;
	qval[1]+=tpp[X_MAT]*z;
	z=peel_famval[5];
	qval[3]+=tp[X_MM_PM]*z;
	qval[2]+=tp[X_MP_PM]*z;
	qval[1]+=tp[X_MM_PP]*z;
	qval[0]+=tp[X_MP_PP]*z;
	z=peel_famval[6];
	qval[2]+=tp[X_MM_PM]*z;
	qval[3]+=tp[X_MP_PM]*z;
	qval[0]+=tp[X_MM_PP]*z;
	qval[1]+=tp[X_MP_PP]*z;
	z=peel_famval[9];
	qval[1]+=tp[X_MM_PM]*z;
	qval[0]+=tp[X_MP_PM]*z;
	qval[3]+=tp[X_MM_PP]*z;
	qval[2]+=tp[X_MP_PP]*z;
	z=peel_famval[10];
	qval[0]+=tp[X_MM_PM]*z;
	qval[1]+=tp[X_MP_PM]*z;
	qval[2]+=tp[X_MM_PP]*z;
	qval[3]+=tp[X_MP_PP]*z;
      }
      for(j=0;j<n_idx;j++) p1+=(qval[j]*=pval[j]);
    }
#ifdef DEBUG
    if(isnan(p1)) ABT_FUNC("Floating point error\n");
#endif
    if(p1<=0.0) {
      if(!(s_flag&1)) {
	return -DBL_MAX;
      }
      ABT_FUNC("Zero probability in peeling operation\n");
    }
    prob+=log(p1);
    k=element->out_index;
    id_array[pivot].rfp=k;
    rf[k].n_ind=2;
    get_rf_memory(rf+k,n_idx,TRT_MBLOCK,loki);
    for(j=0;j<n_idx;j++) rf[k].p[j]=qval[j]/p1;
  }
  return prob;
}

