#include <config.h>
#include <stdlib.h>
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_simple_peel.h"

/* This is when all (non-pruned) children in a family have completely determined genotypes
 * so we can peel to the 2 parents separately */
double peel_to_par(const struct Simple_Element *element,const int locus,pen_func pen,
		   lk_ulong **a_set,struct R_Func *rf,struct loki *loki)
{
  int ids,idd,i,j,k,l,m,n,fsp=0,fsp1=0,n_off,*off,kid,*all,n_all,n_bits,n_idx,comp;
  double prob=0.0,*tpp,p1,pp[2],qp,*qval,*mval,*pval;
  lk_ulong a,cm[2];
  struct fset *peel_fs;
  struct Id_Record *id_array;
  struct peel_mem *work;
  struct Marker *mark;
	
  prt_peel_trace(TRACE_LEVEL_1,"In %s(%p,%d,%p)\n",__func__,(void *)element,locus,pen);
  id_array=loki->pedigree->id_array;
  work=&loki->peel->workspace;
  ids=element->sire-1;
  idd=element->dam-1;
  comp=id_array[ids].comp;
  mark=loki->markers->marker+locus;
  if(mark->ngens[ids]>mark->ngens[idd]) i=mark->ngens[ids];
  else i=mark->ngens[idd];
  off=element->off;
  n_off=element->n_off;
  n_all=mark->n_all1[comp];
  n_bits=num_bits(n_all);
  n_idx=1<<(n_bits+n_bits);
  qval=work->s2;
  pval=qval+n_idx;
  mval=pval+n_idx;
  peel_fs=work->s0;
  /* Construct set of possible parental genotypes */
  for(j=0;j<n_idx;j++) mval[j]=pval[j]=0.0;
  if((k=id_array[ids].rfp)>=0) { /* Insert Previously computed R_Func */
    a=(1L<<n_bits)-1;
    for(j=0;j<rf[k].n_terms;j++) {
      i=(int)rf[k].index[j];
      peel_fs[fsp].pat_gene[X_MAT]=i&a;
      peel_fs[fsp++].pat_gene[X_PAT]=i>>n_bits;
      pval[i]=rf[k].p[j];
    }
  } else for(i=0;i<n_all;i++) {
    a=a_set[ids][i];
    if(a) for(j=0;j<n_all;j++) if(a&(1L<<j)) {
      peel_fs[fsp].pat_gene[X_MAT]=i;
      peel_fs[fsp++].pat_gene[X_PAT]=j;
      pval[(j<<n_bits)|i]=1.0;
    }
  }
  if((k=id_array[idd].rfp)>=0) { /* Insert Previously computed R_Func */
    a=(1L<<n_bits)-1;
    for(j=0;j<rf[k].n_terms;j++) {
      i=(int)rf[k].index[j];
      peel_fs[fsp1].mat_gene[X_MAT]=i&a;
      peel_fs[fsp1++].mat_gene[X_PAT]=i>>n_bits;
      mval[i]=rf[k].p[j];
    }
  } else for(i=0;i<n_all;i++) {
    a=a_set[idd][i];
    if(a) for(j=0;j<n_all;j++) if(a&(1L<<j)) {
      peel_fs[fsp1].mat_gene[X_MAT]=i;
      peel_fs[fsp1++].mat_gene[X_PAT]=j;
      mval[(j<<n_bits)|i]=1.0;
    }
  }
  for(m=0;m<n_off;m++)	{
    kid=abs(off[m])-1;
    if(id_array[kid].rfp>=0)
      ABT_FUNC("Internal error - no R-Functions expected for this operation\n");
    if(mark->ngens[kid]>1)
      ABT_FUNC("Internal error - offspring genotypes should be determined for this operation\n");
    all=id_array[kid].allele;
    for(k=0;k<n_idx;k++) qval[k]=1.0;
    if(pen && off[m]>0) pen(qval,kid,&mark->locus,n_all,n_bits,loki);
    qp=qval[((all[X_PAT]-1)<<n_bits)|(all[X_MAT]-1)];
    cm[0]=mark->req_set[0][kid];
    cm[1]=mark->req_set[1][kid];
    tpp=id_array[kid].tpp[X_PAT];
    for(n=0;n<fsp;n++) {
      i=peel_fs[n].pat_gene[X_MAT];
      j=peel_fs[n].pat_gene[X_PAT];
      if((1L<<i)&cm[X_PAT]) i=n_all-1;
      if((1L<<j)&cm[X_PAT]) j=n_all-1;
      p1=0.0;
      if(all[X_PAT]==(i+1)) p1=tpp[X_MAT];
      if(all[X_PAT]==(j+1)) p1+=tpp[X_PAT];
      pval[(j<<n_bits)|i]*=p1*qp;
    }
    tpp=id_array[kid].tpp[X_MAT];
    for(n=0;n<fsp1;n++) {
      i=peel_fs[n].mat_gene[X_MAT];
      j=peel_fs[n].mat_gene[X_PAT];
      if((1L<<i)&cm[X_MAT]) i=n_all-1;
      if((1L<<j)&cm[X_MAT]) j=n_all-1;
      p1=0.0;
      if(all[X_MAT]==(i+1)) p1=tpp[X_MAT];
      if(all[X_MAT]==(j+1)) p1+=tpp[X_PAT];
      mval[(j<<n_bits)|i]*=p1*qp;
    }
  }
  pp[X_MAT]=pp[X_PAT]=0.0;
  for(i=0;i<n_idx;i++) {
    pp[X_MAT]+=mval[i];
    pp[X_PAT]+=pval[i];
  }
  prob+=log(pp[X_MAT]*pp[X_PAT]);
  k=element->out_index;
  rf[k].n_ind=2;
  rf[k].n_terms=fsp;
  get_rf_memory(rf+k,fsp,MRK_MBLOCK,loki);
  get_rf_memory(rf+k+1,fsp1,MRK_MBLOCK,loki);
  for(n=0;n<fsp;n++) {
    i=peel_fs[n].pat_gene[X_MAT];
    j=peel_fs[n].pat_gene[X_PAT];
    l=(j<<n_bits)|i;
    rf[k].index[n]=(lk_ulong)l;
    rf[k].p[n]=pval[l]>0.0?pval[l]/pp[X_PAT]:0.0;
  }
  id_array[ids].rfp=k++;
  rf[k].n_ind=2;
  rf[k].n_terms=fsp1;
  for(n=0;n<fsp1;n++) {
    i=peel_fs[n].mat_gene[X_MAT];
    j=peel_fs[n].mat_gene[X_PAT];
    l=(j<<n_bits)|i;
    rf[k].index[n]=(lk_ulong)l;
    rf[k].p[n]=mval[l]>0.0?mval[l]/pp[X_MAT]:0.0;
  }
  id_array[idd].rfp=k;
  return prob;
}

