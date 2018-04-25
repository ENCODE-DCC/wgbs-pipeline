/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                  Simon Heath - MSKCC                                     *
 *                                                                          *
 *                       July 2001                                          *
 *                                                                          *
 * pass_founder_genes.c:                                                    *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "seg_pen.h"

static struct cg_stack *stack;
static int stack_size=256;

void pass_founder_genes(struct Locus *loc,const struct loki *loki)
{
  int i,i1,j,s,**genes,**seg,comp,cs;
  struct Id_Record *id_array;
	
  genes=loc->genes;
  seg=loc->seg;
  id_array=loki->pedigree->id_array;
  for(i=comp=0;comp<loki->pedigree->n_comp;comp++) { 
    cs=loki->pedigree->comp_size[comp];
    if(!(loki->pedigree->singleton_flag[comp])) {
      for(i1=j=0;i1<cs;i1++,i++) {
	s=seg[X_PAT][i];
#ifdef DEBUG
	if(s<-1) ABT_FUNC("Shouldn't happen!\n");
#endif
	if(s>=0)	{
	  genes[X_PAT][i]=genes[s][id_array[i].sire-1];
	} else {
#ifdef DEBUG
	  if(id_array[i].sire) {
	    ABT_FUNC("Internal error - not a founder\n");
	  }
#endif
	  genes[X_PAT][i]= ++j;
	}
	s=seg[X_MAT][i];
#ifdef DEBUG
	if(s<-1) ABT_FUNC("Shouldn't happen!\n");
#endif
	if(s>=0) {
	  genes[X_MAT][i]=genes[s][id_array[i].dam-1];
	} else {
#ifdef DEBUG
	  if(id_array[i].dam) {
	    ABT_FUNC("Internal error - not a founder\n");
	  }
#endif
	  genes[X_MAT][i]= ++j;
	}
      }
#ifdef DEBUG
      if(j>loki->pedigree->comp_ngenes[comp]) {
	ABT_FUNC("Mismatch in founder gene number\n");
      }
#endif
    } else {
      for(i1=j=0;i1<cs;i1++,i++) {
	genes[X_MAT][i]=++j;
	genes[X_PAT][i]=++j;
      }
    }
  }
}

#if 0
static void check_founder_genes(struct Locus *loc,const struct loki *loki)
{
  int i,i1,j,s,**genes,**seg,comp,cs;
  struct Id_Record *id_array;
	
  genes=loc->genes;
  seg=loc->seg;
  id_array=loki->pedigree->id_array;
  for(i=comp=0;comp<loki->pedigree->n_comp;comp++) { 
    cs=loki->pedigree->comp_size[comp];
    if(!(loki->pedigree->singleton_flag[comp])) {
      for(i1=j=0;i1<cs;i1++,i++) {
	s=seg[X_PAT][i];
	if(s>=0)	{
	  if(genes[X_PAT][i]!=genes[s][id_array[i].sire-1]) {
	    printf("Err at %s:%d, genes mismatch.  comp=%d, i=%d\n",__FILE__,__LINE__,comp,i);
	  }
	} else {
	  if(genes[X_PAT][i]!=++j) {
	    printf("Err at %s:%d, genes mismatch.  comp=%d, i=%d\n",__FILE__,__LINE__,comp,i);
	  }
	}
	s=seg[X_MAT][i];
	if(s>=0) {
	  if(genes[X_MAT][i]!=genes[s][id_array[i].dam-1]) {
	    printf("Err at %s:%d, genes mismatch.  comp=%d, i=%d\n",__FILE__,__LINE__,comp,i);
	  }
	} else {
	  if(genes[X_MAT][i]!=++j) {
	    printf("Err at %s:%d, genes mismatch.  comp=%d, i=%d\n",__FILE__,__LINE__,comp,i);
	  }
	}
      }
#ifdef DEBUG
      if(j>loki->pedigree->comp_ngenes[comp]) {
	ABT_FUNC("Mismatch in founder gene number\n");
      }
#endif
    } else {
      for(i1=j=0;i1<cs;i1++,i++) {
	genes[X_MAT][i]=++j;
	genes[X_PAT][i]=++j;
      }
    }
  }
}
#endif

void pass_founder_genes_alloc(void) {
  if(!stack) {
    if(!(stack=malloc(sizeof(struct cg_stack)*stack_size))) ABT_FUNC(MMsg);
  }
}

void pass_founder_genes_dealloc(void) {
  if(stack) free(stack);
  stack=0;
}

/* Change genes descended from the par_flag gene of individual i */
int pass_founder_genes1(const int locus,int i,int par_flag,const struct loki *loki)
{
  int **genes,par,fg=0,**seg;
  int g,j=0,k,k1,kid,*hap,stack_ptr=0,locus_type;
  struct Locus *loc;
  struct Marker *mark;
  struct Id_Record *id_array,**kids;
	
  id_array=loki->pedigree->id_array;
  par=(par_flag==X_MAT)?id_array[i].dam-1:id_array[i].sire-1;
  if(locus>=0) {
    mark=loki->markers->marker+locus;
    loc=&mark->locus;
    hap=mark->haplo;
    if(mark->mterm && mark->mterm[0]) locus_type=1;
    else locus_type=0;
    genes=loc->genes;
    seg=loc->seg;
    g=genes[par_flag][i];
    genes[par_flag][i]=genes[seg[par_flag][i]][par];
    if(genes[par_flag][i]!=g) {
      if(hap[i] || (locus_type && id_array[i].res[0])) fg=1;
      if(id_array[i].nkids) {
	for(;;) {
	  k=id_array[i].sex==1?X_PAT:X_MAT;
	  kids=id_array[i].kids;
	  k1=0;
	  for(;j<id_array[i].nkids;j++) {
	    kid=kids[j]->idx;
	    if(seg[k][kid]==par_flag) {
	      g=genes[k][kid];
	      genes[k][kid]=genes[par_flag][i];
	      if(g!=genes[k][kid]) {
		if(hap[kid] || (locus_type && id_array[kid].res[0])) fg=1;
		if(id_array[kid].nkids) {
		  stack[stack_ptr].id=i;
		  stack[stack_ptr].par_flag=par_flag;
		  stack[stack_ptr++].kid_ptr=j+1;
		  i=kid;
		  j=0;
		  par_flag=k;
		  k1=1;
		  break;
		}
	      }
	    }
	  }
	  if(!k1) {
	    if(stack_ptr) {
	      j=stack[--stack_ptr].kid_ptr;
	      par_flag=stack[stack_ptr].par_flag;
	      i=stack[stack_ptr].id;
	    } else break;
	  }
	} 
      }
    }
  } else {
    loc=&loki->models->tlocus[-1-locus];
    genes=loc->genes;
    seg=loc->seg;
    g=genes[par_flag][i];
    genes[par_flag][i]=genes[seg[par_flag][i]][par];
    if(genes[par_flag][i]!=g) {
      if(id_array[i].res[0]) fg=1;
      if(id_array[i].nkids) {
	kids=id_array[i].kids;
	for(;;) {
	  k=id_array[i].sex==1?X_PAT:X_MAT;
	  k1=0;
	  for(;j<id_array[i].nkids;j++) {
	    kid=kids[j]->idx;
	    if(seg[k][kid]==par_flag) {
	      g=genes[k][kid];
	      genes[k][kid]=genes[par_flag][i];
	      if(g!=genes[k][kid]) {
		if(id_array[kid].res[0]) fg=1;
		if(id_array[kid].nkids) {
		  stack[stack_ptr].id=i;
		  stack[stack_ptr].par_flag=par_flag;
		  stack[stack_ptr++].kid_ptr=j+1;
		  i=kid;
		  j=0;
		  par_flag=k;
		  k1=1;
		  break;
		}
	      }
	    }
	  }
	  if(!k1) {
	    if(stack_ptr) {
	      j=stack[--stack_ptr].kid_ptr;
	      par_flag=stack[stack_ptr].par_flag;
	      i=stack[stack_ptr].id;
	    } else break;
	  }
	} 
      }
    }
  }
  return fg;
}

/* Change the genes descended from the par_flag genes of the nkids individuals in kids */
int pass_founder_genes1a(const int locus,const int *kids,const int nkids,int par_flag,const struct loki *loki)
{
  int **genes,fg=0,**seg;
  int i,g,j,j2,k,k1,kid,par,*hap,stack_ptr=0,locus_type;
  struct Locus *loc;
  struct Marker *mark;
  struct Id_Record *id_array,**kds;
	
  id_array=loki->pedigree->id_array;
  i=kids[0];
  par=(par_flag==X_MAT)?id_array[i].dam-1:id_array[i].sire-1;
  if(locus>=0) {
    mark=loki->markers->marker+locus;
    loc=&mark->locus;
    if(mark->mterm && mark->mterm[0]) locus_type=1;
    else locus_type=0;
    genes=loc->genes;
    seg=loc->seg;
    hap=mark->haplo;
    for(j2=0;j2<nkids;j2++) {
      i=kids[j2];
      g=genes[par_flag][i];
      genes[par_flag][i]=genes[seg[par_flag][i]][par];
      if(genes[par_flag][i]!=g) {
	if(hap[i] || (locus_type && id_array[i].res[0])) fg=1;
	if(id_array[i].nkids) {
	  j=0;
	  for(;;) {
	    k=id_array[i].sex==1?X_PAT:X_MAT;
	    kds=id_array[i].kids;
	    k1=0;
	    for(;j<id_array[i].nkids;j++) {
	      kid=kds[j]->idx;
	      if(seg[k][kid]==par_flag) {
		g=genes[k][kid];
		genes[k][kid]=genes[par_flag][i];
		if(g!=genes[k][kid]) {
		  if(hap[kid] || (locus_type && id_array[kid].res[0])) fg=1;
		  if(id_array[kid].nkids) {
		    stack[stack_ptr].id=i;
		    stack[stack_ptr].par_flag=par_flag;
		    stack[stack_ptr++].kid_ptr=j+1;
		    i=kid;
		    j=0;
		    par_flag=k;
		    k1=1;
		    break;
		  }
		}
	      }
	    }
	    if(!k1) {
	      if(stack_ptr) {
		j=stack[--stack_ptr].kid_ptr;
		par_flag=stack[stack_ptr].par_flag;
		i=stack[stack_ptr].id;
	      } else break;
	    }
	  } 
	}
      }
    }
  } else {
    loc=&loki->models->tlocus[-1-locus];
    genes=loc->genes;
    seg=loc->seg;
    for(j2=0;j2<nkids;j2++) {
      i=kids[j2];
      g=genes[par_flag][i];
      genes[par_flag][i]=genes[seg[par_flag][i]][par];
      if(genes[par_flag][i]!=g) {
	if(id_array[i].res[0]) fg=1;
	if(id_array[i].nkids) {
	  j=0;
	  for(;;) {
	    k=id_array[i].sex==1?X_PAT:X_MAT;
	    kds=id_array[i].kids;
	    k1=0;
	    for(;j<id_array[i].nkids;j++) {
	      kid=kds[j]->idx;
	      if(seg[k][kid]==par_flag) {
		g=genes[k][kid];
		genes[k][kid]=genes[par_flag][i];
		if(g!=genes[k][kid]) {
		  if(id_array[kid].res[0]) fg=1;
		  if(id_array[kid].nkids) {
		    stack[stack_ptr].id=i;
		    stack[stack_ptr].par_flag=par_flag;
		    stack[stack_ptr++].kid_ptr=j+1;
		    i=kid;
		    j=0;
		    par_flag=k;
		    k1=1;
		    break;
		  }
		}
	      }
	    }
	    if(!k1) {
	      if(stack_ptr) {
		j=stack[--stack_ptr].kid_ptr;
		par_flag=stack[stack_ptr].par_flag;
		i=stack[stack_ptr].id;
	      } else break;
	    }
	  } 
	}
      }
    }
  }
  return fg;
}

/* Change the genes descended from both genes of the nkids individuals in kids
 * Note that we always go to the grandkids even if the kids genes have not changed */
int pass_founder_genes1b(const int locus,const int *kids,const int nkids,const struct loki *loki)
{
  int **genes,fg=0,**seg;
  int i,g,j,j2,k,k1,kid,par_flag,*hap,par[2],stack_ptr=0,locus_type;
  struct Locus *loc;
  struct Marker *mark;
  struct Id_Record *id_array,**kds;

  id_array=loki->pedigree->id_array;
  i=kids[0];
  par[X_PAT]=id_array[i].sire-1;
  par[X_MAT]=id_array[i].dam-1;
  if(locus>=0) {
    mark=loki->markers->marker+locus;
    loc=&mark->locus;
    if(mark->mterm && mark->mterm[0]) locus_type=1;
    else locus_type=0;
    genes=loc->genes;
    seg=loc->seg;
    hap=mark->haplo;
    for(j2=0;j2<nkids;j2++) {
      i=kids[j2];
      for(par_flag=0;par_flag<2;par_flag++) {
	g=genes[par_flag][i];
	genes[par_flag][i]=genes[seg[par_flag][i]][par[par_flag]];
	if(genes[par_flag][i]!=g && (hap[i] || (locus_type && id_array[i].res[0]))) fg=1;
	if(id_array[i].nkids) {
	  j=0;
	  for(;;) {
	    kds=id_array[i].kids;
	    k=id_array[i].sex==1?X_PAT:X_MAT;
	    k1=0;
	    for(;j<id_array[i].nkids;j++) {
	      kid=kds[j]->idx;
	      if(seg[k][kid]==par_flag) {
		g=genes[k][kid];
		genes[k][kid]=genes[par_flag][i];
		if(g!=genes[k][kid]) {
		  if(hap[kid] || (locus_type && id_array[kid].res[0])) fg=1;
		  if(id_array[kid].nkids) {
		    stack[stack_ptr].id=i;
		    stack[stack_ptr].par_flag=par_flag;
		    stack[stack_ptr++].kid_ptr=j+1;
		    i=kid;
		    j=0;
		    par_flag=k;
		    k1=1;
		    break;
		  }
		}
	      }
	    }
	    if(!k1) {
	      if(stack_ptr) {
		j=stack[--stack_ptr].kid_ptr;
		par_flag=stack[stack_ptr].par_flag;
		i=stack[stack_ptr].id;
	      } else break;
	    }
	  } 
	}
      }
    }
  } else {
    loc=&loki->models->tlocus[-1-locus];
    genes=loc->genes;
    seg=loc->seg;
    for(j2=0;j2<nkids;j2++) {
      i=kids[j2];
      for(par_flag=0;par_flag<2;par_flag++) {
	g=genes[par_flag][i];
	genes[par_flag][i]=genes[seg[par_flag][i]][par[par_flag]];
	if(genes[par_flag][i]!=g && id_array[i].res[0]) fg=1;
	if(id_array[i].nkids) {
	  j=0;
	  for(;;) {
	    kds=id_array[i].kids;
	    k=id_array[i].sex==1?X_PAT:X_MAT;
	    k1=0;
	    for(;j<id_array[i].nkids;j++) {
	      kid=kds[j]->idx;
	      if(seg[k][kid]==par_flag) {
		g=genes[k][kid];
		genes[k][kid]=genes[par_flag][i];
		if(g!=genes[k][kid]) {
		  if(id_array[kid].res[0]) fg=1;
		  if(id_array[kid].nkids) {
		    stack[stack_ptr].id=i;
		    stack[stack_ptr].par_flag=par_flag;
		    stack[stack_ptr++].kid_ptr=j+1;
		    i=kid;
		    j=0;
		    par_flag=k;
		    k1=1;
		    break;
		  }
		}
	      }
	    }
	    if(!k1) {
	      if(stack_ptr) {
		j=stack[--stack_ptr].kid_ptr;
		par_flag=stack[stack_ptr].par_flag;
		i=stack[stack_ptr].id;
	      } else break;
	    }
	  } 
	}
      }
    }
  }
  return fg;
}

int pass_founder_genes2(const int locus,const int comp,int **seg,const struct loki *loki)
{
  int i,i1,j,k,s,**genes,g[2],fg=0,*hap;
  struct Id_Record *id_array;
  struct Marker *mark;
	
  mark=loki->markers->marker+locus;
  if(mark->n_all1[comp]<2) return fg; 
  id_array=loki->pedigree->id_array;
  genes=mark->locus.genes;
  hap=mark->haplo;
  i=loki->pedigree->comp_start[comp];
  for(i1=j=0;i1<loki->pedigree->comp_size[comp];i1++,i++) {
    s=seg[X_PAT][i];
    if(!fg && hap[i]) {
      g[X_PAT]=genes[X_PAT][i];
      g[X_MAT]=genes[X_MAT][i];
      k=1;
    } else k=0;
    if(s>=0)	{
      genes[X_PAT][i]=genes[s][id_array[i].sire-1];
    } else genes[X_PAT][i]= ++j;
    s=seg[X_MAT][i];
    if(s>=0) {
      genes[X_MAT][i]=genes[s][id_array[i].dam-1];
    } else genes[X_MAT][i]= ++j;
    if(k) {
      for(s=0;s<2;s++) if(genes[s][i]!=g[s]) {
	fg=1;
	break;
      }
    }
  }
  return fg;
}
