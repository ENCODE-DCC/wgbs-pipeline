/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - MSKCC                                   *
 *                                                                          *
 *                        July 1997                                         *
 *                                                                          *
 * min_deg.c                                                                *
 *                                                                          *
 * Produce a min-degree ordering for factoring of a sparse matrix           *
 * components on a marker by marker basis.                                  *
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
#include <math.h>
#include <assert.h>

#include "utils.h"
#include "loki_struct.h"
#include "lk_malloc.h"
#include "min_deg.h"

#define BSIZE 1024

static struct mat_elem *free_node_list,**mat;
static struct remember *MinDegBlock,*FirstRemBlock;
static double *score;
static int *degree,*fill_in;

static int calc_degree(int i,int *inv,int *flag,int *nn,struct deg_list *dlist,int fg)
{
  int j,k,x,x1,n=0,n1;
  struct mat_elem *p,*p1;
  struct deg_list *deg,*deg1;
	
  p=mat[i];
  if(p) {
    inv[n++]=i;
    while(p) {
      inv[n++]=p->x;
      p=p->next;
    }
    for(j=1;j<n;j++) {
      x=inv[j];
      if(flag[x]) continue;
      p1=mat[x];
      k=1;
      while(p1) {
	x1=p1->x;
	if(k==n) {
	  if(x1==i) p1=p1->next;
	  break;
	}
	if(inv[k]==x) k++;
	if(x1!=i && x1!=inv[k++]) break;
	p1=p1->next;
	if(!p1 && k<n && inv[k]==x) k++;
      }
      if(!p1 && k==n) {
	flag[x]=3;
	deg=dlist+x;
	deg1=deg;
	while(deg1->abs_list) deg1=deg1->abs_list;
	deg1->abs_list=dlist[i].abs_list;
	dlist[i].abs_list=deg;
      }
    }
  }
  n1=0;
  for(k=(fg&1)?0:1;k<n;k++) {
    x=inv[k];
    if(!flag[x]) {
      n1++;
      deg=dlist[x].abs_list;
      while(deg) {
	n1++;
	deg=deg->abs_list;
      }
    }
  }
  *nn=n;
  return n1;;
}

static int calc_score(int i,double *wt,struct pair_wt *wt1,int *inv,int *inv1,int *flag,int *nn,double *sc,struct deg_list *alist)
{
  int j,k,k1,k2,k3,x,x1,n=0,n1;
  struct mat_elem *p,*p1;
  struct deg_list *deg,*deg1;
  double w=0.0;
	
  p=mat[i];
  if(p) {
    inv[n++]=i;
    while(p) {
      inv[n++]=p->x;
      p=p->next;
    }
    for(j=1;j<n;j++) {
      x=inv[j];
      if(flag[x]) continue;
      p1=mat[x];
      k=1;
      while(p1) {
	x1=p1->x;
	if(k==n) {
	  if(x1==i) p1=p1->next;
	  break;
	}
	if(inv[k]==x) k++;
	if(x1!=i && x1!=inv[k++]) break;
	p1=p1->next;
	if(!p1 && k<n && inv[k]==x) k++;
      }
      if(!p1 && k==n) {
	flag[x]=3;
	deg=alist+x;
	deg1=deg;
	while(deg1->abs_list) deg1=deg1->abs_list;
	deg1->abs_list=alist[i].abs_list;
	alist[i].abs_list=deg;
      }
    }
  }
  n1=0;
  w=wt[i];
  deg=alist[i].abs_list;
  while(deg) {
    w+=wt[deg->node];
    deg=deg->abs_list;
  }
  for(k=1;k<n;k++) {
    x=inv[k];
    if(!flag[x]) {
      inv1[n1++]=x;
      deg=alist[x].abs_list;
      while(deg) {
	inv1[n1++]=deg->node;
	deg=deg->abs_list;
      }
    }
  }
  w=0.0;
  for(k=0;k<n1;k++) {
    k1=inv1[k];
    k2=wt1?wt1[k1].pair_node:-1;
    if(k2>=0) {
      for(k3=0;k3<n1;k3++) if(inv1[k3]==k2) break;
      if(k3<n1) {
	if(k1<k2) w+=wt1[k1].wt;
      } else w+=wt[k1];
    } else w+=wt[k1];
  }
  *nn=n;
  *sc=w;
  return n1;
}

static void delete_nd(struct deg_list *p,struct deg_list **first,int deg)
{
  if(p->last) p->last->next=p->next;
  else first[deg]=p->next;
  if(p->next) p->next->last=p->last;
}

static void add_nd(struct deg_list *p,struct deg_list **first,int deg)
{
  p->next=first[deg];
  if(p->next) p->next->last=p;
  p->last=0;
  first[deg]=p;
}

static void delete_node(int x,int y)
{
  struct mat_elem **p,*p1;
	
  p=mat+y;
  while(*p) {
    if((*p)->x==x) break;
    p=&(*p)->next;
  }
  if(*p) {
    p1=*p;
    *p=p1->next;
    p1->next=free_node_list;
    free_node_list=p1;
  }
}

static void delete_row_col(int y)
{
  struct mat_elem *p,*p1;
	
  p=p1=mat[y];
  while(p) {
    delete_node(y,p->x);
    p1=p;
    p=p->next;
  }
  if(p1) {
    p1->next=free_node_list;
    free_node_list=mat[y];
  }
  mat[y]=0;
}

static struct mat_elem *alloc_new_nodes(int n)
{
  struct mat_elem *p;
  int i;
	
  assert(n);
  p=lk_malloc(sizeof(struct mat_elem)*n);
  MinDegBlock=AddRemem(p,MinDegBlock);
  /* Link blocks together */
  for(i=0;i<n-1;i++) p[i].next=p+i+1;
  p[n-1].next=0;
  return p;
}

static struct mat_elem *new_node(void) 
{
  struct mat_elem *p;

  if(!free_node_list) free_node_list=alloc_new_nodes(BSIZE);
  p=free_node_list;
  free_node_list=p->next;
  return p;
}
 
static void add_node(int x,int y)
{
  struct  mat_elem **p,*p1;
	
  p=mat+y;
  while(*p) {
    if((*p)->x>=x) break;
    p=&(*p)->next;
  }
  p1=*p;
  if(!p1 || p1->x!=x) {
    *p=new_node();
    (*p)->next=p1;
    (*p)->x=x;
  }
}

int *min_deg(int n,int *mm,int *order,int *group,int fg)
{
  int i,j,k,k1,k2,k3,k4,k5,n_nodes;
  int *deg,*flag,*inv,*inv1,ct=0;
  struct deg_list *dlist,**first,*nd,*nd1,*nd2,*pp;
  struct mat_elem *p;
	
  if(n<1 || !mm) {
    if(FirstRemBlock) {
      FreeRemem(FirstRemBlock);
      FirstRemBlock=0;
    }
    return 0;
  }
  if(!FirstRemBlock) {
    FirstRemBlock=lk_malloc(sizeof(struct remember));
    MinDegBlock=FirstRemBlock;
    MinDegBlock->pos=0;
    MinDegBlock->next=0;
    free_node_list=0;
  }
  mat=lk_malloc(sizeof(void *)*n);
  for(i=0;i<n;i++) mat[i]=0;
  j=mm[0];
  for(i=0;i<n;i++) {
    k=mm[i+1];
    for(;j<k;j++) {
      add_node(i,mm[j]);
      add_node(mm[j],i);
    }
  }
  deg=lk_malloc(sizeof(int)*n*4);
  flag=deg+n;
  inv=flag+n;
  inv1=inv+n;
  dlist=lk_malloc(sizeof(struct deg_list)*n);
  first=lk_malloc(sizeof(void *)*n);
  for(i=0;i<n;i++) {
    first[i]=0;
    dlist[i].node=i;
    dlist[i].abs_list=0;
    flag[i]=0;
  }
  for(i=0;i<n;i++) if(!flag[i]) {
      k=calc_degree(i,inv,flag,&j,dlist,fg);
      deg[i]=k;
      for(k1=1;k1<j;k1++) {
	k2=inv[k1];
	if(flag[k2]==3) {
	  flag[k2]=4;
	  delete_row_col(k2);
	}
      }
    }
  n_nodes=0;
  for(i=n-1;i>=0;i--) if(!flag[i]) {
      k=deg[i];
      add_nd(dlist+i,first,k);
      n_nodes++;
    }
  while(n_nodes) {
    for(i=0;i<n;i++) if(first[i]) break;
    assert(i<n);
    nd=first[i];
    nd2=nd->abs_list;
    k2=0;
    while(nd2) {
      k2++;
      nd2=nd2->abs_list;
    }
    nd1=nd->next;
    while(nd1) {
      k3=0;
      nd2=nd1->abs_list;
      while(nd2) {
	k3++;
	nd2=nd2->abs_list;
      }
      if(((fg&2)?k3-k2:k2-k3)>0) {
	nd=nd1;
	k2=k3;
      }
      nd1=nd1->next;
    }
    j=nd->node;
    assert(!flag[j]);
    flag[j]=1;
    n_nodes--;
    if(group) group[ct]=1;
    order[ct++]=j;
    pp=nd->abs_list;
    if(pp) {
      while(pp) {
	k=pp->node;
	flag[k]=1;
	if(group) group[ct]=0;
	order[ct++]=k;
	pp=pp->abs_list;
      }
    }
    p=mat[j];
    k=0;
    while(p) {
      k1=p->x;
      assert(!flag[k1]);
      inv[k++]=k1;
      p=p->next;
    }
    delete_row_col(j);
    delete_nd(nd,first,i);
    for(k1=1;k1<k;k1++) {
      k3=inv[k1];
      for(k2=0;k2<k1;k2++) {
	k4=inv[k2];	
	add_node(k3,k4);
	add_node(k4,k3);
      }
    }
    for(k1=0;k1<k;k1++) {
      k2=inv[k1];
      if(!flag[k2]) {
	k3=calc_degree(k2,inv1,flag,&j,dlist,fg);
	for(k4=1;k4<j;k4++) {
	  k5=inv1[k4];
	  if(flag[k5]==3) {
	    flag[k5]=4;
	    n_nodes--;
	    delete_row_col(k5);
	    delete_nd(dlist+k5,first,deg[k5]);
	  }
	}
	k4=deg[k2];
	if(k3!=k4) {
	  nd=dlist+k2;
	  delete_nd(nd,first,k4);
	  add_nd(nd,first,k3);
	  deg[k2]=k3;
	}
      }
    }
  }
  free(first);
  free(dlist);
  free(deg);
  free(mat);
  mat=0;
  return order;
}

static int cmp_scores_0(const void *s1,const void *s2)
{
  double x1,x2;
  int i1,i2;
	
  x1=score[*((const int *)s1)];
  x2=score[*((const int *)s2)];
  if(x1<x2) return -1;
  if(x1>x2) return 1;
  i1=degree[*((const int *)s1)];
  i2=degree[*((const int *)s2)];
  if(i1<i2) return -1;
  if(i1>i2) return 1;
  i1=fill_in[*((const int *)s1)];
  i2=fill_in[*((const int *)s2)];
  if(i1<i2) return 1;
  if(i1>i2) return -1;
  return 0;
}

static int cmp_scores_1(const void *s1,const void *s2)
{
  double x1,x2;
  int i1,i2;
	
  i1=degree[*((const int *)s1)];
  i2=degree[*((const int *)s2)];
  if(i1<i2) return -1;
  if(i1>i2) return 1;
  x1=score[*((const int *)s1)];
  x2=score[*((const int *)s2)];
  if(x1<x2) return -1;
  if(x1>x2) return 1;
  i1=fill_in[*((const int *)s1)];
  i2=fill_in[*((const int *)s2)];
  if(i1<i2) return 1;
  if(i1>i2) return -1;
  return 0;
}

static int cmp_scores_2(const void *s1,const void *s2)
{
  double x1,x2;
  int i1,i2;
	
  i1=degree[*((const int *)s1)];
  i2=degree[*((const int *)s2)];
  if(i1<i2) return -1;
  if(i1>i2) return 1;
  i1=fill_in[*((const int *)s1)];
  i2=fill_in[*((const int *)s2)];
  if(i1<i2) return 1;
  if(i1>i2) return -1;
  x1=score[*((const int *)s1)];
  x2=score[*((const int *)s2)];
  if(x1<x2) return -1;
  if(x1>x2) return 1;
  return 0;
}

static int cmp_scores_3(const void *s1,const void *s2)
{
  double x1,x2;
  int i1,i2;
	
  i1=degree[*((const int *)s1)];
  i2=degree[*((const int *)s2)];
  if(i1<i2) return -1;
  if(i1>i2) return 1;
  i1=fill_in[*((const int *)s1)];
  i2=fill_in[*((const int *)s2)];
  if(i1<i2) return 1;
  if(i1>i2) return -1;
  x1=score[*((const int *)s1)];
  x2=score[*((const int *)s2)];
  if(x1<x2) return 1;
  if(x1>x2) return -1;
  return 0;
}

static int check_fill_in(int *inv,int *inv1,int *flag,int j,struct deg_list *alist)
{
  int k,k1,k2,k3,fill=0;
  struct mat_elem *p;
  struct deg_list *deg;
	
  for(k3=0,k=1;k<j;k++) {
    k2=inv[k];
    if(!flag[k2]) {
      inv1[k3++]=k2;
      deg=alist[k2].abs_list;
      while(deg) {
	inv1[k3++]=deg->node;
	deg=deg->abs_list;
      }
    }
  }
  for(k=0;k<k3-1;k++) {
    k2=inv1[k];
    p=mat[k2];
    while(p) {
      for(k1=k+1;k1<k3;k1++) if(p->x==inv1[k1]) break;
      if(k1==k3) fill++;
      p=p->next;
    }
  }
  return fill;
}

int *greedy(int n,int *mm,int *order,int *group,double *wt,struct pair_wt *wt1,int fg,double *cost)
{
	
  int i,j,k,k1,k2,k3,k4,k5,n_nodes;
  int *flag,*inv,*inv1,*inv2,*perm,ct=0;
  struct mat_elem *p;
  struct deg_list *alist,*deg;
  int (*cmp_func[])(const void *,const void *)={cmp_scores_0,cmp_scores_1,cmp_scores_2,cmp_scores_3};
  double z;
	
  if(n<1 || !mm) {
    if(FirstRemBlock) {
      FreeRemem(FirstRemBlock);
      FirstRemBlock=0;
    }
    return 0;
  }
  *cost=0;
  if(!FirstRemBlock) {
    FirstRemBlock=lk_malloc(sizeof(struct remember));
    MinDegBlock=FirstRemBlock;
    MinDegBlock->pos=0;
    MinDegBlock->next=0;
    free_node_list=0;
  }
  mat=lk_malloc(sizeof(void *)*n);
  for(i=0;i<n;i++) mat[i]=0;
  j=mm[0];
  for(i=0;i<n;i++) {
    k=mm[i+1];
    for(;j<k;j++) {
      add_node(i,mm[j]);
      add_node(mm[j],i);
    }
  }
  flag=lk_malloc(sizeof(int)*n*7);
  inv=flag+n;
  inv1=inv+n;
  inv2=inv1+n;
  perm=inv2+n;
  degree=perm+n;
  fill_in=degree+n;
  score=lk_malloc(sizeof(double)*n);
  alist=lk_malloc(sizeof(struct deg_list)*n);
  for(i=0;i<n;i++) {
    alist[i].node=i;
    alist[i].abs_list=0;
    flag[i]=0;
  }
  for(i=0;i<n;i++) if(!flag[i]) {
      degree[i]=calc_score(i,wt,wt1,inv,inv1,flag,&j,score+i,alist);
      for(k1=1;k1<j;k1++) {
	k2=inv[k1];
	if(flag[k2]==3) {
	  flag[k2]=4;
	  delete_row_col(k2);
	}
      }
      fill_in[i]=check_fill_in(inv,inv1,flag,j,alist);
      if(degree[i]==1 || !fill_in[i]) score[i]=-1.0;
    }
  n_nodes=0;
  for(i=0;i<n;i++) if(!flag[i]) perm[n_nodes++]=i;
  gnu_qsort(perm,(size_t)n_nodes,sizeof(int),cmp_func[fg]);
  while(n_nodes) {
    j=perm[0];
    if(group) group[ct]=1;
    k1=0;
    order[ct++]=j;
    inv1[k1++]=j;
    flag[j]=1;
    deg=alist[j].abs_list;
    while(deg) {
      k=deg->node;
      if(group) group[ct]=0;
      order[ct++]=k;
      inv1[k1++]=k;
      deg=deg->abs_list;
    }
    p=mat[j];
    k=0;
    while(p) {
      inv[k++]=p->x;
      inv1[k1++]=p->x;
      deg=alist[p->x].abs_list;
      while(deg) {
	inv1[k1++]=deg->node;
	deg=deg->abs_list;
      }
      p=p->next;
    }
    z=0.0;
    for(k2=0;k2<k1;k2++) {
      k3=inv1[k2];
      k4=wt1?wt1[k3].pair_node:-1;
      if(k4>=0) {
	for(k5=0;k5<k1;k5++) if(inv1[k5]==k4) break;
	if(k5<k1) {
	  if(k3<k4) z+=wt1[k3].wt;
	} else z+=wt[k3];
      } else z+=wt[k3];
    }
    *cost+=exp(z);
    delete_row_col(j);
    if(score[j]>=0.0) for(k1=1;k1<k;k1++) {
	k3=inv[k1];
	for(k2=0;k2<k1;k2++) {
	  k4=inv[k2];
	  add_node(k3,k4);
	  add_node(k4,k3);
	}
      }
    for(k1=0;k1<k;k1++) {
      k2=inv[k1];
      if(!flag[k2]) {
	degree[k2]=calc_score(k2,wt,wt1,inv1,inv2,flag,&j,score+k2,alist);
	for(k4=1;k4<j;k4++) {
	  k5=inv1[k4];
	  if(flag[k5]==3) {
	    flag[k5]=4;
	    n_nodes--;
	    delete_row_col(k5);
	  }
	}
	fill_in[k2]=check_fill_in(inv1,inv2,flag,j,alist);
	if(degree[k2]==1 || !fill_in[k2]) score[k2]=-1;
      }
    }
    n_nodes=0;
    for(i=0;i<n;i++) if(!flag[i]) perm[n_nodes++]=i;
    gnu_qsort(perm,(size_t)n_nodes,sizeof(int),cmp_func[fg]);
  }
  free(mat);
  free(alist);
  free(score);
  free(flag);
  return order;
}
