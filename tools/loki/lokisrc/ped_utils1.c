/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                         February 2003                                    *
 *                                                                          *
 * ped_utils.c:                                                             *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
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
#include <string.h>
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "lk_malloc.h"
#include "bin_tree.h"
#include "loki_utils.h"
#include "ped_utils.h"

static struct lk_fam family;

#define HAS_DATA 1
#define IS_PRUNED 16

/* Do a locus by locus prune
 * We must remove:
 *  (a) Untyped individuals with no unpruned descendents
 *  (b) Untyped founder couples with only 1 unpruned child
 *  (c) Marriages with no unpruned children
 */
void prune_ped_for_locus(struct Locus *loc,struct loki *loki)
{
  int i,j,k,k1,ix,ped_size,n_gen,n_fam,n_comp,fg,tl_flag=0,*ff,*prune,**hap=0,n_loci=1,*hap1;
  struct Id_Record *id_array,*father,*mother,**kids,*kid;
  struct Marker *mark=0;
  nuc_fam *fam;
	
  ped_size=loki->pedigree->ped_size;
  k1=loki->models->n_models;
  ix=loc->index;
  id_array=loki->pedigree->id_array;
  n_gen=loki->pedigree->n_genetic_groups;
  n_comp=loki->pedigree->n_comp;
  if(loc->type&ST_MARKER) {
    mark=loki->markers->marker+ix;
    assert(!mark->parent);
    if(mark->n_children) n_loci=mark->n_children;
    if(n_loci>1) {
      hap=lk_malloc(sizeof(void *)*n_loci);
      for(i=0;i<n_loci;i++) hap[i]=mark->children[i]->haplo;
    } else {
      hap=&hap1;
      hap1=mark->haplo;
    }
  } else if(loc->eff || (loc->type&ST_TRAITLOCUS)) tl_flag=1;
  /* See who is oberved for this locus */
  for(i=0;i<ped_size;i++) {
    j=0;
    if(mark) {
      for(k=0;k<n_loci;k++) if(hap[k]) break;
      if(k<n_loci) j=HAS_DATA;
    }
    if(n_gen>1 && id_array[i].group) j=HAS_DATA;
    if(tl_flag) {
      for(k=0;k<k1;k++) if(id_array[i].res[k]) break;
      if(k<k1) j=HAS_DATA;
    }
    id_array[i].flag=j;
    id_array[i].n_fam=0;
  }
  fam=loki->family->families;
  n_fam=loki->family->cp_start[n_comp];
  for(i=0;i<n_fam;i++) {
    fam[i].flag=0;
    if(fam[i].father) {
      fam[i].father->n_fam++;
      fam[i].mother->n_fam++;
    }
  }
  do {
    for(fg=i=0;i<n_fam;i++) if(!(fam[i].flag)) {
      kids=fam[i].kids;
      k=0;
      while((kid=(*kids++))) {
	if(!(kid->flag&IS_PRUNED)) {
	  if((kid->flag&HAS_DATA) || kid->n_fam) k++;
	  else {
	    kid->flag|=IS_PRUNED;
	    fg=1;
	  }
	}
      }
      father=fam[i].father;
      mother=fam[i].mother;
      if(k==1 && father) { /* Only one non-prunable child in family */
	if(!(father->flag&HAS_DATA) && father->n_fam==1) { /* If father has no data and this is only family */
	  if(father->family) {
	    if(father->family->flag) k|=2;
	  } else k|=2;
	}
	if(!(mother->flag&HAS_DATA) && mother->n_fam==1) {
	  if(mother->family) {
	    if(mother->family->flag) k|=4;
	  } else k|=4;
	}
	if(k==7) k=0;  /* If both father & mother belong to no unpruned family, set k to 0 */
      }
      if(!k) { /* If family can be pruned... */
	fam[i].flag=fg=1;
	if(father) {
	  father->n_fam--;
	  mother->n_fam--;
	  if(!father->n_fam && !(father->flag&HAS_DATA)) father->flag|=IS_PRUNED;
	  if(!mother->n_fam && !(mother->flag&HAS_DATA)) mother->flag|=IS_PRUNED;
	}
      }
    }
  } while(fg);
  ff=loc->founder_flag;
  prune=loc->pruned_flag;
  for(i=0;i<ped_size;i++) {
    prune[i]=(id_array[i].flag&IS_PRUNED)?1:0;
    if(prune[i]) ff[i]=2;
    else {
      k=id_array[i].sire;
      if(k) {
	if(prune[k-1]) ff[i]=1;
	else ff[i]=0;
      } else ff[i]=1;
    }
  }
  if(n_loci>1) free(hap);
}

void prune_ped(struct loki *loki)
{
  int i;
  struct Marker *mark;
	
  message(INFO_MSG,"Locus specific pedigree pruning\n");
  mark=loki->markers->marker;
  for(i=0;i<loki->markers->n_markers;i++,mark++) if(!mark->parent) {
    prune_ped_for_locus(&mark->locus,loki);
  }
  if(loki->models->tlocus) prune_ped_for_locus(loki->models->tlocus,loki);
}

static void free_families(void)
{
  if(family.cp_start) free(family.cp_start);
  if(family.families) {
    if(family.families[0].kids) free(family.families[0].kids);
    free(family.families);
  }
}

/* Build up nuclear family structures */
void build_families(struct loki *loki)
{
  int i,j,k,n_kids=0,cs,ped_size,comp,n_comp,*cst,n_fam=0;
  int ids,idd;
  struct Id_Record *id_array,**kid_list,**kd,*kid;
  nuc_fam *fam;
	
  loki->family=&family;
  n_comp=loki->pedigree->n_comp;
  cst=lk_malloc(sizeof(int)*(n_comp+1));
  if(atexit(free_families)) message(WARN_MSG,"Unable to register exit function free_families()\n");
  id_array=loki->pedigree->id_array;
  family.cp_start=cst;
  ped_size=loki->pedigree->ped_size;
  n_kids=0;
  for(i=0;i<ped_size;i++) {
    id_array[i].flag=id_array[i].nkids=0;
    ids=id_array[i].sire;
    idd=id_array[i].dam;
    if(ids) id_array[ids-1].nkids++;
    if(idd) id_array[idd-1].nkids++;
    n_kids+=2;
  }
  if(!n_kids) return;
  kd=lk_malloc(sizeof(void *)*n_kids);
  loki->sys.RemBlock=AddRemem(kd,loki->sys.RemBlock);
  for(i=0;i<ped_size;i++) {
    if(id_array[i].nkids) {
      id_array[i].kids=kd;
      kd+=id_array[i].nkids;
      id_array[i].nkids=0;
    }
    id_array[i].family=0;
    ids=id_array[i].sire;
    idd=id_array[i].dam;
    if(ids) id_array[ids-1].kids[id_array[ids-1].nkids++]=id_array+i;
    if(idd) id_array[idd-1].kids[id_array[idd-1].nkids++]=id_array+i;
  }
  n_kids=0;
  for(i=comp=0;comp<n_comp;comp++) {
    cs=loki->pedigree->comp_size[comp];
    cst[comp]=n_fam;
    if(loki->pedigree->singleton_flag[comp]) { /* Components consisting of singletons we'll lump into a 'pseudo-family' */
      i+=cs;
      n_kids+=cs+1;
      n_fam++;
      continue;
    }
    for(j=0;j<cs;j++,i++) {
      if(id_array[i].sire && !id_array[i].flag) {
	ids=id_array[i].sire;
	idd=id_array[i].dam;
	for(k=0;k<id_array[ids-1].nkids;k++) {
	  kid=id_array[ids-1].kids[k];
	  if(kid->dam==idd) {
	    kid->flag=1;
	    n_kids++;
	  }
	}
	n_kids++;
	n_fam++;
      }
    }
  }
  cst[comp]=n_fam;
  if(!n_fam) return;
  fam=lk_malloc(sizeof(nuc_fam)*n_fam);
  kid_list=lk_malloc(sizeof(void *)*n_kids);
  family.families=fam;
  fam->kids=0;
  n_kids=0;
  for(i=comp=0;comp<n_comp;comp++) {
    cs=loki->pedigree->comp_size[comp];
    if(loki->pedigree->singleton_flag[comp]) {
      fam->father=fam->mother=0;
      fam->kids=kid_list+n_kids;
      fam->comp=comp;
      fam->gtypes=0;
      fam->n_err=0;
      for(j=0;j<cs;j++,i++) {
	kid_list[n_kids++]=id_array+i;
	id_array[i].family=fam;
      }
      kid_list[n_kids++]=0;
      fam++;
      continue;
    }
    for(j=0;j<cs;j++,i++) {
      if(id_array[i].sire && id_array[i].flag==1) {
	ids=id_array[i].sire;
	idd=id_array[i].dam;
	fam->father=id_array+ids-1;
	fam->mother=id_array+idd-1;
	fam->kids=kid_list+n_kids;
	fam->comp=comp;
	fam->gtypes=0;
	fam->n_err=0;
	for(k=0;k<id_array[ids-1].nkids;k++) {
	  kid=id_array[ids-1].kids[k];
	  if(kid->dam==idd) {
	    kid->flag=2;
	    kid->family=fam;
	    kid_list[n_kids++]=kid;
	  }
	}
	kid_list[n_kids++]=0;
	fam++;
      }
    }
  }
}

int check_ped(struct loki *loki)
{
  int i,j=0,k,k1,k2,er=0,ped_size,n_grp;
  struct Id_Record *id_array;
	
  message(DEBUG_MSG,"Checking sex & group information in pedigree\n");
  ped_size=loki->pedigree->ped_size;
  id_array=loki->pedigree->id_array;
  n_grp=loki->pedigree->n_genetic_groups;
  /* First check if we *need* sex information (if we have a sex specific map, 
     or if there is a sexlinked chromosome */
  if(loki->markers->sex_map) j=1;
  else {
    for(k=0;k<loki->markers->n_links;k++) {
      if((loki->markers->linkage[k].type&LINK_TYPES)!=LINK_AUTO) j=1;
    }
  }
  if(j) message(DEBUG_MSG,"Sex information required\n");
  /* Check sex information in the pedigree, filling in missing values for parents */
  for(i=0;i<ped_size;i++) id_array[i].flag=0;
  for(i=ped_size-1;i>=0;i--) {
    k=id_array[i].sire;
    if(k--) {
      if(id_array[k].sex) {
	if(id_array[k].sex!=1 && !id_array[k].flag) {
	  id_array[k].flag=1;
	  print_orig_id(stderr,k+1);
	  fputs(" has inconsistent sex information\n",stderr);
	  er|=1;
	}
      } else id_array[k].sex=1;
    }
    k=id_array[i].dam;
    if(k--) {
      if(id_array[k].sex) {
	if(id_array[k].sex!=2 && !id_array[k].flag) {
	  id_array[k].flag=1;
	  print_orig_id(stderr,k+1);
	  fputs(" has inconsistent sex information\n",stderr);
	  er|=1;
	}
      } else id_array[k].sex=2;
    }
    if(j && !id_array[i].sex) {
      print_orig_id(stderr,i+1);
      fputs(" has missing sex information\n",stderr);
      er|=2;
    }
  }
  for(i=0;i<ped_size;i++) {
    if(n_grp>1 && !id_array[i].group) {
      k1=k2=0;
      if((k=id_array[i].sire)) k1=id_array[k-1].group;
      if((k=id_array[i].dam)) k2=id_array[k-1].group;
      if(k1!=k2 || !k1) {
	print_orig_id(stderr,i+1);
	fputs(" has missing group information\n",stderr);
	er|=1;
      }
      id_array[i].group=k1;
    }
    if(!id_array[i].group) {
      if(n_grp>1) {
	print_orig_id(stderr,i+1);
	fputs(" has missing group information\n",stderr);
	er|=1;
      } else id_array[i].group=1;
    }
  }
  return er;
}

struct loop_link {
  int *ids;
  int n_ids;
  int idx;
};

struct edge {
  struct edge *next;
  unsigned long *link;
};

struct loop_connect {
  struct loop_connect *next;
  struct loop_node *node;
  struct edge *edge;
};

struct loop_node {
  nuc_fam *fam;
  struct loop_connect *conns;
  int order;
};

static struct loop_node *nodes;
static struct loop_link *links;
static int n_nodes,n_links,n_longs,n_bits;
static struct edge *free_edge_list;
static struct loop_connect *free_connect_list;
static unsigned long *tedge;
static struct bin_node *free_bnode_list,*tmp_edge_root;

static struct bin_node *alloc_tmp_edge(const void *s)
{
  struct bin_node *nd;

  if((nd=free_bnode_list)) free_bnode_list=nd->left;
  else {
    nd=lk_malloc(sizeof(struct bin_node));
    nd->data=lk_malloc(sizeof(long)*n_longs);
  }
  nd->left=nd->right=0;
  nd->balance=0;
  memcpy(nd->data,s,sizeof(long)*n_longs);
  return nd;
}

static void free_tmp_edge(struct bin_node *nd)
{
  nd->left=free_bnode_list;
  free_bnode_list=nd;
}

static int cmp_edge(const void *s1,const void *s2)
{
  return memcmp(s1,s2,sizeof(long)*n_longs);
}

static struct edge *new_edge(int nl)
{
  struct edge *e;
  
  assert(nl);
  if((e=free_edge_list)) {
    free_edge_list=e->next;
  } else {
    e=lk_malloc(sizeof(struct edge));
    e->link=lk_malloc(sizeof(long)*nl);
  }
  e->next=0;
  return e;
}

static struct loop_connect *new_connect(void)
{
  struct loop_connect *c;

  if((c=free_connect_list)) {
    free_connect_list=c->next;
  } else c=lk_malloc(sizeof(struct loop_connect));
  c->next=0;
  return c;
}

static void free_connects(struct loop_node *n1)
{
  struct loop_node *n2;
  struct loop_connect *c,*c1,**cp;
  struct edge *e;

  c=n1->conns;
  while(c) {
    e=c->edge;
    assert(e);
    while(e->next) e=e->next;
    e->next=free_edge_list;
    free_edge_list=c->edge;
    n2=c->node;
    cp=&n2->conns;
    c1=0;
    while((c1=*cp)) {
      if(c1->node==n1) {
	*cp=c1->next;
	c1->next=free_connect_list;
	free_connect_list=c1;
	break;
      }
      cp=&c1->next;
    }
    n2->order--;
    assert(c1);
    c=c->next;
  }
  n1->conns=0;
}

static void output_edge(FILE *fptr,int nl,unsigned long *t)
{
  int i,i1,i2,j,k;
  unsigned long a;

  for(i=i1=i2=0;i<nl;i++) {
    a=t[i];
    j=0;
    while(a) {
      if(a&1) {
	if(i2++) fputc(',',fptr);
	if(links[i1+j].n_ids>1) fputc('[',stdout);
	for(k=0;k<links[i1+j].n_ids;k++) {
	  if(k) fputc(',',fptr);
	  print_orig_id(fptr,links[i1+j].ids[k]+1);
	}
	if(k>1) fputc(']',stdout);
      }
      a>>=1;
      j++;
    }
    i1+=n_bits;
  }
}

static struct loop_connect **get_connects(struct loop_node *n1,struct loop_node *n2) 
{
  struct loop_connect *cn;
  static struct loop_connect *cnp[2];

  cn=n1->conns;
  while(cn) {
    if(cn->node==n2) break;
    cn=cn->next;
  }
  if(cn) {
    cnp[0]=cn;
    cn=n2->conns;
    while(cn) {
      if(cn->node==n1) break;
      cn=cn->next;
    }
    assert(cn);
    cnp[1]=cn;
  } else {
    cn=new_connect();
    cn->node=n2;
    cn->edge=0;
    cn->next=n1->conns;
    n1->conns=cn;
    cnp[0]=cn;
    /* and from n2 */
    cn=new_connect();
    cn->node=n1;
    cn->edge=0;
    cn->next=n2->conns;
    n2->conns=cn;
    cnp[1]=cn;
    n1->order++;
    n2->order++;
  }
  return cnp;
}

static void add_edges(struct bin_node *nd,struct edge **ep)
{
  int i=1;
  struct edge *e,*e1;

  if(nd->left) add_edges(nd->left,ep);
  while((e=*ep)) {
    i=memcmp(e->link,nd->data,sizeof(long)*n_longs);
    if(i<=0) break;
    ep=&e->next;
  }
  if(i) { /* New node required */
    e1=new_edge(n_longs);
    memcpy(e1->link,nd->data,sizeof(long)*n_longs);
    e1->next=e;
    *ep=e1;
  }
  free_tmp_edge(nd);
  if(nd->right) add_edges(nd->right,ep);
}

#if 0
static struct loop_connect *add_edge(struct loop_connect *cn,struct loop_node *n1,struct loop_node *n2,int nl,unsigned long *ed) 
{
  int i;
  struct edge *e,*e1;
  unsigned long *lp1,*lp2;
 
  if(cn) {
    /* Check if edge already exists (if so, ignore) */
    e=cn->edge;
    while(e) {
      lp1=e->link;
      lp2=ed;
      for(i=0;i<nl;i++) if(*lp1++!=*lp2++) break;
      if(i==nl) break;
      e=e->next;
    }
    if(!e) { /* New edge for existing connection */
      e=cn->edge;
      assert(e);
      while(e->next) e=e->next;
      e1=new_edge(nl);
      memcpy(e1->link,ed,sizeof(long)*nl);
      e->next=e1;
    }
  } else {
    /* Create new edge */
    e=new_edge(nl);
    memcpy(e->link,ed,sizeof(long)*nl);
    /* Set up connection from n1 */
    cn=new_connect();
    cn->node=n2;
    cn->edge=e;
    cn->next=n1->conns;
    n1->conns=cn;
    /* and from n2 */
    cn=new_connect();
    cn->node=n1;
    cn->edge=e;
    cn->next=n2->conns;
    n2->conns=cn;
    n1->order++;
    n2->order++;
  }
  return cn;
}

#endif

static void add_connects(struct loop_node *n1,struct loop_node *n2,int nl,unsigned long *ed) 
{
  int i;
  struct loop_connect *cn;
  struct edge *e,*e1;
  unsigned long *lp1,*lp2;

  cn=n1->conns;
  /*  fputs("add_connects: (",stdout);
  print_orig_id(stdout,n1->fam->father->idx+1);
  fputc(',',stdout);
  print_orig_id(stdout,n1->fam->mother->idx+1);
  printf(") <-> (");
  print_orig_id(stdout,n2->fam->father->idx+1);
  fputc(',',stdout);
  print_orig_id(stdout,n2->fam->mother->idx+1);
  printf(") ");
  output_edge(stdout,nl,ed);
  fputc('\n',stdout); */
  while(cn) {
    if(cn->node==n2) break;
    cn=cn->next;
  }
  if(cn) { /* Connection exists */
    /* Check if edge already exists (if so, ignore) */
    e=cn->edge;
    while(e) {
      lp1=e->link;
      lp2=ed;
      #if 0
      for(i=j=0;i<nl;i++) {
	a=*lp1++;
	b=*lp2++;
	if(a!=b) {
	  c=a&b;
	  if(c==a) j|=1;
	  else if(c==b) j|=2;
	  else j=3;
	  if(j==3) break;
	}
      }
      /*      if(j<3) {
	printf("AK: ");
	output_edge(stdout,nl,e->link);
	printf("  ");
	output_edge(stdout,nl,ed);
	printf(" j=%d\n",j);
	} */
      if(j<3) break;
      #endif
      for(i=0;i<nl;i++) if(*lp1++!=*lp2++) break;
      if(i==nl) break;
      e=e->next;
    }
    if(!e) { /* New edge for existing connection */
      e=cn->edge;
      assert(e);
      while(e->next) e=e->next;
      e1=new_edge(nl);
      memcpy(e1->link,ed,sizeof(long)*nl);
      e->next=e1;
    } else {
      printf("OOOK!\n");
    }
#if 0
    else if(j==2) {
      memcpy(e->link,ed,sizeof(long)*nl);
    }
#endif
  } else { /* New connection */
    /* Create new edge */
    e=new_edge(nl);
    memcpy(e->link,ed,sizeof(long)*nl);
    /* Set up connection from n1 */
    cn=new_connect();
    cn->node=n2;
    cn->edge=e;
    cn->next=n1->conns;
    n1->conns=cn;
    /* and from n2 */
    cn=new_connect();
    cn->node=n1;
    cn->edge=e;
    cn->next=n2->conns;
    n2->conns=cn;
    n1->order++;
    n2->order++;
  }
}

static int add_to_flist(nuc_fam *s,int i,nuc_fam **flist,int *sz)
{
  int j;

  for(j=0;j<i;j++) if(flist[j]==s) break;
  if(j==i) {
    if(i==*sz) {
      (*sz)<<=1;
      flist=lk_realloc(flist,sizeof(int)*(*sz));
    }
    flist[i++]=s;
  }
  return i;
}

static void get_links(nuc_fam *fam,struct loop_node *node,nuc_fam *fam1,int *path,int path_len) 
{
  int i,j,k,size,*tpath=0;
  struct loop_node *nd=0;
  struct loop_link *link;
  struct loop_connect **cnp;
  struct Id_Record *par,*kid,*kid1;
  nuc_fam **flist;

  if(fam->flag>2) {
    nd=nodes+(n_nodes++);
    nd->fam=fam;
    nd->conns=0;
    nd->order=0;
  }
  if(node) {
    if(abs(fam->flag)>2) {
      /* Check for duplicates - at this stage this can
	 only occur in the case of multiple marriages, and
	 will only involve a single individual */
      link=0;
      if(path_len==1) {
	for(i=0;i<n_links;i++) if(links[i].n_ids==1 && links[i].ids[0]==path[0]) {
	  link=links+i;
	  break;
	}
      }
      /* Create new link */
      if(!link) {
	link=links+(n_links++);
	link->ids=lk_malloc(sizeof(int)*path_len);
	link->n_ids=path_len;
	link->idx=n_links-1;
	memcpy(link->ids,path,sizeof(int)*path_len);
	printf("Adding path %d:",link->idx);
	for(i=0;i<path_len;i++) {
	  fputc(' ',stdout);
	  print_orig_id(stdout,path[i]+1);
	}
	fputc('\n',stdout);
      }
      path_len=0;
      if(!nd) {
	for(i=0;i<n_nodes;i++) if(nodes[i].fam==fam) {
	  nd=nodes+i;
	  break;
	}
	assert(nd);
      }
      i=link->idx;
      j=i%n_bits;
      i/=n_bits;
      memset(tedge,0,sizeof(long)*n_longs);
      tedge[i]=1UL<<j;
      tmp_edge_root=alloc_tmp_edge(tedge);
      cnp=get_connects(node,nd);
      add_edges(tmp_edge_root,&cnp[0]->edge);
      cnp[1]->edge=cnp[0]->edge;
      /*      add_connects(node,nd,n_longs,tedge); */
      node=nd;
    }
  } else node=nd;
  assert(node);
  if(fam->flag>0) {
    fam->flag*=-1;
    size=8;
    flist=lk_malloc(sizeof(void *)*size);
    if(tpath) {
      tpath=lk_malloc(sizeof(int)*path_len);
      memcpy(tpath,path,sizeof(int)*path_len);
    }
    par=fam->father;
    i=0;
    if(par->family && par->family!=fam1 && par->family!=fam && par->family->flag!=1) i=add_to_flist(par->family,0,flist,&size);
    if(par->n_fam>1) {
      for(j=0;j<par->nkids;j++) {
	kid=par->kids[j];
	if(kid->family!=fam && kid->family!=fam1 && kid->family->flag!=1) i=add_to_flist(kid->family,i,flist,&size);
      }
    }
    if(i) {
      for(j=0;j<i;j++) {
	if(tpath) memcpy(path,tpath,sizeof(int)*path_len);
	path[path_len]=par->idx;
	get_links(flist[j],node,fam,path,path_len+1);
      }
    }
    par=fam->mother;
    i=0;
    if(par->family && par->family!=fam1 && par->family!=fam && par->family->flag!=1) i=add_to_flist(par->family,0,flist,&size);
    if(par->n_fam>1) {
      for(j=0;j<par->nkids;j++) {
	kid=par->kids[j];
	if(kid->family!=fam && kid->family!=fam1 && kid->family->flag!=1) i=add_to_flist(kid->family,i,flist,&size);
      }
    }
    if(i) {
      for(j=0;j<i;j++) {
	if(tpath) memcpy(path,tpath,sizeof(int)*path_len);
	path[path_len]=par->idx;
	get_links(flist[j],node,fam,path,path_len+1);
      }
    }
    j=0;
    kid=fam->kids[j++];
    while(kid) {
      if(kid->n_fam) {
	i=0;
	for(k=0;k<kid->nkids;k++) {
	  kid1=kid->kids[k];
	  if(kid1->family!=fam1 && kid1->family->flag!=1) i=add_to_flist(kid1->family,i,flist,&size);
	}
	if(i) {
	  for(k=0;k<i;k++) {
	    if(tpath) memcpy(path,tpath,sizeof(int)*path_len);
	    path[path_len]=kid->idx;
	    get_links(flist[k],node,fam,path,path_len+1);
	  }
	}
      }
      kid=fam->kids[j++];
    }
    free(flist);
    if(tpath) free(tpath);
  }
}

void get_loops(struct loki *loki)
{
  int i,j,k,kk,comp,*cst,n_comp,n_fam,*path,n,n_loops;
  nuc_fam *fam;
  struct Id_Record *par,*kid,*pivot;
  struct loop_node *nd,**nd_list;
  struct loop_connect *c,*c1,**cnp;
  struct edge *e,*e1;
  struct bin_node *bin_nd;

  cst=loki->family->cp_start;
  fam=loki->family->families;
  n_comp=loki->pedigree->n_comp;
  n_fam=cst[n_comp];
  for(i=0;i<loki->pedigree->ped_size;i++) loki->pedigree->id_array[i].n_fam=0;
  for(i=0;i<n_fam;i++) {
    fam[i].flag=-1;
    if(fam[i].father) {
      fam[i].father->n_fam++;
      fam[i].mother->n_fam++;
    }
  }
  /* Get order (no. connections) for each nuclear family */
  for(comp=0;comp<n_comp;comp++) {
    if(loki->pedigree->singleton_flag[comp]) continue; /* Skip over singleton component */
    for(i=cst[comp];i<cst[comp+1];i++) {
      if(fam[i].flag>=0 && fam[i].flag<2) continue;
      j=0;
      pivot=0;
      par=fam[i].father;
      assert(par);
      /* Add connections from father */
      k=par->n_fam-1+((par->family && (par->family->flag<0 || par->family->flag>1))?1:0); 
      if(k) {
	pivot=par;
	j++;
      }
      par=fam[i].mother;
      assert(par);
      /* and connections from mother */
      k=par->n_fam-1+((par->family && (par->family->flag<0 || par->family->flag>1))?1:0); 
      if(k) {
	pivot=par;
	j++;
      }
      k=0;
      kid=fam[i].kids[k++];
      while(kid) {
	if(kid->n_fam) j++;
	kid=fam[i].kids[k++];
      }
      fam[i].flag=j;
      if(j==1) {
	if(pivot) pivot->n_fam--;
	i=cst[comp]-1;
      }
    }
    n_nodes=n_links=0;
    for(i=cst[comp];i<cst[comp+1];i++) {
      if(fam[i].flag>=0 && fam[i].flag<2) continue;
      j=0;
      pivot=0;
      par=fam[i].father;
      assert(par);
      /* Add connections from father */
      k=par->n_fam-1+((par->family && (par->family->flag<0 || par->family->flag>1))?1:0); 
      if(k) {
	pivot=par;
	j+=k;
      }
      par=fam[i].mother;
      assert(par);
      /* and connections from mother */
      k=par->n_fam-1+((par->family && (par->family->flag<0 || par->family->flag>1))?1:0); 
      if(k) {
	pivot=par;
	j+=k;
      }
      k=0;
      kid=fam[i].kids[k++];
      while(kid) {
	j+=kid->n_fam;
	kid=fam[i].kids[k++];
      }
      fam[i].flag=j;
      if(j>2) {
	n_nodes++;
	n_links+=j;
      }
      if(j>1) {
	printf("Node: ");
	print_orig_id(stdout,fam[i].father->idx+1);
	fputc(' ',stdout);
	print_orig_id(stdout,fam[i].mother->idx+1);
	printf(" %d\n",j);
      }
    }
    assert(n_links && !(n_links&1));
    n_links>>=1;
    printf("n_nodes=%d, n_links=%d\n",n_nodes,n_links);
    links=lk_malloc(sizeof(struct loop_link)*n_links);
    nodes=lk_malloc(sizeof(struct loop_node)*n_nodes);
    path=lk_malloc(sizeof(int)*loki->pedigree->comp_size[comp]);
    n_bits=sizeof(long)*8;
    n_longs=n_links/n_bits;
    if(n_links%n_bits) n_longs++;
    n_nodes=n_links=0;
    tedge=lk_malloc(sizeof(long)*n_longs);
    for(i=cst[comp];i<cst[comp+1];i++) {
      if(fam[i].flag>2) {
	get_links(fam+i,0,0,path,0);
      }
    }
    free(path);
    nd_list=lk_malloc(sizeof(void *)*n_nodes);
    for(i=0;i<n_nodes;i++) {
      nd_list[i]=nodes+i;
      c=nodes[i].conns;
      n=0;
      while(c) {
	fputs("NODE: (",stdout);
	print_orig_id(stdout,nodes[i].fam->father->idx+1);
	fputc(',',stdout);
	print_orig_id(stdout,nodes[i].fam->mother->idx+1);
	printf(") -> ");
	nd=c->node;
	e=c->edge;
	kk=0;
	while(e) {
	  if(kk) fputs(" |",stdout);
	  output_edge(stdout,n_longs,e->link);
	  kk++;
	  e=e->next;
	}
	printf(" -> (");
	print_orig_id(stdout,nd->fam->father->idx+1);
	fputc(',',stdout);
	print_orig_id(stdout,nd->fam->mother->idx+1);
	printf(")\n");
	n++;
	c=c->next;
      }
      printf("Order %d\n",n);
    }
    n_loops=0;
    while(n_nodes) {
      /* Pick minimum order node */
      j=nd_list[0]->order;
      k=0;
      for(i=1;i<n_nodes;i++) {
	if(nd_list[i]->order<j) {
	  j=nd_list[i]->order;
	  k=i;
	}
      }
      printf("n_nodes=%d, n_loops=%d\n",n_nodes,n_loops);
      fputs("Peeling: (",stdout);
      print_orig_id(stdout,nd_list[k]->fam->father->idx+1);
      fputc(',',stdout);
      print_orig_id(stdout,nd_list[k]->fam->mother->idx+1);
      printf(") %d\n",nd_list[k]->order);
      printf("Check 1\n");
      /* Check if any loops between this node and it's neighbours */
      c=nd_list[k]->conns;
      while(c) {
	e=c->edge;
	while(e) {
	  e1=e->next;
	  while(e1) {
	    for(i=0;i<n_longs;i++) if(e1->link[i]&e->link[i]) break;
	    if(i==n_longs) { /* Valid loop found - output */
	      n_loops++;
	      /*	      printf("Loop: ");
	      for(i=0;i<n_longs;i++) tedge[i]=e1->link[i]|e->link[i];
	      output_edge(stdout,n_longs,tedge);
	      fputc('\n',stdout);  */
	    }
	    e1=e1->next;
	  }
	  e=e->next;
	}
	c=c->next;
	printf("n_loops=%d\n",n_loops);
      }
      printf("Check 2\n");
      /* Create valid paths between neighbours */
      c=nd_list[k]->conns;
      j=0;
      while(c) {
	printf("At: %d/%d\n",++j,nd_list[k]->order);
	c1=c->next;
	while(c1) {
	  fputc('.',stdout);
	  fflush(stdout);
	  e=c->edge;
	  tmp_edge_root=0;
	  while(e) {
	    e1=c1->edge;
	    while(e1) {
	      for(i=0;i<n_longs;i++) if(e1->link[i]&e->link[i]) break;
	      if(i==n_longs) { /* Valid path found */
		for(i=0;i<n_longs;i++) tedge[i]=e1->link[i]|e->link[i];
		if(tmp_edge_root) tmp_edge_root=insert_bin_node(tmp_edge_root,tedge,&bin_nd,&kk,cmp_edge,alloc_tmp_edge);
		else tmp_edge_root=alloc_tmp_edge(tedge);
		/*		fputs("NE: (",stdout);
		print_orig_id(stdout,c->node->fam->father->idx+1);
		fputc(',',stdout);
		print_orig_id(stdout,c->node->fam->mother->idx+1);
		printf(") , (");
		print_orig_id(stdout,c1->node->fam->father->idx+1);
		fputc(',',stdout);
		print_orig_id(stdout,c1->node->fam->mother->idx+1);
		printf(") -> ");
		output_edge(stdout,n_longs,tedge);
		fputc('\n',stdout);  */
	      }
	      e1=e1->next;
	    }
	    e=e->next;
	  }
	  if(tmp_edge_root) {
	    cnp=get_connects(c->node,c1->node);
	    add_edges(tmp_edge_root,&cnp[0]->edge);
	    cnp[1]->edge=cnp[0]->edge;
	  }
	  c1=c1->next;
	}
	fputc('\n',stdout);
	c=c->next;
      }
      /* Remove old connections */
      free_connects(nd_list[k]);
      nd_list[k]=nd_list[--n_nodes];
    }
  }
  exit(0);
}
