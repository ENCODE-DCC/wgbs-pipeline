/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG                                    *
 *                                                                          *
 *                       April 2002                                         *
 *                                                                          *
 * count_relationships.c:                                                   *
 *                                                                          *
 * Counts different types of relationships                                  *
 *                                                                          *
 * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *
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
#include "y.tab.h"
#include "scan.h"
#include "count_relationships.h"

static struct relate_off **rel_mat;
static struct relate_atom *rel_buf;
static int rel_buf_size=16;

static double kin(int a,int b) 
{
	int i,j,ids,ids1;
	
	if(!(a&&b)) return 0.0;
	i=ped_recode1[a-1];
	j=ped_recode1[b-1];
	if(!(i&&j)) return 0.0;
	if(i==j) return .5*(1.0+kin(id_array[a-1].sire,id_array[a-1].dam));
	ids=id_array[a-1].sire;
	if(ids) ids=ped_recode1[ids-1];
	ids1=id_array[b-1].sire;
	if(ids1) ids1=ped_recode1[ids1-1];
	if(!ids) {
		if(j<i || !ids1) return 0.0;
		return .5*(kin(a,id_array[b-1].sire)+kin(a,id_array[b-1].dam));
	}
	if(!ids1) {
		if(i<j) return 0.0;
		return .5*(kin(b,id_array[a-1].sire)+kin(b,id_array[a-1].dam));
	}
	if(i<j) return .5*(kin(a,id_array[b-1].sire)+kin(a,id_array[b-1].dam));
	return .5*(kin(b,id_array[a-1].sire)+kin(b,id_array[a-1].dam));
}

void insert_rel(int i,struct relate_off *p)
{
	int j;
	struct relate_off *p1,**p2;
	
	j=p->x;
	p2=rel_mat+i-1;
	p1=*p2;
	while(p1) {
		if(p1->x>j) break;
		p2=&(p1->next);
		p1=*p2;
	}
	*p2=p;
	p->next=p1;
}

struct relate_off *find_rel(int a,int b)
{
	int i,j;
	struct relate_off *p;
	
	if(!(a&&b)) return 0;
	i=ped_recode1[a-1];
	j=ped_recode1[b-1];
	if(!(i&&j)) return 0;
	p=rel_mat[i-1];
	while(p) {
		if(p->x>=j) break;
		p=p->next;
	}
	if(p && p->x==j) return p;
	return 0;
}

double calc_rel(struct relate_off *p)
{
	struct relate_atom *p1;
	double z=0.0;
	int i,j;
	
	p1=p->atoms;
	for(i=0;i<p->n;i++,p1++) {
		j=p1->deg[0]+p1->deg[1];
		if(p1->type==REL_FULL) j--;
		z+=(double)p1->n/(double)(1<<j);
	}
	return z;
}

void print_off(struct relate_off *p)
{
	int i;
	struct relate_atom *p1;
	
	p1=p->atoms;
	for(i=0;i<p->n;i++) {
		printf("%d",p1->n);
		fputc(p1->type==REL_HALF?'H':'F',stdout);
		printf("%d,%d",p1->deg[0],p1->deg[1]);
		fputc(' ',stdout);
	}
	printf("%g\n",p->p);
}

struct relate_off *get_relate(int a,int b)
{
	int i,j,k=0,k1,k2,l,ids,idd,ids1,idd1;
	struct relate_off *o1=0,*o2;
	
	if(!(a&&b)) return 0;
	i=ped_recode1[a-1];
	j=ped_recode1[b-1];
	if(!(i&&j)) return 0;
	if(i==j) {
		o1=find_rel(a,id_array[a-1].sire);
		o2=find_rel(a,id_array[a-1].dam);
		l=1;
		if(o1) l+=o1->n;
		if(o2) l+=o2->n;
		if(l>rel_buf_size) {
			rel_buf_size=l;
			free(rel_buf);
			if(!(rel_buf=malloc(sizeof(struct relate_atom)*rel_buf_size))) ABT_FUNC(MMsg);
		}
		rel_buf[k].n=1;
		rel_buf[k].type=REL_HALF;
		rel_buf[k].deg[0]=0;
		rel_buf[k++].deg[1]=0;
		if(o1) for(k1=0;k1<o1->n;k1++,k++) memcpy(rel_buf+k,o1->atoms+k1,sizeof(struct relate_atom));
		if(o2) {
			for(k1=0;k1<o2->n;k1++) {
				for(k2=0;k2<k1;k2++) {
					if(rel_buf[k2].type==o2->atoms[k1].type && rel_buf[k2].deg[0]==o2->atoms[k1].deg[0] && 
						rel_buf[k2].deg[1]==o2->atoms[k1].deg[1]) {
						rel_buf[k2].n++;
						break;
					}
				}
				if(k2==k1) memcpy(rel_buf+(k++),o2->atoms+k1,sizeof(struct relate_atom));
			}
		}
		if(!(o1=malloc(sizeof(struct relate_off)))) ABT_FUNC(MMsg);
		o1->n=k;
		o1->x=i;
		if(!(o1->atoms=malloc(sizeof(struct relate_atom)*k))) ABT_FUNC(MMsg);
		for(k1=0;k1<k;k1++) memcpy(o1->atoms+k1,rel_buf+k1,sizeof(struct relate_atom));
		o1->p=calc_rel(o1);
		insert_rel(i,o1);
	} else {
		ids1=idd1=0;
		ids=id_array[b-1].sire;
		if(ids) ids1=ped_recode1[ids-1];
		idd=id_array[b-1].dam;
		if(idd) idd1=ped_recode1[idd-1];
		l=0;
		if(ids1==i) {
			o1=find_rel(ids,ids);
			if(!o1) o1=get_relate(ids,ids);
			l=o1->n;
		} else o1=0;
		if(idd1==i) {
			o2=find_rel(idd,idd);
			if(!o2) o2=get_relate(idd,idd);
			l+=o2->n;
		} else o2=0;
		if(!(o1 || o2)) {
			k1=0;
			if(ids1 && id_array[a-1].sire==ids) k1++;
			if(idd1 && id_array[a-1].dam==idd) k1++;
			if(k1) {
				print_orig_id(stdout,a,1);
				print_orig_id(stdout,b,1);
				printf("Ha!\n");
			}
		}
		if(l) {
			if(l>rel_buf_size) {
				rel_buf_size=l;
				free(rel_buf);
				if(!(rel_buf=malloc(sizeof(struct relate_atom)*rel_buf_size))) ABT_FUNC(MMsg);
			}
			if(o1) for(k1=0;k1<o1->n;k1++,k++) memcpy(rel_buf+k,o1->atoms+k1,sizeof(struct relate_atom));
			if(o2) {
				for(k1=0;k1<o2->n;k1++) {
					for(k2=0;k2<k1;k2++) {
						if(rel_buf[k2].type==o2->atoms[k1].type && rel_buf[k2].deg[0]==o2->atoms[k1].deg[0] && 
							rel_buf[k2].deg[1]==o2->atoms[k1].deg[1]) {
							rel_buf[k2].n++;
							break;
						}
					}
					if(k2==k1) memcpy(rel_buf+(k++),o2->atoms+k1,sizeof(struct relate_atom));
				}
			}
			for(k1=0;k1<k;k1++) {
				rel_buf[k1].deg[1]++;
			}
			if(!(o1=malloc(sizeof(struct relate_off)))) ABT_FUNC(MMsg);
			o1->n=k;
			o1->x=j;
			if(!(o1->atoms=malloc(sizeof(struct relate_atom)*k))) ABT_FUNC(MMsg);
			for(k1=0;k1<k;k1++) memcpy(o1->atoms+k1,rel_buf+k1,sizeof(struct relate_atom));
			o1->p=calc_rel(o1);
			insert_rel(i,o1);
		}
	}
	if(o1) {
		print_orig_id(stdout,a,1);
		print_orig_id(stdout,b,1);
		print_off(o1);
	}
	return o1;
}

void count_relationships(void)
{
	int i,j,j1,nr,i1,comp,cs;
	int *counts,*counts1,*perm;
	struct relate_off *off;
	double z;
	char *relate[]={"full-sibs","half-sibs","parent:offspring","grandparent:grandchild","great-grandparent:great-grandchild",
		  "avuncular","half-avuncular",
		  "half first-cousins","first-cousins","first-cousins","double first-cousins",0};

	nr=0;
	while(relate[nr]) nr++;
	if(!nr) return;
	if(!(counts=calloc((size_t)(2*nr+1),sizeof(int)))) ABT_FUNC(MMsg);
	counts1=counts+nr;
	if(!(perm=malloc(sizeof(int)*pruned_ped_size))) ABT_FUNC(MMsg);
	for(i=0;i<ped_size;i++)	{
		j=ped_recode1[i];
		if(j) perm[j-1]=i;
	}
	if(!(rel_mat=malloc(sizeof(void *)*pruned_ped_size))) ABT_FUNC(MMsg);
	if(!(rel_buf=malloc(sizeof(struct relate_atom)*rel_buf_size))) ABT_FUNC(MMsg);
	comp=id_array[perm[0]].component;
	cs=0;
	if(!comp_sflag || comp<n_comp) for(i1=1;i1<pruned_ped_size;i1++) {
		relate[i1]=0;
		i=perm[i1];
		if(id_array[i].component!=comp) {
			cs=i1;
			comp=id_array[i].component;
			if(comp_sflag && comp==n_comp) break;
			continue;
		}
		if(!id_array[i].sire || !ped_recode1[id_array[i].sire-1]) continue;
		for(j1=cs;j1<i1;j1++) {
			j=perm[j1];
			off=get_relate(j+1,i+1);
			if(off) {
				z=off->p;
			}
		}
	}
	free(perm);
	free(counts);
}
