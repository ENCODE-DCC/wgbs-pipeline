/****************************************************************************
*                                                                          *
*     Loki - Programs for genetic analysis of complex traits using MCMC    *
*                                                                          *
*                     Simon Heath - CNG, Evry                              *
*                                                                          *
*                         February 2003                                    *
*                                                                          *
* oldpeel_compat.c:                                                        *
*                                                                          *
* Routines to convert from the new genotype elimination routines to the    *
* old style peeling routine.  Temporary (hopefully) until the peeling      *
* routines can be updated
*                                                                          *
* Copyright (C) Simon C. Heath 2003                                        *
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
#include <assert.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "loki_utils.h"
#include "lk_malloc.h"
#include "oldpeel_compat.h"

void alloc_reqset(struct Marker *mark,int ped_size)
{
	mark->req_set[0]=lk_malloc(sizeof(lk_ulong)*ped_size*2);
	mark->req_set[1]=mark->req_set[0]+ped_size;
}

void set_ind_allsets(const int i,const int ltype,struct Marker *mark,struct Id_Record *id_array)
{
	int j,k,k1,k2,k3,g[2],n_all,sex;
	gtype **gtypes;
	int *ngens,*nhaps[2],*prune,sire,dam,diploid=0,comp;
	lk_ulong *req_set[2],mask,mask1,a,a1,b,b1,c,c1,**all_set;
	
	prune=mark->locus.pruned_flag;
	if(prune[i]) return;
	comp=id_array[i].comp;
	sex=id_array[i].sex;
	gtypes=mark->gtypes;
	ngens=mark->ngens;
	all_set=mark->all_set;
	for(j=0;j<2;j++) {
		req_set[j]=mark->req_set[j];
		nhaps[j]=mark->nhaps[j];
	}
	n_all=mark->n_all1[comp];
	if(n_all>(int)(sizeof(lk_ulong)<<3)) ABT_FUNC("Too many alleles for current version\n");
	mask=(1UL<<n_all)-1UL;
	mask1=1UL<<(n_all-1);
	k=ngens[i];
	sex=id_array[i].sex;
	if(ltype==LINK_AUTO || sex==2) diploid=1;
	else if(sex==1) diploid=0;
	else ABT_FUNC("Unknown sex\n");
	if(diploid) {
		a=(mask^req_set[0][i])|mask1;
		a1=(mask^req_set[1][i])|mask1;
		if(k) {
			b=b1=0;
			for(k1=0;k1<k;k1++) {
				g[0]=gtypes[i][k1].mat;
				g[1]=gtypes[i][k1].pat;
				if(g[0] && (req_set[0][i]&(1UL<<(g[0]-1)))) g[0]=n_all;
				if(g[1] && (req_set[1][i]&(1UL<<(g[1]-1)))) g[1]=n_all;
				if(g[0]) {
					if(g[1]) {
						c=1UL<<(g[1]-1);
						all_set[i][g[0]-1]|=c;
						b|=c;
					} else {
						all_set[i][g[0]-1]|=a1;
						b|=a1;
					}
					b1|=(1<<(g[0]-1));
				} else {
					if(g[1]) c1=1UL<<(g[1]-1);
					else c1=mask1;
					c=a;
					j=0;
					while(c) {
						if(c&1) {
							all_set[i][j]|=c1;
							b|=c1;
							b1|=(1UL<<j);
						}
						j++;
						c>>=1;
					}
				}
			}
			k1=k2=k3=0;
			for(k=0;k<n_all;k++) {
				a=all_set[i][k];
				if(!a) continue;
				k1++;
				while(a)	{
					if(a&1) k3++;
					a>>=1;
				}
			}
			c=b;
			while(c)	{
				if(c&1) k2++;
				c>>=1;
			}
		} else {
			j=0;
			b1=a;
			b=a1;
			k1=0;
			while(a) {
				if(a&1) {
					all_set[i][j]=a1;
					k1++;
				}
				j++;
				a>>=1;
			}
			k2=0;
			while(a1) {
				if(a1&1) k2++;
				a1>>=1;
			}
			k3=k1*k2;
		}
		nhaps[X_MAT][i]=k1;
		nhaps[X_PAT][i]=k2;
		ngens[i]=k3;
		mark->temp[X_PAT][i]=b;
		mark->temp[X_MAT][i]=b1;
	} else {
		if(k) {
			b1=0;
			for(k1=0;k1<k;k1++) {
				g[0]=gtypes[i][k1].mat;
				assert(g[0]);
				all_set[i][g[0]-1]=1;
				b1|=(1UL<<(g[0]-1));
			}
		} else {
			a=(mask^req_set[0][i])|mask1;
			b1=a;
			j=0;
			k1=0;
			while(a) {
				if(a&1) {
					all_set[i][j]=1;
					k1++;
				}
				j++;
				a>>=1;
			}
		}
		nhaps[X_MAT][i]=k1;
		nhaps[X_PAT][i]=0;
		ngens[i]=k1;
		mark->temp[X_MAT][i]=b1;
		mark->temp[X_PAT][i]=0;
	}
	mark->m_flag[i]=0;
	if(diploid && ngens[i]==1) {
		for(k1=k=0;k<n_all;k++) {
			a=all_set[i][k];
			if(a) {
				if(a==(lk_ulong)(1<<k)) mark->m_flag[i]|=1;
				break;
			}
		}
	}
	sire=id_array[i].sire;
	if(sire && prune[sire-1]) sire=0;
	dam=id_array[i].dam;
	if(dam && prune[dam-1]) dam=0;
	if(dam) {
		b=req_set[X_MAT][i];
		b1=mark->temp[X_MAT][i];
		if(b&b1) b1|=b;
		b=mark->temp[X_PAT][dam-1];
		k1=0;
		if(b1&b) k1=1;
		b=mark->temp[X_MAT][dam-1];
		if(b1&b) k1|=2;
		assert(k1);
		if(k1!=3) {
			mark->m_flag[i]|=2;
			mark->locus.seg[X_MAT][i]=(k1==1?X_PAT:X_MAT);
		}
	}
	if(diploid && sire) {
		b=req_set[X_PAT][i];
		b1=mark->temp[X_PAT][i];
		if(b&b1) b1|=b;
		b=mark->temp[X_PAT][sire-1];
		k1=0;
		if(b1&b) k1=1;
		b=mark->temp[X_MAT][sire-1];
		if(b1&b) k1|=2;
		assert(k1);
		if(k1!=3) {
			mark->m_flag[i]|=4;
			mark->locus.seg[X_PAT][i]=(k1==1?X_PAT:X_MAT);
		}
	}
}

void setup_allset(struct Marker *mark,struct loki *loki)
{
	int i,i1,ped_size,n_all,ltype;
	int cs,comp;
	struct Id_Record *id_array;
	
	message(DEBUG_MSG,"Setting up allele sets for peeling for locus %s\n",mark->name);
	ltype=loki->markers->linkage[mark->locus.link_group].type&LINK_TYPES_MASK;
	assert(ltype==LINK_AUTO || ltype==LINK_X);
	ped_size=loki->pedigree->ped_size;
	id_array=loki->pedigree->id_array;
	n_all=mark->locus.n_alleles;
	mark->all_set=lk_malloc(sizeof(void *)*ped_size);
	mark->all_set[0]=lk_calloc((size_t)(ped_size*n_all),sizeof(lk_ulong));
	for(i=1;i<ped_size;i++) mark->all_set[i]=mark->all_set[i-1]+n_all;
	mark->temp[0]=lk_malloc(sizeof(lk_ulong)*2*ped_size);
	mark->temp[1]=mark->temp[0]+ped_size;
	mark->nhaps[0]=lk_calloc((size_t)(ped_size*2),sizeof(int));
	mark->nhaps[1]=mark->nhaps[0]+ped_size;
	i=0;
	for(comp=0;comp<loki->pedigree->n_comp;comp++) {
		cs=loki->pedigree->comp_size[comp];
		for(i1=0;i1<cs;i1++,i++) set_ind_allsets(i,ltype,mark,id_array);
	}
}
